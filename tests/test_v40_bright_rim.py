"""V40：亮圈独立定位、原图内沿尺寸、混合/缺失阵列及真实 worker。"""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'SEM_0830_0831_compare_V40_bundle'))
import sem_before_after_compare_0831_V40 as sem
from test_v40_recovery import process
from test_v38_slit_fusion import slit_scene
from test_v38_morphology import contour as irregular_contour
from test_v32 import analytic_trenches

np, cv = sem.np, sem.cv2


def rim_scene(pattern='slot', thickness=3, amplitude=70., contrast=0., texture=4., noise=.3,
              missing=(), conventional=(), gap=0, shape='normal', only_sites=None):
    """Interior and background coincide by default; the outer annulus is bright."""
    image_shape = (400,360) if pattern=='slot' else (360,430)
    yy, xx = np.indices(image_shape)
    background = 150.+.03*xx+.02*yy+texture*np.sin(xx/29.)*np.cos(yy/37.)
    total = np.zeros(image_shape,np.float32)
    ys = (55,145,235,325) if pattern=='slot' else (65,140,215,290)
    xs = (60,180,300) if pattern=='slot' else (65,140,215,290,365)
    for r,y in enumerate(ys):
        for c,x in enumerate(xs):
            if (r,c) in missing or (only_sites is not None and (r,c) not in only_sites):
                continue
            mask = np.zeros(image_shape,np.uint8)
            if pattern=='slot':
                hh = (55,43,31,21)[r]//2
                cv.rectangle(mask,(x-8,y-hh),(x+8,y+hh),1,-1)
            elif shape=='ellipse':
                cv.ellipse(mask,(x,y),(28,10),30.,0,360,1,-1)
            elif shape=='star':
                points = irregular_contour('star').reshape(-1,2)-[60,60]+[x,y]
                cv.fillPoly(mask,[np.rint(points).astype(np.int32)],1)
            else:
                cv.circle(mask,(x,y),20,1,-1)
            if (r,c) in conventional:
                total -= 80.*mask
                continue
            outer = cv.dilate(mask,cv.getStructuringElement(cv.MORPH_ELLIPSE,
                                                          (2*thickness+1,2*thickness+1)))
            ring = outer-mask
            if gap:
                edge_x = x+(8 if pattern=='slot' else 20)
                ring[y-gap//2:y+(gap+1)//2,edge_x:]=0
            total += amplitude*ring-contrast*mask
    raw = background+cv.GaussianBlur(total,(0,0),.8)
    raw += np.random.default_rng(40).normal(0,noise,image_shape)
    return np.clip(raw,0,255).astype(np.uint8)


def run(raw, pattern='slot', scale=1., pool=None, **options):
    if pool is None:
        pool,_ = sem.extract_descriptors_with_rejections(sem.segment_features(raw))
    key = 'slot210' if pattern=='slot' else 'via60'
    func = sem._select_and_measure_bottom_slots_v19 if pattern=='slot' else sem._select_and_measure_vias_v19
    with sem._v33_detection_context(options), contextlib.redirect_stdout(io.StringIO()):
        return func(raw,pool,scale,'image.png','single','007',key)


class BrightRimTests(unittest.TestCase):
    def assert_rim_rows(self, rows, count):
        self.assertEqual(len(rows),count)
        self.assertTrue(all(r['v40_bright_rim_used'] for r in rows))
        self.assertTrue(all(r['edge_method']=='V40_bright_rim_inner_gray50_normal_tracking' for r in rows))
        self.assertTrue(all(r['edge_signal_source']=='original_image_gray' and not r['edge_template_forced'] for r in rows))
        self.assertTrue(all(r['edge_support_fraction'] >= .85 and r['edge_max_missing_fraction'] <= .08 for r in rows))

    def test_equal_interior_and_background_is_located_from_all_bright_rims(self):
        for pattern, count in (('slot',12),('via',20)):
            with self.subTest(pattern=pattern):
                rims, _, meta = sem._v40_bright_rim_candidates(rim_scene(pattern),pattern,1.)
                self.assertEqual(len(rims),count)
                self.assertEqual(meta['v40_rim_candidate_count'],count)
                self.assertTrue(all(d._v32_edge_audit['v40_rim_prominence_gray']>10. for d in rims))

    def test_slit_rims_work_without_any_intensity_locator_and_keep_bottom_rule(self):
        raw = rim_scene()
        with patch.object(sem,'_v32_slot_flatfield_candidates',return_value=[]), \
             patch.object(sem,'_v37_feature_candidates',return_value=[]):
            rows,selected,bottom,_,meta = run(raw,pool=[])
        self.assert_rim_rows(rows,3)
        self.assertEqual((len(selected),len(bottom)),(12,3))
        self.assertTrue(all(r['slit_locator_path']=='V40_bright_rim' for r in rows))
        np.testing.assert_allclose([[r['x_width_nm'],r['y_height_nm']] for r in rows],[[17.,21.]]*3,atol=.5)
        self.assertEqual(meta['v40_slit_rim_observed_count'],12)

    def test_weak_one_pixel_rim_recovers_three_bottom_slits_after_v39_strict_zero(self):
        import importlib
        sys.path.insert(0,str(ROOT/'SEM_0830_0831_compare_V39_bundle'))
        old=importlib.import_module('sem_before_after_compare_0831_V39')
        yy,xx=np.indices((400,360))
        background=150.+.03*xx+.02*yy+20.*np.sin(xx/29.)*np.cos(yy/37.)
        source=rim_scene(thickness=1,texture=20.,noise=0.)
        raw=background+(source-background)*(12./70.)+np.random.default_rng(41).normal(0,.3,source.shape)
        raw=np.clip(raw,0,255).astype(np.uint8)
        pool,_=old.extract_descriptors_with_rejections(old.segment_features(raw))
        with contextlib.redirect_stdout(io.StringIO()):
            previous,*_=old._select_and_measure_bottom_slots_v19(raw,pool,1.,'slot.png','single','007','slot210')
        self.assertEqual(previous,[])
        rows,*_=run(raw)
        self.assertEqual(len(rows),3)
        np.testing.assert_allclose([[r['x_width_nm'],r['y_height_nm']] for r in rows],[[17.,21.]]*3,atol=1.3)

    def test_rim_grid_removes_large_curved_background_candidates(self):
        yy,xx=np.indices((400,360))
        background=150.+.03*xx+.02*yy+20.*np.sin(xx/29.)*np.cos(yy/37.)
        source=rim_scene(thickness=1,texture=20.,noise=0.)
        raw=background+(source-background)*(20./70.)+np.random.default_rng(41).normal(0,1.,source.shape)
        rows,_,_,rejects,meta=run(np.clip(raw,0,255).astype(np.uint8))
        self.assertEqual(len(rows),3)
        self.assertTrue(meta['v40_slit_rim_grid_locked'])
        self.assertGreater(meta['v40_slit_rim_grid_rejected_count'],0)
        self.assertTrue(any(r['reason_code']=='slot_outside_bright_rim_grid' for r in rejects))
        np.testing.assert_allclose([r['center_image_x_px'] for r in rows],[60.,180.,300.],atol=1.)

    def test_halo_thickness_and_brightness_do_not_become_the_aperture_size(self):
        measurements=[]
        for thickness,amplitude in ((3,35.),(7,75.)):
            rims,_,_ = sem._v40_bright_rim_candidates(rim_scene(thickness=thickness,amplitude=amplitude),'slot',1.)
            bottom = [d for d in rims if d.center[1]>300]
            self.assertEqual(len(bottom),3)
            values=np.array([np.ptp(d.contour[:,0],axis=0) for d in bottom])
            np.testing.assert_allclose(values,[[17.,21.]]*3,atol=.5)
            measurements.append(values)
        np.testing.assert_allclose(*measurements,atol=.5)

    def test_via_inner_contour_is_used_once_per_lattice_site(self):
        rows,selected,_,bundle = run(rim_scene('via',thickness=4),'via',scale=2.)
        self.assert_rim_rows(rows,20)
        self.assertEqual(len(selected),20)
        self.assertEqual(bundle['meta']['v40_via_rim_measured_count'],20)
        np.testing.assert_allclose([[r['x_width_nm'],r['y_height_nm']] for r in rows],[[80.,80.]]*20,atol=2.)
        sites={(round(r['center_image_x_px']/5),round(r['center_image_y_px']/5)) for r in rows}
        self.assertEqual(len(sites),20)

    def test_one_isolated_rim_does_not_replace_a_repeated_normal_via_scale(self):
        from test_v38_via_recovery import via_scene
        raw=via_scene().astype(np.float32)
        small=np.zeros(raw.shape,np.uint8)
        cv.circle(small,(82,82),5,1,-1)
        outer=cv.dilate(small,cv.getStructuringElement(cv.MORPH_ELLIPSE,(7,7)))
        raw+=cv.GaussianBlur(50.*(outer-small).astype(np.float32),(0,0),.8)
        raw=np.clip(raw,0,255).astype(np.uint8)
        rows,_,_,bundle=run(raw,'via',pool=[])
        self.assertEqual(len(rows),20)
        self.assertFalse(bundle['meta']['v40_via_rim_scale_used'])
        self.assertGreater(bundle['meta']['via_search_diameter_px'],38.)

    def test_bright_rim_is_absolute_brightness_in_both_legacy_polarities(self):
        raw = rim_scene()
        for polarity in ('dark','bright'):
            with self.subTest(polarity=polarity),patch.object(sem,'FEATURE_POLARITY',polarity):
                rows,*_ = run(raw,pool=[])
                self.assert_rim_rows(rows,3)
                self.assertEqual(sem.FEATURE_POLARITY,polarity)
        self.assertEqual(sem.FEATURE_POLARITY,'dark')

    def test_original_pixels_are_required_even_for_previously_valid_rim_locators(self):
        raw = rim_scene()
        rims,_,_ = sem._v40_bright_rim_candidates(raw,'slot',1.)
        blank = np.full(raw.shape,170,np.float32)
        with patch.object(sem,'_V32_MEASUREMENT_GRAY',blank):
            candidates,_,_ = sem._v40_bright_rim_candidates(raw,'slot',1.)
            self.assertEqual(candidates,[])
            for desc in rims[:3]:
                refined,audit = sem._v32_refine_desc(raw,desc,'slot')
                self.assertIsNone(refined)
                self.assertEqual(audit['edge_support_fraction'],0.)

    def test_blanks_tilts_and_noise_cannot_form_valid_rims(self):
        yy,xx=np.indices((200,240))
        for raw in (np.full(xx.shape,170.),130.+.05*xx+.03*yy,
                    170.+np.random.default_rng(40).normal(0,1.5,xx.shape)):
            for pattern in ('slot','via'):
                with self.subTest(pattern=pattern):
                    candidates,_,_ = sem._v40_bright_rim_candidates(raw.astype(np.float32),pattern,1.)
                    self.assertEqual(candidates,[])

    def test_filled_bright_regions_and_plain_dark_slots_are_not_bright_rings(self):
        for raw in (slit_scene(),255-slit_scene()):
            candidates,_,_ = sem._v40_bright_rim_candidates(raw,'slot',1.)
            self.assertEqual(candidates,[])

    def test_wide_open_rim_gaps_do_not_get_filled_as_measurements(self):
        raw = rim_scene(gap=25)
        candidates,_,_ = sem._v40_bright_rim_candidates(raw,'slot',1.)
        self.assertEqual(candidates,[])
        closed,_,_ = sem._v40_bright_rim_candidates(rim_scene(),'slot',1.)
        for desc in closed[:3]:
            refined,audit = sem._v32_refine_desc(raw,desc,'slot')
            self.assertIsNone(refined)
            self.assertTrue(audit['edge_support_fraction'] < .85 or audit['edge_max_missing_fraction'] > .08)

    def test_mixed_dark_and_rim_slits_supplement_partial_success(self):
        ordinary=tuple((r,0) for r in range(4))
        rows,selected,_,_,_ = run(rim_scene(conventional=ordinary))
        self.assertEqual((len(rows),len(selected)),(3,12))
        self.assertEqual(sum(bool(r.get('v40_bright_rim_used')) for r in rows),2)
        np.testing.assert_allclose([[r['x_width_nm'],r['y_height_nm']] for r in rows],[[17.,21.]]*3,atol=.6)

    def test_missing_upper_sites_are_allowed_but_missing_bottom_is_not_measured(self):
        raw = rim_scene(missing=((0,0),(1,1),(2,2),(3,1)))
        rows,_,bottom,_,_ = run(raw)
        self.assert_rim_rows(rows,2)
        self.assertEqual(len(bottom),2)
        self.assertTrue(all(r['slot_is_partial_column'] for r in rows))
        self.assertLess(np.min(np.abs(np.array([r['center_image_x_px'] for r in rows])-60)),1.)
        self.assertTrue(all(abs(r['center_image_x_px']-180)>20 for r in rows))
        strict,*_ = run(raw,slit_min_column_rows=4)
        self.assertEqual(strict,[])

    def test_a_whole_missing_fourth_row_is_never_promoted_from_row_three(self):
        rows,*_ = run(rim_scene(missing=((3,0),(3,1),(3,2))))
        self.assertEqual(rows,[])

    def test_explicit_size_limits_apply_to_the_inner_aperture(self):
        raw = rim_scene(thickness=7)
        rows,*_ = run(raw,scale=2.,slit_width_nm=34.,slit_length_nm=42.,width_tolerance=.1,length_tolerance=.1)
        self.assert_rim_rows(rows,3)
        wrong,*_ = run(raw,slit_width_nm=5.,width_tolerance=.1)
        self.assertEqual(wrong,[])
        via_rows,*_ = run(rim_scene('via'),'via',via_diameter_nm=15.,width_tolerance=.1)
        self.assertEqual(via_rows,[])

    def test_oval_rims_remain_observed_while_star_boundaries_fail_shape_qc(self):
        ovals,_,_ = sem._v40_bright_rim_candidates(rim_scene('via',shape='ellipse'),'via',1.)
        self.assertEqual(len(ovals),20)
        self.assertTrue(all(d.aspect>2. for d in ovals))
        # Touching star rims can enclose small interstitial proposals. The
        # complete shape/scale/lattice pipeline must not report these as vias.
        rows,*_ = run(rim_scene('via',shape='star'),'via')
        self.assertEqual(rows,[])

    def test_four_ring_holes_keep_strict_sparse_bootstrap_and_three_are_insufficient(self):
        sites=((0,0),(0,1),(1,0),(1,1))
        for only,expected in ((sites,4),(sites[:-1],0)):
            raw=rim_scene('via',only_sites=only)
            with sem._v33_detection_context(),contextlib.redirect_stdout(io.StringIO()):
                rows,_,_,_,_,audit=sem._v39_recover_zero_result(raw,1.,'via.png','single','007','via60')
            self.assertEqual(len(rows),expected)
            self.assertEqual(audit['v39_zero_recovery_used'],bool(expected))
            if rows:self.assert_rim_rows(rows,4)

    def test_real_image_status_exports_rim_method_scale_and_review_flag(self):
        rows,status,_ = process(rim_scene(texture=20.),'slot',scale=2.,
                                 slit_width_nm=34.,slit_length_nm=42.,width_tolerance=.15,length_tolerance=.15)
        self.assert_rim_rows(rows,3)
        self.assertEqual(status['v40_bright_rim_measured_count'],3)
        self.assertEqual(status['source_algorithm'],'V40_bright_rim_on_V39_recovery')
        self.assertTrue(status['annotation_exists'])
        self.assertTrue(all(r['confidence_flag']=='BRIGHT_RIM_REVIEW' for r in rows))
        self.assertTrue(all(r['measurement_confidence_flag']=='BRIGHT_RIM_REVIEW' for r in rows))
        np.testing.assert_allclose([[r['x_width_nm'],r['y_height_nm']] for r in rows],[[34.,42.]]*3,atol=1.)

    def test_actual_worker_single_auto_writes_all_patterns_and_rim_audit(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); sample=root/'input'/'007'; sample.mkdir(parents=True)
            for i,raw in enumerate((analytic_trenches(),rim_scene(),rim_scene('via')),1):
                sem._imwrite_unicode_v19(sample/f'scan_{i}.png',raw)
                (sample/f'scan_{i}.txt').write_text('PixelSize=2\n')
            with contextlib.redirect_stdout(io.StringIO()):
                outputs=sem.main_v40(input_root=sample.parent,output_root=root/'output',
                                     via_pattern_mode='auto',image_timeout_seconds=120.)
            frame=sem.pd.read_csv(outputs['all_objects'])
            self.assertEqual(frame.groupby('pattern').size().to_dict(),{'slot':3,'trench':9,'via':20})
            rims=frame[frame.pattern.isin(['slot','via'])]
            self.assertTrue(rims.v40_bright_rim_used.all())
            self.assertEqual(set(rims.confidence_flag),{'BRIGHT_RIM_REVIEW'})
            status=sem.pd.read_csv(outputs['image_status'])
            self.assertEqual(status.v40_bright_rim_measured_count.tolist(),[0,3,20])
            self.assertTrue(all(Path(p).is_file() for p in status.annotated_path))
            settings=json.loads(outputs['settings'].read_text())
            self.assertTrue(settings['v40_bright_rim_enabled'])
            self.assertIn('-V40-',settings['script_version'])
            self.assertEqual(settings['v40_bright_rim_limits'],sem.V40_RIM_LIMITS)
            self.assertTrue(sem.pd.read_csv(outputs['errors']).empty)


if __name__=='__main__':
    unittest.main()
