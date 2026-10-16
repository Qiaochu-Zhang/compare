"""V38 R2: recover V35 via recall while retaining weak and elongated holes."""
import contextlib
import importlib
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
MODULES={}
for v in (35,37,38):
    sys.path.insert(0,str(ROOT/f'SEM_0830_0831_compare_V{v}_bundle'))
    MODULES[v]=importlib.import_module(f'sem_before_after_compare_0831_V{v}')
sem=MODULES[38]
np,cv=sem.np,sem.cv2
from test_v38_edges import partial_feature


def via_scene(texture=0.,pitch=75):
    yy,xx=np.indices((3*pitch+90,4*pitch+90))
    background=170.+.015*xx+.02*yy+texture*np.sin(xx/29.)*np.cos(yy/37.)
    signal=np.zeros(xx.shape,np.float32)
    for r in range(4):
        for c in range(5):
            cv.circle(signal,(45+c*pitch,45+r*pitch),20,80.,-1)
    raw=background-cv.GaussianBlur(signal,(0,0),.8)
    raw+=np.random.default_rng(39).normal(0,.2,xx.shape)
    return np.clip(raw,0,255).astype(np.uint8)


def run(raw,version=38,scale=1.,empty_pool=False,**options):
    module=MODULES[version]
    pool,_=module.extract_descriptors_with_rejections(module.segment_features(raw))
    with contextlib.redirect_stdout(io.StringIO()), module._v33_detection_context(options):
        return module._select_and_measure_vias_v19(raw,[] if empty_pool else pool,
                                                   scale,'via.png','single','007','via60')


class ViaRecoveryTests(unittest.TestCase):
    def assert_twenty_real_holes(self,rows,pitch=75,scale=1.):
        self.assertEqual(len(rows),20)
        centers=np.array([[r['center_image_x_px'],r['center_image_y_px']] for r in rows])
        expected=np.array([[45+c*pitch,45+r*pitch] for r in range(4) for c in range(5)])
        distance=np.linalg.norm(centers[:,None,:]-expected[None,:,:],axis=2)
        self.assertEqual(len(set(np.argmin(distance,axis=1))),20)
        self.assertLess(float(distance.min(axis=1).max()),1.)
        np.testing.assert_allclose([[r['x_width_nm'],r['y_height_nm']] for r in rows],
                                   [[40.*scale,40.*scale]]*20,atol=1.5*scale)
        self.assertTrue(all(r['edge_signal_source']=='original_image_gray' for r in rows))
        self.assertTrue(all(not r['edge_template_forced'] for r in rows))

    def test_curved_background_recovers_twenty_instead_of_v37_six(self):
        raw=via_scene(texture=8.)
        old,_,_,_=run(raw,35)
        lost,_,rejected,_=run(raw,37)
        rows,_,_,bundle=run(raw)
        self.assertEqual((len(old),len(lost)),(20,6))
        self.assertTrue(any(r.get('edge_status')=='UNRESOLVED_OUTER_BOUNDARY' for r in rejected))
        self.assert_twenty_real_holes(rows)
        self.assertTrue(bundle['meta']['v38_via_primary_scale_locked'])
        self.assertTrue(all(r['via_boundary_rule'].startswith('V38_bracketed_') for r in rows))

    def test_dense_candidates_do_not_replace_the_primary_diameter_with_gap_fragments(self):
        raw=via_scene(pitch=50)
        old,_,_,_=run(raw,35)
        lost,_,_,old_bundle=run(raw,37)
        rows,_,_,bundle=run(raw)
        self.assertEqual((len(old),len(lost)),(20,0))
        self.assertLess(old_bundle['meta']['via_search_diameter_px'],15.)
        self.assertGreater(bundle['meta']['via_search_diameter_px'],38.)
        self.assertEqual(bundle['meta']['via_scale_locator'],'V35_primary_near_round')
        self.assert_twenty_real_holes(rows,pitch=50)

    def test_supplemental_search_still_works_with_no_primary_candidates(self):
        rows,_,_,bundle=run(via_scene(),empty_pool=True)
        self.assert_twenty_real_holes(rows)
        self.assertFalse(bundle['meta']['v38_via_primary_scale_locked'])
        self.assertEqual(bundle['meta']['via_scale_locator'],'V37_supplemental_scale')

    def test_via_specific_boundary_path_keeps_weak_lower_halves_in_both_polarities(self):
        for weak in (6.,12.):
            for polarity in ('dark','bright'):
                raw,seed=partial_feature('via',weak=weak)
                with self.subTest(weak=weak,polarity=polarity),patch.object(sem,'FEATURE_POLARITY',polarity):
                    desc,audit=sem._v32_refine_desc(raw if polarity=='dark' else 255-raw,
                                                   sem.describe_contour(seed),'via')
                self.assertIsNotNone(desc,audit)
                np.testing.assert_allclose(np.ptp(desc.contour[:,0],axis=0),[45.,77.],atol=1.)
                self.assertTrue(audit['boundary_relocalized'])

    def test_bright_polarity_scale_hints_and_explicit_filters_stay_effective(self):
        with patch.object(sem,'FEATURE_POLARITY','bright'):
            rows,_,_,bundle=run(255-via_scene(texture=8.),scale=2.,via_diameter_nm=80.,
                                via_min_circularity=.8,via_max_size_deviation=.1,width_tolerance=.15)
        self.assert_twenty_real_holes(rows,scale=2.)
        self.assertEqual(bundle['meta']['via_scale_locator'],'user_size_hint')
        rejected,_,reasons,_=run(via_scene(),via_min_circularity=1.)
        self.assertEqual(rejected,[])
        self.assertTrue(any(r['reason_code']=='via_circularity_below_threshold' for r in reasons))

    def test_locators_on_blank_measurement_pixels_cannot_create_holes(self):
        raw=via_scene(pitch=50)
        with patch.object(sem,'_V32_MEASUREMENT_GRAY',np.full(raw.shape,180,np.float32)):
            rows,_,rejected,_=run(raw)
        self.assertEqual(rows,[])
        self.assertTrue(any(r['reason_code']=='original_image_edge_not_supported' for r in rejected))

    def test_real_worker_preserves_dense_holes_scale_and_audit_fields(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);image=root/'密集 via.png';ann=root/'ann';ann.mkdir()
            sem._imwrite_unicode_v19(image,via_scene(pitch=50))
            image.with_suffix('.txt').write_text('PixelSize=2.0\n')
            with contextlib.redirect_stdout(io.StringIO()):
                rows,status,rejected,key,trials=sem._v32_process_image_with_timeout(
                    image,'single','007',2,'via','via60',ann,timeout_seconds=120.,
                    detection_options=dict(via_min_circularity=.8,via_max_size_deviation=.1))
            self.assert_twenty_real_holes(rows,pitch=50,scale=2.)
            self.assertEqual((key,status['n_measurements']),('via60',20))
            self.assertEqual(status['measure_via_scale_locator'],'V35_primary_near_round')
            self.assertTrue(status['measure_v38_via_primary_scale_locked'])
            self.assertTrue(Path(status['annotated_path']).is_file())

    def test_v37_elongated_holes_and_abnormal_shape_rejection_are_preserved(self):
        from test_v38_morphology import via_scene as ellipses
        raw=ellipses(35.)
        old,_,_,_=run(raw,35)
        v37,_,_,_=run(raw,37)
        rows,selected,rejected,_=run(raw)
        self.assertLess(len(old),18)
        self.assertEqual((len(v37),len(rows),len(selected)),(18,18,18))
        self.assertTrue(all(np.linalg.norm(d.center-[230,230])>30 for d in selected))
        self.assertTrue(all(np.linalg.norm(d.center-[145,145])>30 for d in selected))
        self.assertTrue(rejected)


if __name__=='__main__':
    unittest.main()
