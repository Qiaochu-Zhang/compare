"""V39 zero-result recovery must retain measured edges and structural evidence."""
import contextlib
import importlib
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
for version in (38, 39):
    sys.path.insert(0, str(ROOT/f'SEM_0830_0831_compare_V{version}_bundle'))
old = importlib.import_module('sem_before_after_compare_0831_V38')
sem = importlib.import_module('sem_before_after_compare_0831_V39')
np, cv = sem.np, sem.cv2
from test_v38_edges import shaded_slots, shaded_vias, shaded_trenches
from test_v38_slit_fusion import slit_scene
from test_v38_via_recovery import via_scene
from test_v38_morphology import via_scene as ellipses, polluted_slots


def sparse_vias(centers=((60,60),(140,60),(60,140),(140,140)), radius=20, shape=(200,200)):
    yy, xx = np.indices(shape)
    signal = np.zeros(shape, np.float32)
    for center in centers:
        cv.circle(signal, tuple(center), radius, 80., -1)
    return (180.+.02*xx+.02*yy-cv.GaussianBlur(signal,(0,0),.8)).astype(np.uint8)


def process(raw, pattern, module=sem, scale=1., **options):
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp); image = root/'图案 image.png'; ann = root/'ann'; ann.mkdir()
        module._imwrite_unicode_v19(image, raw)
        image.with_suffix('.txt').write_text(f'PixelSize={scale}\n')
        with module._v33_detection_context(options), contextlib.redirect_stdout(io.StringIO()):
            rows, status, rejects = module._v32_process_v25_assigned(
                image,'single','007',2,pattern,
                dict(via='via60',slot='slot210',trench='trench160')[pattern],ann)
        status['annotation_exists'] = Path(status['annotated_path']).is_file()
        return rows, status, rejects


def recover(raw, pattern='via', **options):
    with sem._v33_detection_context(options), contextlib.redirect_stdout(io.StringIO()):
        return sem._v39_recover_zero_result(raw,1.,'image.png','single','007',
                    dict(via='via60',slot='slot210',trench='trench160')[pattern])


class ZeroResultRecoveryTests(unittest.TestCase):
    def assert_raw_measurements(self, rows):
        self.assertTrue(rows)
        self.assertTrue(all(r['edge_signal_source']=='original_image_gray' for r in rows))
        self.assertTrue(all(not r['edge_template_forced'] for r in rows))

    def test_opposite_polarity_recovers_via_slit_and_trench_after_v38_zero(self):
        for pattern, raw, expected, dimensions in (
                ('via',255-shaded_vias(),20,(40.,40.)),
                ('slot',255-shaded_slots(),3,(17.,21.)),
                ('trench',255-shaded_trenches('weak_lower'),9,(20.,60.))):
            with self.subTest(pattern=pattern):
                previous, status38, _ = process(raw,pattern,module=old)
                self.assertEqual(previous,[])
                self.assertEqual(status38['auto_fallback_attempts'],'L1,L2')
                rows,status,_ = process(raw,pattern,scale=2.)
                self.assertEqual(len(rows),expected)
                self.assertEqual(status['v39_feature_polarity'],'bright')
                self.assertEqual(status['v39_zero_recovery_method'],'V39_opposite_polarity')
                self.assertEqual(status['auto_fallback_level_used'],4)
                self.assertTrue(status['annotation_exists'])
                self.assertTrue(all(r['confidence_flag']=='ZERO_RESULT_RECOVERY_REVIEW' for r in rows))
                second = 'tip_gap_y_nm' if pattern=='trench' else 'y_height_nm'
                np.testing.assert_allclose([[r['x_width_nm'],r[second]] for r in rows],
                    [[2*dimensions[0],2*dimensions[1]]]*expected,atol=2.)
                if pattern != 'trench':
                    self.assert_raw_measurements(rows)
                self.assertEqual(sem.FEATURE_POLARITY,'dark')
                self.assertFalse(sem._V39_ZERO_RECOVERY_ACTIVE)

    def test_four_holes_bootstrap_only_after_the_existing_passes_fail(self):
        raw = sparse_vias()
        previous,_,_ = process(raw,'via',module=old)
        self.assertEqual(previous,[])
        rows,status,_ = process(raw,'via',scale=2.,via_diameter_nm=80.,width_tolerance=.1)
        self.assertEqual(len(rows),4)
        self.assertEqual(status['auto_fallback_attempts'],'L1,L2,V39_multiscale_structure')
        self.assertTrue(status['measure_v39_via_sparse_bootstrap'])
        self.assertEqual(status['v39_feature_polarity'],'dark')
        self.assert_raw_measurements(rows)
        np.testing.assert_allclose([[r['x_width_nm'],r['y_height_nm']] for r in rows],[[80.,80.]]*4,atol=2.)
        centers=np.array([[r['center_image_x_px'],r['center_image_y_px']] for r in rows])
        actual=np.array([[60,60],[140,60],[60,140],[140,140]])
        self.assertLess(np.max(np.min(np.linalg.norm(centers[:,None]-actual,axis=2),axis=1)),1.)

    def test_sparse_and_large_canvas_does_not_depend_on_half_percent_occupancy(self):
        raw=sparse_vias(centers=((200,200),(280,200),(200,280),(280,280)),radius=7,shape=(800,800))
        with contextlib.redirect_stdout(io.StringIO()):
            candidates=sem._v39_zero_candidates(raw,'via',1.)
        self.assertEqual(len(candidates),4)
        self.assertTrue(all(abs(2*np.sqrt(d.area/np.pi)-14.)<2 for d in candidates))
        rows,_,_,_,level,audit=recover(raw,via_diameter_nm=14.,width_tolerance=.2)
        self.assertEqual((len(rows),level),(4,3))
        self.assertTrue(audit['v39_zero_recovery_used'])
        self.assert_raw_measurements(rows)
        np.testing.assert_allclose([[r['x_width_nm'],r['y_height_nm']] for r in rows],[[14.,14.]]*4,atol=1.5)

    def test_three_holes_single_row_and_nonrectangular_points_cannot_bootstrap(self):
        cases=(((60,60),(140,60),(60,140)),
               ((40,90),(80,90),(120,90),(160,90)),
               ((50,50),(140,50),(65,145),(125,125)))
        for centers in cases:
            with self.subTest(centers=centers):
                rows,_,_,_,_,audit=recover(sparse_vias(centers=centers))
                self.assertEqual(rows,[])
                self.assertFalse(audit['v39_zero_recovery_used'])

    def test_one_shape_invalid_site_or_user_filtered_site_cannot_complete_two_by_two(self):
        raw=sparse_vias()
        yy,xx=np.indices(raw.shape)
        raw[110:171,110:171]=(180.+.02*xx+.02*yy)[110:171,110:171].astype(np.uint8)
        angles=np.linspace(0,2*np.pi,64,endpoint=False)
        radii=20*(1+.35*np.cos(8*angles))
        points=np.column_stack([140+radii*np.cos(angles),140+radii*np.sin(angles)]).astype(np.int32)
        cv.fillPoly(raw,[points],100)
        rows,_,_,_,_,_=recover(raw)
        self.assertEqual(rows,[])
        rows,_,_,_,_,_=recover(sparse_vias(),via_min_circularity=1.)
        self.assertEqual(rows,[])
        rows,_,_,_,_,_=recover(sparse_vias(),via_diameter_nm=15.,width_tolerance=.1)
        self.assertEqual(rows,[])

    def test_blank_tilt_and_noise_do_not_create_measurements_in_either_polarity(self):
        yy,xx=np.indices((200,240))
        for raw in (np.full(xx.shape,170,np.uint8),
                    (130.+.05*xx+.03*yy).astype(np.uint8),
                    np.clip(170.+np.random.default_rng(39).normal(0,1.5,xx.shape),0,255).astype(np.uint8)):
            for pattern in ('via','slot','trench'):
                with self.subTest(pattern=pattern,std=float(raw.std())):
                    rows,_,_,_,_,audit=recover(raw,pattern)
                    self.assertEqual(rows,[])
                    self.assertFalse(audit['v39_zero_recovery_used'])

    def test_locator_proposals_on_blank_raw_pixels_still_fail_edge_validation(self):
        raw=sparse_vias()
        pool,_=sem.extract_descriptors_with_rejections(sem.segment_features(raw))
        with patch.object(sem,'_v39_zero_candidates',return_value=pool), \
             patch.object(sem,'_V32_MEASUREMENT_GRAY',np.full(raw.shape,180,np.float32)), \
             contextlib.redirect_stdout(io.StringIO()):
            rows,_,_,_,_,_=recover(raw)
        self.assertEqual(rows,[])

    def test_existing_success_skips_all_new_trials_and_preserves_weak_holes(self):
        with patch.object(sem,'_v39_recover_zero_result',side_effect=AssertionError('must not run')) as recovery:
            rows,status,_=process(shaded_vias(),'via')
        recovery.assert_not_called()
        self.assertEqual(len(rows),20)
        self.assertFalse(status['v39_zero_recovery_attempted'])
        self.assertFalse(any(r['v39_zero_recovery_used'] for r in rows))
        self.assert_raw_measurements(rows)

    def test_failure_restores_polarity_sparse_gate_and_threshold_context(self):
        previous=(sem.FEATURE_POLARITY,sem._V39_ZERO_RECOVERY_ACTIVE,sem._V25_ZERO_FALLBACK_LEVEL)
        seen=[]
        def trial(*args):
            seen.append(sem.FEATURE_POLARITY)
            if len(seen)==1:
                return [],None,[],{},0
            raise RuntimeError('trial failed')
        with patch.object(sem,'_v39_zero_trial',side_effect=trial),self.assertRaisesRegex(RuntimeError,'trial failed'):
            recover(sparse_vias())
        self.assertEqual(seen,['dark','bright'])
        self.assertEqual((sem.FEATURE_POLARITY,sem._V39_ZERO_RECOVERY_ACTIVE,sem._V25_ZERO_FALLBACK_LEVEL),previous)

    def test_opposite_trial_respects_an_explicit_bright_default_and_restores_it(self):
        with patch.object(sem,'FEATURE_POLARITY','bright'):
            rows,_,_,_,level,audit=recover(sparse_vias())
            self.assertEqual(sem.FEATURE_POLARITY,'bright')
        self.assertEqual((len(rows),level),(4,4))
        self.assertEqual(audit['v39_feature_polarity'],'dark')
        self.assertEqual(sem.FEATURE_POLARITY,'dark')

    def test_filtering_one_of_four_valid_contours_invalidates_sparse_bootstrap(self):
        raw=sparse_vias(centers=((60,60),(140,60),(60,140)))
        signal=np.zeros(raw.shape,np.float32)
        cv.ellipse(signal,(140,140),(24,10),0,0,360,80.,-1)
        raw=(raw.astype(float)-cv.GaussianBlur(signal,(0,0),.8)).astype(np.uint8)
        rows,_,rejects,_,_,_=recover(raw,via_min_circularity=.85)
        self.assertEqual(rows,[])
        self.assertTrue(any(r['reason_code']=='via_sparse_array_not_fully_supported' for r in rejects))

    def test_actual_worker_exports_four_holes_and_recovery_audit(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); image=root/'稀疏 image.png'; ann=root/'ann'; ann.mkdir()
            sem._imwrite_unicode_v19(image,sparse_vias())
            image.with_suffix('.txt').write_text('PixelSize=2\n')
            with contextlib.redirect_stdout(io.StringIO()):
                rows,status,_,_,_=sem._v32_process_image_with_timeout(
                    image,'single','007',2,'via','via60',ann,timeout_seconds=120.,
                    detection_options={'via_diameter_nm':80.,'width_tolerance':.1})
            self.assertEqual(len(rows),4)
            self.assertTrue(status['v39_zero_recovery_used'])
            self.assertTrue(status['measure_v39_via_sparse_bootstrap'])
            self.assertTrue(Path(status['annotated_path']).is_file())
            self.assertEqual(status['source_algorithm'],'V39_zero_result_recovery_on_V38_fusion')


class PreservedRecognitionTests(unittest.TestCase):
    def test_v38_curved_dense_vias_and_smooth_ovals_keep_their_counts(self):
        for raw,expected in ((via_scene(texture=8.),20),(via_scene(pitch=50),20),(ellipses(35.),18)):
            with self.subTest(expected=expected),contextlib.redirect_stdout(io.StringIO()):
                pool,_=sem.extract_descriptors_with_rejections(sem.segment_features(raw))
                rows,_,_,_=sem._select_and_measure_vias_v19(raw,pool,1.,'via.png','single','007','via60')
            self.assertEqual(len(rows),expected)
            self.assertTrue(all(r['edge_signal_source']=='original_image_gray' for r in rows))

    def test_v38_curved_dense_and_incomplete_slit_columns_are_preserved(self):
        for raw,expected in ((slit_scene(texture=8.),3),(slit_scene(width=64,pitch=80),3),
                             (polluted_slots(missing=((0,0),(1,1),(2,2))),3)):
            with self.subTest(expected=expected),contextlib.redirect_stdout(io.StringIO()):
                pool,_=sem.extract_descriptors_with_rejections(sem.segment_features(raw))
                rows,_,_,_,_=sem._select_and_measure_bottom_slots_v19(raw,pool,1.,'slot.png','single','007','slot210')
            self.assertEqual(len(rows),expected)
            self.assertTrue(all(r['edge_signal_source']=='original_image_gray' for r in rows))


if __name__=='__main__':
    unittest.main()
