"""V33 user controls: real contours, weak images, calibration and worker propagation."""
import contextlib
import io
import json
import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from test_v33 import sem, analytic_vias, analytic_trenches
from test_v30 import gradient_slots

np, cv = sem.np, sem.cv2


def contour_desc(radius=20., center=(50., 50.), y_factor=1.):
    angle = np.linspace(0, 2*math.pi, 128, endpoint=False)
    points = np.column_stack([center[0]+radius*np.cos(angle), center[1]+radius*y_factor*np.sin(angle)])
    return sem.describe_contour(points.astype(np.float32).reshape(-1, 1, 2))


def varied_vias():
    yy, xx = np.mgrid[:360, :430]
    distance = np.full(xx.shape, 1000.)
    for iy, y in enumerate((65, 140, 215, 290)):
        for ix, x in enumerate((65, 140, 215, 290, 365)):
            radius = 28. if (iy, ix) == (2, 2) else 20.
            distance = np.minimum(distance, np.hypot(xx-x, yy-y)-radius)
    return (180.-120./(1.+np.exp(np.clip(distance/.8, -100, 100)))).astype(np.uint8)


class OptionsTests(unittest.TestCase):
    def test_cli_exposes_all_parameters_and_slot_aliases(self):
        parser = sem._v32_build_parser()
        args = parser.parse_args(['--via-min-circularity', '.82', '--via-max-size-deviation', '.2',
                                  '--via-diameter-nm', '80', '--slot-width-nm', '17',
                                  '--slot-length-nm', '21', '--trench-width-nm', '20',
                                  '--trench-gap-nm', '60', '--size-tolerance', '.4'])
        options = sem.V33DetectionOptions(**{key: getattr(args, key) for key in sem.asdict(sem.V33DetectionOptions())})
        self.assertEqual(options.via_min_circularity, .82)
        self.assertEqual(options.via_max_size_deviation, .2)
        self.assertEqual(options.via_diameter_nm, 80.)
        self.assertEqual((options.slit_width_nm, options.slit_length_nm), (17., 21.))
        self.assertEqual((options.trench_width_nm, options.trench_gap_nm), (20., 60.))
        self.assertEqual(options.size_tolerance, .4)

    def test_invalid_options_fail_before_creating_outputs(self):
        cases = [('via_min_circularity', -.01), ('via_min_circularity', 1.1),
                 ('via_max_size_deviation', -1), ('via_max_size_deviation', 1),
                 ('size_tolerance', 0), ('size_tolerance', 1)]
        for key in sem.asdict(sem.V33DetectionOptions()):
            cases.extend((key, value) for value in (float('nan'), float('inf'), True))
            if key.endswith('_nm'):
                cases.extend(((key, 0), (key, -4)))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for key, value in cases:
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    sem.main_v33(root/'before', root/'after', root/'out', detection_options={key: value})
            self.assertFalse((root/'out').exists())

    def test_context_and_calibration_restore_after_exception(self):
        original = sem._V33_OPTIONS
        with self.assertRaises(ValueError), sem._v33_detection_context({'slit_width_nm': 1000.}):
            sem._select_and_measure_bottom_slots_v19(np.full((80, 80), 150, np.uint8), [], 1.,
                                                    'bad.png', 'before', '1', 'slot210')
        self.assertIs(sem._V33_OPTIONS, original)
        self.assertIsNone(sem._V33_PIXEL_SIZE_NM)


class ViaQualityTests(unittest.TestCase):
    def test_final_circularity_threshold_is_adjustable(self):
        circle, oval = contour_desc(), contour_desc(y_factor=.55)
        self.assertGreater(circle.circularity, .99)
        self.assertLess(oval.circularity, .95)
        for limit, expected in ((.95, 1), (.7, 2)):
            rejects = []
            with sem._v33_detection_context({'via_min_circularity': limit, 'via_max_size_deviation': 0}):
                selected, rows, meta = sem._v33_filter_vias([circle, oval], [{}, {}], 1., rejects)
            self.assertEqual(len(rows), expected)
            self.assertEqual(len(selected), expected)
            self.assertEqual(meta['v33_via_circularity_rejected'], 2-expected)
            if rejects:
                self.assertEqual(rejects[0]['reason_code'], 'via_circularity_below_threshold')
                self.assertAlmostEqual(rejects[0]['circularity'], oval.circularity)

    def test_same_image_median_removes_one_large_via_without_clipping_between_stages(self):
        for scale in (1., 2.5):
            descs = [contour_desc(radius=r*scale) for r in (20., 20., 20., 28.)]
            rejects = []
            with sem._v33_detection_context({'via_max_size_deviation': .3}):
                selected, rows, meta = sem._v33_filter_vias(descs, [{} for _ in descs], 1., rejects)
            self.assertEqual(len(selected), 3)
            self.assertEqual(rejects[0]['reason_code'], 'via_size_outlier')
            self.assertAlmostEqual(rejects[0]['via_size_deviation_fraction'], .4, places=5)
            self.assertAlmostEqual(meta['via_image_median_diameter_nm'], 40.*scale, delta=.03)
            self.assertTrue(all(r['via_size_filter_applied'] for r in rows))

    def test_size_filter_disabled_or_insufficient_peers(self):
        descs = [contour_desc(radius=r) for r in (20., 20., 20., 28.)]
        for fraction, items in ((0., descs), (.05, [descs[0], descs[-1]])):
            with sem._v33_detection_context({'via_max_size_deviation': fraction}):
                selected, rows, meta = sem._v33_filter_vias(items, [{} for _ in items], 1., [])
            self.assertEqual(len(selected), len(items))
            self.assertFalse(meta['via_size_filter_applied'])

    def test_real_via_pipeline_removes_outlier_and_refits_annotations(self):
        gray = varied_vias()
        pool, _ = sem.extract_descriptors_with_rejections(sem.segment_features(gray))
        results = []
        for fraction in (0., .2):
            with sem._v33_detection_context({'via_max_size_deviation': fraction}), contextlib.redirect_stdout(io.StringIO()):
                results.append(sem._select_and_measure_vias_v19(gray, pool, 1., 'via.png', 'before', '1', 'via60'))
        self.assertEqual(len(results[0][0]), 20)
        rows, selected, rejects, bundle = results[1]
        self.assertEqual(len(rows), 19)
        self.assertEqual(len(selected), 19)
        self.assertEqual([r['object_id'] for r in rows], list(range(1, 20)))
        self.assertEqual(sum(r['reason_code'] == 'via_size_outlier' for r in rejects), 1)
        self.assertEqual(bundle['meta']['via_v24_final_count'], 19)
        self.assertEqual(sem._annotate_vias_v19(gray, selected, rows, bundle).shape[:2], gray.shape)
        normalized = sem.pd.DataFrame(rows).assign(pattern='via', region=2, pair_key='C1_R2_via')
        metrics = sem._v32_build_image_metrics(normalized)
        self.assertFalse(metrics.empty)
        self.assertTrue((metrics['n'] == 19).all())


class SizeHintTests(unittest.TestCase):
    def test_blurred_slots_gain_supported_columns_without_forcing_dimensions(self):
        gray = cv.GaussianBlur(gradient_slots(5.).astype(np.float32), (0, 0), 3.5).astype(np.uint8)
        with contextlib.redirect_stdout(io.StringIO()):
            baseline = sem._select_and_measure_bottom_slots_v19(gray, [], 1., 'blur.png', 'before', '1', 'slot210')
            with sem._v33_detection_context({'slit_width_nm': 17., 'slit_length_nm': 21.}):
                rows, selected, bottom, _, meta = sem._select_and_measure_bottom_slots_v19(
                    gray, [], 1., 'blur.png', 'before', '1', 'slot210')
        self.assertGreater(len(rows), len(baseline[0]))
        self.assertGreater(meta['v33_slot_hint_candidates'], 0)
        self.assertEqual(len(selected), 4*len(rows))
        self.assertEqual(len(bottom), len(rows))
        for r in rows:
            self.assertGreaterEqual(r['edge_support_fraction'], .85)
            self.assertFalse(r['edge_template_forced'])
            self.assertAlmostEqual(r['x_width_nm'], 17., delta=2.)
            self.assertAlmostEqual(r['y_height_nm'], 21., delta=3.)
            self.assertEqual(r['edge_signal_source'], 'original_image_gray')

    def test_slit_hints_apply_only_to_bottom_length_and_via_circularity_does_not_affect_it(self):
        gray = gradient_slots(40.)
        opts = dict(slit_width_nm=13.6, slit_length_nm=16.8, size_tolerance=.1, via_min_circularity=1.)
        with sem._v33_detection_context(opts), contextlib.redirect_stdout(io.StringIO()):
            rows, selected, bottom, _, _ = sem._select_and_measure_bottom_slots_v19(gray, [], .8, 'slot.png', 'before', '1', 'slot210')
        self.assertEqual(len(rows), 3)
        self.assertEqual(len(selected), 12)
        self.assertTrue(any(d.length_px*.8 > 2.*opts['slit_length_nm'] for d in selected))
        np.testing.assert_allclose([r['x_width_nm'] for r in rows], 13.6, atol=.8)
        with sem._v33_detection_context(dict(slit_length_nm=50., size_tolerance=.1)), contextlib.redirect_stdout(io.StringIO()):
            rows, selected, bottom, rejects, _ = sem._select_and_measure_bottom_slots_v19(gray, [], .8, 'slot.png', 'before', '1', 'slot210')
        self.assertEqual((rows, selected, bottom), ([], [], []))
        self.assertTrue(any(r['reason_code'] == 'slit_size_outside_hint_range' for r in rejects))

    def test_trench_priors_guide_edges_but_do_not_set_measured_width_and_gap(self):
        opts = dict(trench_width_nm=44., trench_gap_nm=130., via_min_circularity=1.)
        with sem._v33_detection_context(opts), contextlib.redirect_stdout(io.StringIO()):
            rows, pairs, _, _ = sem._select_and_measure_trenches_v19(
                analytic_trenches(), [], 2., 'trench.png', 'before', '1', 'trench160')
        self.assertEqual(len(rows), 9)
        self.assertEqual(len(pairs), 9)
        np.testing.assert_allclose([r['x_width_nm'] for r in rows], 40., atol=1.5)
        np.testing.assert_allclose([r['tip_gap_y_nm'] for r in rows], 120., atol=2.)
        self.assertTrue(all(r['input_trench_width_nm'] == 44 for r in rows))

    def test_size_hints_cannot_create_objects_in_blank_images(self):
        gray = np.full((400, 430), 160, np.uint8)
        opts = dict(via_diameter_nm=40., slit_width_nm=17., slit_length_nm=21.,
                    trench_width_nm=20., trench_gap_nm=60.)
        with sem._v33_detection_context(opts), contextlib.redirect_stdout(io.StringIO()):
            for key, function in [('via60', sem._select_and_measure_vias_v19),
                                  ('slot210', sem._select_and_measure_bottom_slots_v19),
                                  ('trench160', sem._select_and_measure_trenches_v19)]:
                with self.subTest(pattern=key):
                    result = function(gray, [], 1., 'blank.png', 'after', '1', key)
                    self.assertEqual(result[0], [])


class WorkerOptionsTests(unittest.TestCase):
    def test_real_worker_preserves_user_options_and_measures_actual_via_diameter(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root/'via.png'
            sem._imwrite_unicode_v19(path, analytic_vias()[0])
            path.with_suffix('.txt').write_text('PixelSize=2.0')
            rows, status, _, _, _ = sem._v32_process_image_with_timeout(
                path, 'before', '007', 2, 'via', 'via60', root/'ann', 60.,
                dict(via_diameter_nm=90., via_min_circularity=.85, via_max_size_deviation=.2))
            self.assertEqual(len(rows), 20)
            np.testing.assert_allclose([r['x_width_nm'] for r in rows], 80., atol=2.)
            for row in rows:
                self.assertEqual(row['via_search_scale_source'], 'user_size_hint')
                self.assertEqual(row['via_search_diameter_nm'], 90.)
                self.assertEqual(row['input_via_min_circularity'], .85)
                self.assertEqual(row['input_via_max_size_deviation'], .2)
            self.assertTrue(Path(status['annotated_path']).exists())

    def test_fallback_never_weakens_user_circularity_threshold(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root/'via.png'
            sem._imwrite_unicode_v19(path, analytic_vias()[0])
            path.with_suffix('.txt').write_text('PixelSize=1.0')
            rows, status, rejects, _, _ = sem._v32_process_image_with_timeout(
                path, 'after', '1', 2, 'via', 'via60', root/'ann', 60.,
                dict(via_min_circularity=1., via_diameter_nm=40.))
            self.assertEqual(rows, [])
            self.assertTrue(status['auto_fallback_triggered'])
            failures = [r for r in rejects if r['reason_code'] == 'via_circularity_below_threshold']
            self.assertGreaterEqual(len(failures), 20)
            self.assertTrue(all(r['input_via_min_circularity'] == 1. for r in failures))
            self.assertEqual(status['n_measurements'], 0)

    def test_batch_passes_options_to_every_pattern_and_exports_settings(self):
        opts = sem.V33DetectionOptions(via_diameter_nm=90., slit_width_nm=14., trench_gap_nm=60.)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for stage in ('before', 'after'):
                folder = root/stage/'01'
                folder.mkdir(parents=True)
                for i in range(3):
                    (folder/f'{i}.png').touch()
                    (folder/f'{i}.txt').write_text('PixelSize=2')
            calls = []
            def worker(**kw):
                calls.append(kw)
                return [], dict(stage=kw['stage'], condition='01', region=2, pattern=kw['public_pattern'],
                                image=kw['image_path'].name, status='OK'), [], kw['v25_pattern_key'], []
            with patch.object(sem, '_v32_process_image_with_timeout', side_effect=worker), \
                    patch.object(sem, '_v32_save_pair_plots'), contextlib.redirect_stdout(io.StringIO()):
                outputs = sem.main_v33(root/'before', root/'after', root/'out', detection_options=opts)
            self.assertEqual(len(calls), 6)
            self.assertTrue(all(c['detection_options'] == opts for c in calls))
            settings = json.loads(outputs['settings'].read_text())
            self.assertEqual(settings['detection_options'], sem.asdict(opts))
            self.assertIn('V33', outputs['excel'].name)


if __name__ == '__main__':
    unittest.main()
