"""V34 split tolerances: measured edges, independent limits and worker settings."""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'SEM_0830_0831_compare_V34_bundle'))
import sem_before_after_compare_0831_V34 as sem
from test_v30 import gradient_slots
from test_v33 import analytic_trenches
from test_v33_controls import contour_desc

np = sem.np


class ToleranceOptionsTests(unittest.TestCase):
    def options_from_cli(self, *arguments):
        args = sem._v32_build_parser().parse_args(arguments)
        return sem.V34DetectionOptions(**{name: getattr(args, name) for name in sem.asdict(sem.V34DetectionOptions())})

    def test_defaults_and_explicit_cli_overrides(self):
        for options in (sem.V34DetectionOptions(), self.options_from_cli()):
            self.assertEqual((options.width_tolerance, options.length_tolerance), (.5, 1.))
            self.assertIsNone(options.size_tolerance)
        options = self.options_from_cli('--width-tolerance', '.25', '--gap-tolerance', '.8')
        self.assertEqual((options.width_tolerance, options.length_tolerance), (.25, .8))

    def test_legacy_uniform_tolerance_and_specific_precedence_survive_serialization(self):
        cases = [(['--size-tolerance', '.3'], (.3, .3)),
                 (['--width-tolerance', '.4', '--size-tolerance', '.3'], (.4, .3)),
                 (['--size-tolerance', '.3', '--length-tolerance', '1'], (.3, 1.)),
                 (['--width-tolerance', '.5', '--length-tolerance', '1', '--size-tolerance', '.3'], (.5, 1.))]
        for arguments, expected in cases:
            with self.subTest(arguments=arguments):
                options = self.options_from_cli(*arguments)
                self.assertEqual((options.width_tolerance, options.length_tolerance), expected)
                restored = sem._v33_options(json.loads(json.dumps(sem.asdict(options))))
                self.assertEqual(restored, options)

    def test_invalid_tolerances_fail_before_output_creation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for field in ('width_tolerance', 'length_tolerance', 'size_tolerance'):
                for value in (0, -.1, 1.01, float('nan'), float('inf'), True, '0.5'):
                    with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                        sem.main_v34(root/'before', root/'after', root/'out', detection_options={field: value})
            self.assertFalse((root/'out').exists())

    def test_ranges_include_boundaries_but_never_zero_or_nonfinite_sizes(self):
        with sem._v33_detection_context():
            for value, accepted in ((49.99, False), (50., True), (150., True), (150.01, False)):
                self.assertEqual(sem._v33_size_in_range(value, 100.), accepted)
            for value, accepted in ((0., False), (-1., False), (.01, True), (200., True), (200.01, False),
                                    (float('nan'), False), (float('inf'), False)):
                self.assertEqual(sem._v33_size_in_range(value, 100., tolerance=sem._V33_OPTIONS.length_tolerance), accepted)


class MeasuredToleranceTests(unittest.TestCase):
    def test_slit_length_can_differ_by_over_50_percent_in_both_directions(self):
        gray = gradient_slots(40.)
        for length_hint in (12., 50.):
            for tolerance, expected in ((1., 3), (.5, 0)):
                with self.subTest(hint=length_hint, tolerance=tolerance), contextlib.redirect_stdout(io.StringIO()), \
                        sem._v33_detection_context(dict(slit_width_nm=17., slit_length_nm=length_hint,
                                                        length_tolerance=tolerance)):
                    rows, selected, bottom, _, _ = sem._select_and_measure_bottom_slots_v19(
                        gray, [], 1., 'slit.png', 'before', '1', 'slot210')
                self.assertEqual(len(rows), expected)
                self.assertEqual(len(selected), 4*expected)
                self.assertEqual(len(bottom), expected)
                if rows:
                    np.testing.assert_allclose([r['y_height_nm'] for r in rows], 21., atol=1.)
                    np.testing.assert_allclose([r['x_width_nm'] for r in rows], 17., atol=1.)

    def test_trench_gap_can_differ_by_over_50_percent_in_both_directions(self):
        gray = analytic_trenches()
        for gap_hint in (35., 140.):
            for tolerance, expected in ((1., 9), (.5, 0)):
                with self.subTest(hint=gap_hint, tolerance=tolerance), contextlib.redirect_stdout(io.StringIO()), \
                        sem._v33_detection_context(dict(trench_width_nm=20., trench_gap_nm=gap_hint,
                                                        length_tolerance=tolerance)):
                    rows, pairs, _, _ = sem._select_and_measure_trenches_v19(
                        gray, [], 1., 'trench.png', 'before', '1', 'trench160')
                self.assertEqual(len(rows), expected)
                self.assertEqual(len(pairs), expected)
                if rows:
                    np.testing.assert_allclose([r['tip_gap_y_nm'] for r in rows], 60., atol=1.)
                    np.testing.assert_allclose([r['x_width_nm'] for r in rows], 20., atol=1.)

    def test_default_length_tolerance_does_not_relax_via_or_slit_width(self):
        desc = contour_desc(radius=20.)
        for width_tolerance, expected in ((.5, 0), (.8, 1)):
            with sem._v33_detection_context(dict(via_diameter_nm=24., width_tolerance=width_tolerance)):
                _, rows, _ = sem._v33_filter_vias([desc], [{}], 1., [])
            self.assertEqual(len(rows), expected)
        with sem._v33_detection_context(dict(slit_width_nm=10., slit_length_nm=12.)), \
                contextlib.redirect_stdout(io.StringIO()):
            rows, _, _, rejects, _ = sem._select_and_measure_bottom_slots_v19(
                gradient_slots(40.), [], 1., 'slit.png', 'before', '1', 'slot210')
        self.assertEqual(rows, [])
        self.assertTrue(any(r['reason_code'] == 'slit_size_outside_hint_range' for r in rejects))

    def test_blank_images_remain_rejected_with_100_percent_length_tolerance(self):
        gray = np.full((400, 430), 160, np.uint8)
        with sem._v33_detection_context(dict(slit_width_nm=17., slit_length_nm=50.,
                                            trench_width_nm=20., trench_gap_nm=140.)), \
                contextlib.redirect_stdout(io.StringIO()):
            for key, function in [('slot210', sem._select_and_measure_bottom_slots_v19),
                                  ('trench160', sem._select_and_measure_trenches_v19)]:
                with self.subTest(pattern=key):
                    self.assertEqual(function(gray, [], 1., 'blank.png', 'after', '1', key)[0], [])


class WorkerAndExportTests(unittest.TestCase):
    def test_real_worker_uses_independent_tolerances_and_preserves_nm_calibration(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root/'trench.png'
            sem._imwrite_unicode_v19(path, analytic_trenches())
            path.with_suffix('.txt').write_text('PixelSize=2.0')
            rows, status, _, _, _ = sem._v32_process_image_with_timeout(
                path, 'before', '007', 2, 'trench', 'trench160', root/'ann', 60.,
                dict(trench_width_nm=40., trench_gap_nm=70., width_tolerance=.25, length_tolerance=.8))
            self.assertEqual(len(rows), 9)
            np.testing.assert_allclose([r['tip_gap_y_nm'] for r in rows], 120., atol=2.)
            self.assertTrue(all(r['input_width_tolerance'] == .25 and r['input_length_tolerance'] == .8 for r in rows))
            self.assertTrue(Path(status['annotated_path']).exists())

    def test_batch_records_resolved_defaults_in_settings_and_each_job(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for stage in ('before', 'after'):
                folder = root/stage/'01'
                folder.mkdir(parents=True)
                for i in range(3):
                    (folder/f'{i}.png').touch()
                    (folder/f'{i}.txt').write_text('PixelSize=1')
            calls = []
            def worker(**kwargs):
                calls.append(kwargs)
                return [], dict(stage=kwargs['stage'], condition='01', region=2, pattern=kwargs['public_pattern'],
                                image=kwargs['image_path'].name, status='OK'), [], kwargs['v25_pattern_key'], []
            with patch.object(sem, '_v32_process_image_with_timeout', side_effect=worker), \
                    patch.object(sem, '_v32_save_pair_plots'), contextlib.redirect_stdout(io.StringIO()):
                outputs = sem.main_v34(root/'before', root/'after', root/'out')
            self.assertEqual(len(calls), 6)
            self.assertTrue(all((c['detection_options'].width_tolerance, c['detection_options'].length_tolerance) == (.5, 1.) for c in calls))
            settings = json.loads(outputs['settings'].read_text())
            self.assertEqual(settings['detection_options']['width_tolerance'], .5)
            self.assertEqual(settings['detection_options']['length_tolerance'], 1.)
            self.assertIn('V34', settings['script_version'])
            self.assertIn('V34', outputs['excel'].name)


if __name__ == '__main__':
    unittest.main()
