"""V31 signal accuracy, shape preservation, original-gray gates and hard deadlines."""
import contextlib
import io
import json
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'SEM_0830_0831_compare_V31_bundle'))
import sem_before_after_compare_0831_V31 as sem
from test_v30 import gradient_slots
np, cv = sem.np, sem.cv2


def analytic_shape(power=2, angle=0., gradient=True):
    yy, xx = np.mgrid[:160, :180]
    radians = np.deg2rad(angle)
    xa = np.array([np.cos(radians), np.sin(radians)])
    ya = np.array([-np.sin(radians), np.cos(radians)])
    x = (xx-90)*xa[0] + (yy-80)*xa[1]
    y = (xx-90)*ya[0] + (yy-80)*ya[1]
    distance = ((np.abs(x)/22)**power + (np.abs(y)/32)**power)**(1/power)
    gray = 180.-120./(1.+np.exp(np.clip((distance-1)*22/.8, -100, 100)))
    if gradient:
        gray += .10*xx + .08*yy
    contours, _ = cv.findContours((distance<=1).astype(np.uint8)*255, cv.RETR_EXTERNAL, cv.CHAIN_APPROX_SIMPLE)
    return gray.astype(np.float32), contours[0], xa, ya


def analytic_vias():
    yy, xx = np.mgrid[:360, :430]
    distance = np.full(xx.shape, 100.)
    for y in (65, 140, 215, 290):
        for x in (65, 140, 215, 290, 365):
            distance = np.minimum(distance, np.sqrt(((xx-x)/20)**2+((yy-y)/20)**2)-1)
    gray = 170.+.03*xx+.02*yy-110./(1+np.exp(np.clip(distance*20/.8, -100, 100)))
    contours, _ = cv.findContours((distance<=0).astype(np.uint8)*255, cv.RETR_EXTERNAL, cv.CHAIN_APPROX_SIMPLE)
    return gray.astype(np.uint8), [sem.describe_contour(c) for c in contours]


def analytic_trenches():
    yy, xx = np.mgrid[:400, :720]
    sigmoid = lambda z: 1/(1+np.exp(np.clip(-z, -100, 100)))
    body = np.zeros(xx.shape, float)
    for x in range(40, 720, 80):
        body += sigmoid((xx-(x-10))/.8)*sigmoid(((x+10)-xx)/.8)*(sigmoid((170-yy)/.8)+sigmoid((yy-230)/.8))
    return (170.+.03*xx+.02*yy-110.*body).astype(np.uint8)


class EdgeAccuracyTests(unittest.TestCase):
    def test_measurement_reader_preserves_8bit_and_unclipped_16bit_signal(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)/'signal.tif'
            original = np.arange(256, dtype=np.uint8).reshape(16, 16)
            sem._imwrite_unicode_v19(path, original)
            np.testing.assert_array_equal(sem._v31_read_measurement_gray(path), original)
            high_depth = original.astype(np.uint16)*173
            sem._imwrite_unicode_v19(path, high_depth)
            decoded = sem._v31_read_measurement_gray(path)
            self.assertLess(decoded.min(), 0.)
            self.assertGreater(decoded.max(), 255.)
            # Every step remains linear, including values beyond the percentiles.
            np.testing.assert_allclose(np.diff(decoded.ravel()), np.diff(decoded.ravel()).mean(), atol=5e-5)

    def test_half_level_accuracy_and_shape_are_preserved(self):
        for power in (2, 6):
            for angle in (0., 23.):
                with self.subTest(power=power, angle=angle):
                    gray, coarse, xa, ya = analytic_shape(power, angle)
                    contour, audit = sem._v31_supported_contour(gray, coarse)
                    self.assertIsNotNone(contour, audit)
                    self.assertFalse(audit['edge_template_forced'])
                    self.assertGreater(audit['edge_support_fraction'], .95)
                    points = contour[:, 0]
                    self.assertAlmostEqual(np.ptp(points @ xa), 44., delta=.7)
                    self.assertAlmostEqual(np.ptp(points @ ya), 64., delta=.7)
                    if power == 6:
                        # A rounded rectangle must not collapse into an ellipse.
                        self.assertGreater(cv.contourArea(contour)/(np.pi*22*32), 1.15)

    def test_uniform_and_tilted_backgrounds_do_not_create_edges(self):
        _, coarse, _, _ = analytic_shape()
        yy, xx = np.mgrid[:160, :180]
        for gray in (np.full(xx.shape, 120.), 80.+.24*xx+.12*yy):
            contour, audit = sem._v31_supported_contour(gray, coarse)
            self.assertIsNone(contour)
            self.assertEqual(audit['edge_status'], 'INSUFFICIENT_ORIGINAL_IMAGE_SUPPORT')

    def test_enhanced_locator_cannot_replace_original_measurement_pixels(self):
        enhanced, coarse, _, _ = analytic_shape()
        desc = sem.describe_contour(coarse)
        original = np.full(enhanced.shape, 120., np.float32)
        with patch.object(sem, '_V31_MEASUREMENT_GRAY', original):
            contour, _ = sem._v31_refine_desc(enhanced, desc, 'via')
        self.assertIsNone(contour)

    def test_original_gray_context_restored_even_after_failure(self):
        old = sem._V31_MEASUREMENT_GRAY
        image = np.ones((5, 5), np.uint8)
        def fail(*args):
            self.assertIs(sem._V31_MEASUREMENT_GRAY, image)
            raise RuntimeError('measurement failure')
        with patch.object(sem, '_v31_read_measurement_gray', return_value=image), \
             patch.object(sem, '_V31_PROCESS_BASE', side_effect=fail):
            with self.assertRaises(RuntimeError):
                sem._process_one_image_v19(Path('image.png'), 'before', '1', Path('ann'))
        self.assertIs(sem._V31_MEASUREMENT_GRAY, old)


class FusionPipelineTests(unittest.TestCase):
    def test_file_to_measurement_uses_raw_gray_with_v27_locator(self):
        with tempfile.TemporaryDirectory() as temporary, contextlib.redirect_stdout(io.StringIO()):
            root = Path(temporary)
            image = root/'slot.png'
            sem._imwrite_unicode_v19(image, gradient_slots(4.))
            image.with_suffix('.txt').write_text('PixelSize=0.8')
            rows, status, _ = sem._v31_process_v25_assigned(image, 'before', '007', 2, 'slot', 'slot210', root/'ann')
            self.assertEqual(len(rows), 3)
            self.assertTrue(Path(status['annotated_path']).is_file())
            np.testing.assert_allclose([r['x_width_nm'] for r in rows], 17*.8, atol=1.)
            np.testing.assert_allclose([r['y_height_nm'] for r in rows], 21*.8, atol=1.)

    def test_real_segmentation_candidates_do_not_bias_weak_slot_size(self):
        for contrast in (4., 8., 40.):
            with self.subTest(contrast=contrast), contextlib.redirect_stdout(io.StringIO()):
                gray = gradient_slots(contrast)
                initial, _ = sem.extract_descriptors_with_rejections(sem.segment_features(gray))
                pool, _ = sem.build_candidate_pool_v12(gray, initial, 'slot', 'weak.png')
                rows, _, _, _, _ = sem._select_and_measure_bottom_slots_v19(
                    gray, pool, .8, 'weak.png', 'before', '1', 'slot210')
                self.assertEqual(len(rows), 3)
                np.testing.assert_allclose([r['x_width_nm'] for r in rows], 17*.8, atol=1.)
                np.testing.assert_allclose([r['y_height_nm'] for r in rows], 21*.8, atol=1.)

    def test_full_slot_pipeline_weak_bottom_displacement_and_polarity(self):
        for contrast, bottom_y, bottom_contrast, polarity in (
                (8., 325, None, 'dark'), (4., 325, None, 'dark'),
                (8., 345, None, 'dark'), (40., 325, 5., 'dark'), (8., 325, None, 'bright')):
            with self.subTest(contrast=contrast, bottom_y=bottom_y, bottom_contrast=bottom_contrast, polarity=polarity), \
                 patch.object(sem, 'FEATURE_POLARITY', polarity), contextlib.redirect_stdout(io.StringIO()):
                gray = gradient_slots(contrast, bottom_y, bottom_contrast)
                if polarity == 'bright':
                    gray = 255-gray
                rows, selected, bottom, _, meta = sem._select_and_measure_bottom_slots_v19(
                    gray, [], .8, 'gradient.png', 'before', '007', 'slot210')
                self.assertEqual(len(rows), 3)
                self.assertEqual(len(selected), 12)
                self.assertEqual(len(bottom), 3)
                self.assertEqual(meta['v23_complete_slot_columns'], 3)
                # cv.rectangle is inclusive: true half-gray dimensions are 17 x 21 px.
                np.testing.assert_allclose([r['x_width_nm'] for r in rows], 17*.8, atol=.8)
                np.testing.assert_allclose([r['y_height_nm'] for r in rows], 21*.8, atol=.8)
                for row, desc in zip(rows, bottom):
                    self.assertEqual(row['slot_row_id_from_v18'], 4)
                    self.assertEqual(row['edge_signal_source'], 'original_image_gray')
                    self.assertAlmostEqual(row['x_width_nm'], row['x_width_px']*.8)
                    self.assertAlmostEqual(row['x_width_px'], sem._projection_span_v19(desc.contour, desc._v22_slot_axis['x_axis']))

    def test_valid_v27_columns_are_not_replaced_by_v30(self):
        gray = gradient_slots(40.)
        pool = sem._v31_slot_flatfield_candidates(gray, [])
        with contextlib.redirect_stdout(io.StringIO()), \
             patch.object(sem, '_v31_slot_flatfield_candidates', side_effect=AssertionError('unnecessary V30 replacement')):
            rows, _, _, _, _ = sem._select_and_measure_bottom_slots_v19(gray, pool, .8, 'good.png', 'before', '1', 'slot210')
        self.assertEqual(len(rows), 3)
        self.assertEqual({r['v31_locator'] for r in rows}, {'V27'})

    def test_missing_bottom_does_not_produce_a_complete_column(self):
        gray = gradient_slots(40.)
        pool = sem._v31_slot_flatfield_candidates(gray, [])
        with contextlib.redirect_stdout(io.StringIO()):
            original = sem._V31_SLOT_BASE(gray, pool, .8, 'missing.png', 'before', '1', 'slot210')
        yy, xx = np.mgrid[:400, :360]
        background = (80.+.24*xx+.12*yy).astype(np.uint8)
        gray[305:345, 160:200] = background[305:345, 160:200]
        columns, rejects, _ = sem._v31_valid_slot_columns(gray, original, 'V27')
        self.assertEqual(len(columns), 2)
        self.assertTrue(any(r.get('edge_status') == 'INSUFFICIENT_ORIGINAL_IMAGE_SUPPORT' for r in rejects))

    def test_via_lattice_keeps_twenty_actual_gray_boundaries(self):
        gray, pool = analytic_vias()
        with contextlib.redirect_stdout(io.StringIO()):
            rows, selected, _, bundle = sem._select_and_measure_vias_v19(gray, pool, 1., 'via.png', 'before', '1', 'via40')
        self.assertEqual(len(rows), 20)
        np.testing.assert_allclose([r['x_width_nm'] for r in rows], 40., atol=.7)
        np.testing.assert_allclose([r['y_height_nm'] for r in rows], 40., atol=.7)
        self.assertEqual(bundle['meta']['via_v24_final_count'], 20)
        self.assertEqual(sum(r['is_axis_bottom_row'] for r in rows), 5)
        self.assertEqual(sem._annotate_vias_v19(gray, selected, rows, bundle).shape[:2], gray.shape)

    def test_trench_mean_rectangle_width_and_projected_gap(self):
        gray = analytic_trenches()
        rows, pairs, _, _ = sem._select_and_measure_trenches_v19(gray, [], 1., 'trench.png', 'before', '1', 'trench160')
        self.assertEqual(len(rows), 9)
        np.testing.assert_allclose([r['x_width_nm'] for r in rows], 20., atol=.7)
        np.testing.assert_allclose([r['tip_gap_y_nm'] for r in rows], 60., atol=.7)
        for row, pair in zip(rows, pairs):
            self.assertAlmostEqual(row['x_width_nm'], .5*(row['top_trench_width_nm']+row['bottom_trench_width_nm']))
            self.assertAlmostEqual(pair['top_width_px_v22'], row['top_trench_width_nm'])
        self.assertEqual(sem._annotate_trenches_v19(gray, pairs, rows).shape[:2], gray.shape)


class TimeoutTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def image(self, name='slot.png'):
        image = self.root / name
        sem._imwrite_unicode_v19(image, gradient_slots())
        image.with_suffix('.txt').write_text('PixelSize=0.8\n')
        return image

    def test_hard_timeout_kills_worker_and_discards_partial_output(self):
        fake_worker = self.root / 'blocked_worker.py'
        marker = self.root / 'started'
        fake_worker.write_text(
            'import json, sys, time\nfrom pathlib import Path\n'
            'request=json.loads(Path(sys.argv[2]).read_text())\n'
            f'Path({str(marker)!r}).write_text("started")\n'
            'output=Path(request["ann_dir"]); output.mkdir()\n'
            '(output/"partial.png").write_bytes(b"incomplete")\n'
            'while True: time.sleep(1)\n')
        processes = []
        original_popen = sem._subprocess_v31.Popen
        def start(*args, **kwargs):
            child = original_popen(*args, **kwargs)
            processes.append(child)
            return child
        started = time.monotonic()
        with patch.object(sem, '__file__', str(fake_worker)), \
                patch.object(sem._subprocess_v31, 'Popen', side_effect=start):
            with self.assertRaises(sem.V31ImageTimeout):
                sem._v31_process_image_with_timeout(
                    self.root / 'input.png', 'before', '1', 2, 'slot', 'slot210',
                    self.root / 'annotations', timeout_seconds=.5)
        self.assertLess(time.monotonic() - started, 4.)
        self.assertTrue(marker.exists())
        self.assertIsNotNone(processes[0].poll())
        self.assertFalse((self.root / 'annotations/partial.png').exists())

    def test_real_standalone_worker_publishes_valid_result(self):
        rows, status, _, key, trials = sem._v31_process_image_with_timeout(
            self.image(), 'before', '007', 2, 'slot', 'slot210',
            self.root / 'annotations', timeout_seconds=30.)
        self.assertEqual(len(rows), 3)
        self.assertEqual(key, 'slot210')
        self.assertEqual(trials, [])
        self.assertEqual(status['source_algorithm'], 'V31_V27_Yakun_fusion')
        self.assertTrue(Path(status['annotated_path']).is_file())
        self.assertTrue(Path(status['annotated_path']).is_relative_to(self.root / 'annotations'))
        self.assertLess(status['processing_elapsed_seconds'], 30.)

    def test_both_via_trials_share_job_and_winner_is_not_reprocessed(self):
        def measure(image, stage, condition, region, pattern, key, ann_dir):
            rows = [dict(x_width_nm=40., y_height_nm=40.) for _ in range(4)]
            return rows, {'measure_array_occupancy_seed': .9 if key == 'via60' else .5}, []
        job = dict(image_path='test.png', ann_dir=str(self.root / 'ann'), stage='before',
                   condition='1', region=2, public_pattern='via', v25_pattern_key=None)
        with patch.object(sem, '_v31_process_v25_assigned', side_effect=measure) as process:
            _, _, _, key, trials = sem._v31_measure_image_job(job)
        self.assertEqual(key, 'via60')
        self.assertEqual(process.call_count, 2)
        self.assertEqual([t['selected'] for t in trials], [False, True])

    def test_timeout_limit_cannot_exceed_five_minutes(self):
        for limit in (0, -1, 301, float('inf'), float('nan')):
            with self.subTest(limit=limit), self.assertRaises(ValueError):
                sem._v31_process_image_with_timeout(
                    Path('unused.png'), 'before', '1', 2, 'slot', 'slot210', self.root, limit)

    def test_batch_exports_timeout_and_continues_later_images(self):
        for stage in ('before', 'after'):
            for condition in ('007', '10'):
                folder = self.root / stage / condition
                folder.mkdir(parents=True)
                for index in (1, 2):
                    (folder / f'{index}.png').touch()
                    (folder / f'{index}.txt').write_text('PixelSize=0.8')
        calls = []
        def raw(image_path, stage, condition, ann_dir):
            row = dict(stage=stage, stage_zh=sem.V31_STAGE_ZH[stage], condition=condition,
                       condition_label=f'Condition {condition}', image=image_path.name,
                       pattern=sem._classify_pattern_key_v19(image_path), object_id=1,
                       x_width_nm=10., y_height_nm=20., tip_gap_y_nm=30.,
                       center_image_x_px=50., center_image_y_px=50.)
            return [row], dict(row, status='OK'), []
        def process(**kw):
            calls.append((kw['condition'], kw['stage'], kw['public_pattern']))
            if len(calls) == 1:
                raise sem.V31ImageTimeout('test image exceeded 300s')
            args = {k: v for k, v in kw.items() if k != 'timeout_seconds'}
            rows, status, rejects = sem._v31_process_v25_assigned(**args)
            return rows, status, rejects, kw['v25_pattern_key'], []
        output = self.root / 'output'
        with patch.object(sem, '_process_one_image_v19', side_effect=raw), \
                patch.object(sem, '_v31_process_image_with_timeout', side_effect=process), \
                patch.object(sem, '_v31_save_pair_plots'), contextlib.redirect_stdout(io.StringIO()):
            paths = sem.main_v31(self.root / 'before', self.root / 'after', output)
        self.assertEqual(len(calls), 8)
        status = sem.pd.read_csv(paths['image_status'], dtype={'condition': str})
        self.assertEqual(status['status'].tolist().count('SKIPPED_TIMEOUT'), 1)
        self.assertEqual(status['status'].tolist().count('OK'), 7)
        errors = sem.pd.read_csv(paths['errors'])
        self.assertEqual(errors['error_type'].tolist(), ['IMAGE_TIMEOUT'])
        settings = json.loads(paths['settings'].read_text())
        self.assertEqual(settings['image_timeout_seconds'], 300.)
        self.assertIn('V31', settings['script_version'])
        self.assertEqual(paths['excel'].name, 'SEM_0830_0831_before_after_V31_results.xlsx')
        pairs = sem.pd.read_csv(paths['image_pair_comparison'], dtype={'condition': str})
        partial = pairs[(pairs['condition'] == '007') & (pairs['pattern'] == 'trench')]
        self.assertTrue(partial['delta_after_minus_before_nm'].isna().all())

    def test_auto_via_timeout_skips_dependent_after_and_continues_next_region(self):
        for stage in ('before', 'after'):
            folder = self.root / stage / '1'
            folder.mkdir(parents=True)
            for index in range(1, 5):
                (folder / f'{index}.png').touch()
                (folder / f'{index}.txt').write_text('PixelSize=0.8')
        calls = []
        def process(**kw):
            calls.append((kw['region'], kw['public_pattern'], kw['stage']))
            if kw['public_pattern'] == 'via':
                self.assertIsNone(kw['v25_pattern_key'])
                raise sem.V31ImageTimeout('via trials exhausted the shared image deadline')
            status = dict(stage=kw['stage'], condition='1', region=kw['region'],
                          pattern=kw['public_pattern'], image=kw['image_path'].name, status='OK')
            return [], status, [], kw['v25_pattern_key'], []
        with patch.object(sem, '_v31_process_image_with_timeout', side_effect=process), \
                patch.object(sem, '_v31_save_pair_plots'), contextlib.redirect_stdout(io.StringIO()):
            paths = sem.main_v31(self.root / 'before', self.root / 'after', self.root / 'output')
        self.assertEqual(len(calls), 7)
        self.assertEqual(calls[-2:], [(3, 'trench', 'before'), (3, 'trench', 'after')])
        status = sem.pd.read_csv(paths['image_status'])
        via = status[status['pattern'] == 'via']
        self.assertEqual(via['status'].tolist(), ['SKIPPED_TIMEOUT', 'SKIPPED_VIA_MODEL_UNAVAILABLE'])
        errors = sem.pd.read_csv(paths['errors'])
        self.assertEqual(errors['error_type'].tolist(), ['IMAGE_TIMEOUT', 'VIA_MODEL_UNAVAILABLE'])


if __name__ == '__main__':
    unittest.main()
