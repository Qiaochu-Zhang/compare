"""V30 recovery accuracy, hard timeouts, and batch continuation."""
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
sys.path.insert(0, str(ROOT / 'SEM_0830_0831_compare_V30_bundle'))
import sem_before_after_compare_0831_V30 as sem


def gradient_slots(contrast=8., bottom_y=325, bottom_contrast=None):
    yy, xx = sem.np.mgrid[:400, :360]
    gray = 80. + .24 * xx + .12 * yy
    mask = sem.np.zeros(gray.shape, sem.np.uint8)
    for row, y in enumerate((55, 145, 235, bottom_y)):
        height = (55, 43, 31, 21)[row]
        for x in (60, 180, 300):
            sem.cv2.rectangle(mask, (x - 8, y - height // 2), (x + 8, y + height // 2), 1, -1)
    amplitudes = sem.np.full(gray.shape, contrast, dtype=sem.np.float32)
    if bottom_contrast is not None:
        amplitudes[bottom_y - 30:bottom_y + 31] = bottom_contrast
    gray -= sem.cv2.GaussianBlur(mask.astype(sem.np.float32) * amplitudes, (0, 0), .8)
    gray += sem.np.random.default_rng(30).normal(0, .45, gray.shape)
    return sem.np.clip(gray, 0, 255).astype(sem.np.uint8)


class RecoveryTests(unittest.TestCase):
    def test_recovers_visible_slot_missed_by_legacy_under_gradient(self):
        image = gradient_slots(8.)
        args = (image, (180., 325.), 120., 90., 8., 21., (30., 2000.), 0)
        self.assertIsNone(sem._recover_local_slot_v18_legacy(*args))
        result = sem._recover_local_slot_v18_deep(*args)
        self.assertIsNotNone(result)
        x, y, width, height = sem.cv2.boundingRect(result.contour)
        self.assertAlmostEqual(width, 17, delta=2)
        self.assertAlmostEqual(height, 21, delta=2)
        self.assertAlmostEqual(x + width / 2, 180., delta=2)
        self.assertAlmostEqual(y + height / 2, 325., delta=2)

    def test_weak_partial_contour_is_replaced_by_full_slot(self):
        image = gradient_slots(4.)
        args = (image, (180., 325.), 120., 90., 8., 21., (30., 2000.), 0)
        old = sem._recover_local_slot_v18_legacy(*args)
        new = sem._recover_local_slot_v18_deep(*args)
        self.assertIsNotNone(new)
        self.assertLess(sem.cv2.boundingRect(old.contour)[2], 12)
        self.assertAlmostEqual(sem.cv2.boundingRect(new.contour)[2], 17, delta=2)

    def test_full_pipeline_weak_slots_displaced_row_and_bright_polarity(self):
        for contrast, bottom_y, polarity in ((8., 325, 'dark'), (4., 325, 'dark'),
                                             (8., 345, 'dark'), (8., 325, 'bright')):
            with self.subTest(contrast=contrast, bottom_y=bottom_y, polarity=polarity), \
                    patch.object(sem, 'FEATURE_POLARITY', polarity):
                image = gradient_slots(contrast, bottom_y)
                if polarity == 'bright':
                    image = 255 - image
                with contextlib.redirect_stdout(io.StringIO()):
                    rows, selected, bottom, _, meta = sem._select_and_measure_bottom_slots_v19(
                        image, [], .8, 'gradient.png', 'before', '007', 'slot210')
                self.assertEqual(len(rows), 3)
                self.assertEqual(len(selected), 12)
                self.assertEqual(len(bottom), 3)
                self.assertEqual(meta['v23_complete_slot_columns'], 3)
                self.assertGreater(meta['slot_v30_flatfield_candidates'], 0)
                sem.np.testing.assert_allclose([r['x_width_nm'] for r in rows], 12.8, atol=1.6)
                sem.np.testing.assert_allclose([r['y_height_nm'] for r in rows], 16., atol=1.6)
                sem.np.testing.assert_allclose([r['center_image_y_px'] for r in rows], bottom_y, atol=3.)

    def test_background_gradient_is_not_accepted_as_a_slot(self):
        yy, xx = sem.np.mgrid[:400, :360]
        image = (80 + .24 * xx + .12 * yy).astype(sem.np.uint8)
        self.assertEqual(sem._v30_slot_flatfield_candidates(image, []), [])
        result = sem._v30_recover_flatfield_slot(image, (180., 325.), 120., 90., 17., 21., (30., 2000.))
        self.assertIsNone(result)

    def test_bottom_row_is_recovered_when_upper_rows_are_much_darker(self):
        image = gradient_slots(contrast=40., bottom_contrast=5.)
        with contextlib.redirect_stdout(io.StringIO()) as log:
            rows, _, _, _, meta = sem._select_and_measure_bottom_slots_v19(
                image, [], .8, 'weak_bottom.png', 'before', '1', 'slot210')
        self.assertEqual(len(rows), 3)
        self.assertGreaterEqual(meta['array_recovered_count'], 3)
        self.assertIn('row=4', log.getvalue())
        sem.np.testing.assert_allclose([r['x_width_nm'] for r in rows], 12.8, atol=1.6)


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
        original_popen = sem._subprocess_v30.Popen
        def start(*args, **kwargs):
            child = original_popen(*args, **kwargs)
            processes.append(child)
            return child
        started = time.monotonic()
        with patch.object(sem, '__file__', str(fake_worker)), \
                patch.object(sem._subprocess_v30, 'Popen', side_effect=start):
            with self.assertRaises(sem.V30ImageTimeout):
                sem._v30_process_image_with_timeout(
                    self.root / 'input.png', 'before', '1', 2, 'slot', 'slot210',
                    self.root / 'annotations', timeout_seconds=.5)
        self.assertLess(time.monotonic() - started, 4.)
        self.assertTrue(marker.exists())
        self.assertIsNotNone(processes[0].poll())
        self.assertFalse((self.root / 'annotations/partial.png').exists())

    def test_real_standalone_worker_publishes_valid_result(self):
        rows, status, _, key, trials = sem._v30_process_image_with_timeout(
            self.image(), 'before', '007', 2, 'slot', 'slot210',
            self.root / 'annotations', timeout_seconds=30.)
        self.assertEqual(len(rows), 3)
        self.assertEqual(key, 'slot210')
        self.assertEqual(trials, [])
        self.assertEqual(status['source_algorithm'], 'V30_slot_recovery')
        self.assertTrue(Path(status['annotated_path']).is_file())
        self.assertTrue(Path(status['annotated_path']).is_relative_to(self.root / 'annotations'))
        self.assertLess(status['processing_elapsed_seconds'], 30.)

    def test_both_via_trials_share_job_and_winner_is_not_reprocessed(self):
        def measure(image, stage, condition, region, pattern, key, ann_dir):
            rows = [dict(x_width_nm=40., y_height_nm=40.) for _ in range(4)]
            return rows, {'measure_array_occupancy_seed': .9 if key == 'via60' else .5}, []
        job = dict(image_path='test.png', ann_dir=str(self.root / 'ann'), stage='before',
                   condition='1', region=2, public_pattern='via', v25_pattern_key=None)
        with patch.object(sem, '_v30_process_v25_assigned', side_effect=measure) as process:
            _, _, _, key, trials = sem._v30_measure_image_job(job)
        self.assertEqual(key, 'via60')
        self.assertEqual(process.call_count, 2)
        self.assertEqual([t['selected'] for t in trials], [False, True])

    def test_timeout_limit_cannot_exceed_five_minutes(self):
        for limit in (0, -1, 301, float('inf'), float('nan')):
            with self.subTest(limit=limit), self.assertRaises(ValueError):
                sem._v30_process_image_with_timeout(
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
            row = dict(stage=stage, stage_zh=sem.V30_STAGE_ZH[stage], condition=condition,
                       condition_label=f'Condition {condition}', image=image_path.name,
                       pattern=sem._classify_pattern_key_v19(image_path), object_id=1,
                       x_width_nm=10., y_height_nm=20., tip_gap_y_nm=30.,
                       center_image_x_px=50., center_image_y_px=50.)
            return [row], dict(row, status='OK'), []
        def process(**kw):
            calls.append((kw['condition'], kw['stage'], kw['public_pattern']))
            if len(calls) == 1:
                raise sem.V30ImageTimeout('test image exceeded 300s')
            args = {k: v for k, v in kw.items() if k != 'timeout_seconds'}
            rows, status, rejects = sem._v30_process_v25_assigned(**args)
            return rows, status, rejects, kw['v25_pattern_key'], []
        output = self.root / 'output'
        with patch.object(sem, '_process_one_image_v19', side_effect=raw), \
                patch.object(sem, '_v30_process_image_with_timeout', side_effect=process), \
                patch.object(sem, '_v30_save_pair_plots'), contextlib.redirect_stdout(io.StringIO()):
            paths = sem.main_v30(self.root / 'before', self.root / 'after', output)
        self.assertEqual(len(calls), 8)
        status = sem.pd.read_csv(paths['image_status'], dtype={'condition': str})
        self.assertEqual(status['status'].tolist().count('SKIPPED_TIMEOUT'), 1)
        self.assertEqual(status['status'].tolist().count('OK'), 7)
        errors = sem.pd.read_csv(paths['errors'])
        self.assertEqual(errors['error_type'].tolist(), ['IMAGE_TIMEOUT'])
        settings = json.loads(paths['settings'].read_text())
        self.assertEqual(settings['image_timeout_seconds'], 300.)
        self.assertIn('V30', settings['script_version'])
        self.assertEqual(paths['excel'].name, 'SEM_0830_0831_before_after_V30_results.xlsx')
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
                raise sem.V30ImageTimeout('via trials exhausted the shared image deadline')
            status = dict(stage=kw['stage'], condition='1', region=kw['region'],
                          pattern=kw['public_pattern'], image=kw['image_path'].name, status='OK')
            return [], status, [], kw['v25_pattern_key'], []
        with patch.object(sem, '_v30_process_image_with_timeout', side_effect=process), \
                patch.object(sem, '_v30_save_pair_plots'), contextlib.redirect_stdout(io.StringIO()):
            paths = sem.main_v30(self.root / 'before', self.root / 'after', self.root / 'output')
        self.assertEqual(len(calls), 7)
        self.assertEqual(calls[-2:], [(3, 'trench', 'before'), (3, 'trench', 'after')])
        status = sem.pd.read_csv(paths['image_status'])
        via = status[status['pattern'] == 'via']
        self.assertEqual(via['status'].tolist(), ['SKIPPED_TIMEOUT', 'SKIPPED_VIA_MODEL_UNAVAILABLE'])
        errors = sem.pd.read_csv(paths['errors'])
        self.assertEqual(errors['error_type'].tolist(), ['IMAGE_TIMEOUT', 'VIA_MODEL_UNAVAILABLE'])


if __name__ == '__main__':
    unittest.main()
