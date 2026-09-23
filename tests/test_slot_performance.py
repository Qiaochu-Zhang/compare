"""Regression coverage for V29's silent slot row-fitting slowdown."""
import contextlib
import io
import json
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'SEM_0830_0831_compare_V29_bundle'))
import sem_before_after_compare_0831_V29_V25only as sem


def slot_image():
    gray = sem.np.full((400, 360), 220, sem.np.uint8)
    for row, y in enumerate((55, 145, 235, 325)):
        height = (55, 43, 31, 21)[row]
        for x in (60, 180, 300):
            sem.cv2.rectangle(gray, (x - 8, y - height // 2),
                              (x + 8, y + height // 2), 30, -1)
    return gray


class SlotPerformanceTests(unittest.TestCase):
    def test_row_fit_matches_frozen_exhaustive_v25_results(self):
        # Produced by the original exhaustive fitter before the V29 fix. Includes
        # missing rows, duplicate centers, boundary ties and randomized noise, at
        # three tolerances (including overlapping support intervals).
        cases = json.loads((ROOT / 'tests/data/slot_row_fit_v29_legacy.json').read_text())
        for case in cases:
            with self.subTest(name=case['name'], fraction=case['residual_fraction']):
                candidates = [SimpleNamespace(center=sem.np.array([0., y])) for y in case['ys']]
                with patch.object(sem, 'V17_SLOT_ROW_MODEL_RESIDUAL_FRAC', case['residual_fraction']):
                    actual = sem._fit_four_row_centers_v17(candidates, case['image_h'])
                expected = case['expected']
                if expected is None:
                    self.assertIsNone(actual)
                else:
                    self.assertEqual(set(actual), set(expected))
                    for key in expected:
                        sem.np.testing.assert_equal(actual[key], expected[key])

    def test_dense_noisy_four_rows_keep_all_seed_support(self):
        rng = sem.np.random.default_rng(29)
        ys = sem.np.repeat([100., 300., 500., 700.], 45) + rng.uniform(-2., 2., 180)
        candidates = [SimpleNamespace(center=sem.np.array([0., y])) for y in ys]
        result = sem._fit_four_row_centers_v17(candidates, 800)
        self.assertEqual(result['seed_inliers'], 180)
        self.assertEqual(result['occupied_seed_rows'], 4)
        sem.np.testing.assert_allclose(result['row_centers'], [100., 300., 500., 700.], atol=1.)

    def test_lowgray_cap_stops_mask_generation_and_keeps_candidate_order(self):
        gray = slot_image()
        for polarity in ('dark', 'bright'):
            with self.subTest(polarity=polarity), patch.object(sem, 'FEATURE_POLARITY', polarity):
                image = gray if polarity == 'dark' else 255 - gray
                complete = sem._slot_lowgray_model_proposals_v18(image)
                self.assertGreater(len(complete), 1)
                with patch.object(sem, 'V18_SLOT_MODEL_LOWGRAY_MAX_PROPOSALS', 1), \
                        patch.object(sem.np, 'percentile', wraps=sem.np.percentile) as percentile:
                    capped = sem._slot_lowgray_model_proposals_v18(image)
                self.assertEqual(len(capped), 1)
                sem.np.testing.assert_array_equal(capped[0].contour, complete[0].contour)
                self.assertLess(percentile.call_count, 18)

    def test_actual_bottom_row_measurement_and_sparse_seed_recovery(self):
        gray = slot_image()
        pool, _ = sem.extract_descriptors_with_rejections(sem.segment_features(gray))
        for sparse in (False, True):
            with self.subTest(sparse=sparse):
                log = io.StringIO()
                with contextlib.redirect_stdout(log):
                    rows, selected, bottom, _, meta = sem._select_and_measure_bottom_slots_v19(
                        gray, [] if sparse else pool.copy(), 0.8,
                        '210_synthetic.png', 'before', '1', 'slot210')
                self.assertEqual(len(rows), 3)
                self.assertEqual(len(selected), 12)
                self.assertEqual(len(bottom), 3)
                self.assertEqual(meta['v23_complete_slot_columns'], 3)
                self.assertTrue(all(row['slot_row_id_from_v18'] == 4 for row in rows))
                self.assertTrue(all(row['center_image_y_px'] > 300 for row in rows))
                self.assertIn('slot row fit done:', log.getvalue())
                if sparse:
                    self.assertGreater(meta['slot_lowgray_locator_proposals'], 0)
                    self.assertIn('slot low-gray proposals=', log.getvalue())
                else:
                    sem.np.testing.assert_allclose([row['x_width_nm'] for row in rows], 12.8)
                    sem.np.testing.assert_allclose([row['y_height_nm'] for row in rows], 16.)


if __name__ == '__main__':
    unittest.main()
