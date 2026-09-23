"""Regression tests for variable-size before/after image inventories."""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'SEM_0830_0831_compare_V27_bundle'))
import sem_before_after_compare_0831_V27_V25only as sem


class InventoryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.before = self.root / 'before'
        self.after = self.root / 'after'

    def make_images(self, root, condition, count, prefix='image', names=None):
        folder = root / condition
        folder.mkdir(parents=True, exist_ok=True)
        if names is None:
            names = [f'{prefix}_{i}.tif' for i in range(1, count + 1)]
        for name in reversed(names):
            image = folder / name
            image.touch()
            image.with_suffix('.txt').write_text('PixelSize=0.8\n', encoding='utf-8')
        return folder

    def test_equal_counts_have_no_fixed_limit_or_multiple_of_three_requirement(self):
        for count in (0, 1, 2, 3, 8, 9, 10, 11, 12, 37):
            with self.subTest(count=count):
                condition = str(100 + count)
                self.make_images(self.before, condition, count)
                self.make_images(self.after, condition, count, prefix='different')
                inv, errors = sem._v27_build_inventory(self.before, self.after, [condition])
                self.assertTrue(errors.empty, errors.to_dict('records'))
                self.assertEqual(len(inv), count * 2)
                if count:
                    self.assertEqual(inv['pair_key'].nunique(), count)
                    self.assertTrue((inv.groupby('pair_key')['stage'].nunique() == 2).all())
                    for stage in ('before', 'after'):
                        rows = inv[inv['stage'] == stage].sort_values('sequence_index')
                        self.assertEqual(rows['sequence_index'].tolist(), list(range(1, count + 1)))

    def test_original_nine_positions_are_preserved_and_extended(self):
        self.make_images(self.before, '10', 11)
        self.make_images(self.after, '10', 11)
        inv, errors = sem._v27_build_inventory(self.before, self.after, ['10'])
        before = inv[inv['stage'] == 'before'].sort_values('sequence_index')
        self.assertEqual(list(zip(before['region'], before['pattern'])), [
            (2, 'trench'), (2, 'slot'), (2, 'via'),
            (3, 'trench'), (3, 'slot'), (3, 'via'),
            (4, 'trench'), (4, 'slot'), (4, 'via'),
            (5, 'trench'), (5, 'slot'),
        ])

    def test_natural_order_mixed_tif_extensions_and_different_filenames(self):
        before_names = ['scan_1.TIF', 'scan_2.tiff', 'scan_10.TIFF']
        after_names = ['result_4.tif', 'result_8.TIFF', 'result_20.tiff']
        before_folder = self.make_images(self.before, '10', 3, names=before_names)
        self.make_images(self.after, '10', 3, names=after_names)
        (before_folder / 'preview.png').touch()
        (before_folder / 'notes.txt').touch()
        nested = before_folder / 'nested'
        nested.mkdir()
        (nested / 'extra.tif').touch()
        inv, errors = sem._v27_build_inventory(self.before, self.after, ['10'])
        self.assertTrue(errors.empty)
        self.assertEqual(len(inv), 6)
        for stage, expected in [('before', before_names), ('after', after_names)]:
            self.assertEqual(inv[inv['stage'] == stage].sort_values('sequence_index')['image'].tolist(), expected)

    def test_mismatches_skip_both_stages_without_affecting_valid_conditions(self):
        conditions = []
        for i, (bcount, acount) in enumerate([(10, 9), (9, 10), (0, 1), (1, 0)]):
            condition = str(10 + i)
            conditions.append(condition)
            self.make_images(self.before, condition, bcount)
            self.make_images(self.after, condition, acount)
        self.make_images(self.before, '14', 12)
        self.make_images(self.after, '14', 12)
        inv, errors = sem._v27_build_inventory(self.before, self.after, conditions + ['14'])
        self.assertEqual(set(inv['condition']), {'14'})
        self.assertEqual(len(inv), 24)
        self.assertEqual(errors['error_type'].tolist(), ['IMAGE_COUNT_MISMATCH'] * 4)
        self.assertEqual(list(zip(errors['before_count'], errors['after_count'])), [(10, 9), (9, 10), (0, 1), (1, 0)])

    def test_missing_folder_does_not_leave_one_sided_inventory(self):
        self.make_images(self.before, '10', 10)
        self.make_images(self.after, '11', 0)
        inv, errors = sem._v27_build_inventory(self.before, self.after, ['10', '11'])
        self.assertTrue(inv.empty)
        self.assertEqual(errors['error_type'].tolist(), ['MISSING_CONDITION_FOLDER'] * 2)

    def test_metadata_failure_does_not_change_the_image_count(self):
        folder = self.make_images(self.before, '10', 10)
        self.make_images(self.after, '10', 10)
        (folder / 'image_10.txt').unlink()
        inv, errors = sem._v27_build_inventory(self.before, self.after, ['10'])
        self.assertEqual(len(inv), 20)
        self.assertEqual(errors['error_type'].tolist(), ['PIXELSIZE_PREFLIGHT_ERROR'])
        self.assertEqual(errors.iloc[0]['image'], 'image_10.tif')

    def test_all_invalid_inputs_still_export_count_diagnostics(self):
        self.make_images(self.before, '10', 10)
        self.make_images(self.after, '10', 9)
        output = self.root / 'output'
        with contextlib.redirect_stdout(io.StringIO()), patch.object(sem, '_process_one_image_v19') as process:
            with self.assertRaisesRegex(RuntimeError, '照片数量是否一致'):
                sem.main_v27(self.before, self.after, output, conditions=['10'])
            process.assert_not_called()
        errors = sem.pd.read_csv(output / 'processing_errors.csv')
        self.assertEqual(errors.iloc[0]['error_type'], 'IMAGE_COUNT_MISMATCH')
        self.assertEqual(errors.iloc[0]['before_count'], 10)
        self.assertEqual(errors.iloc[0]['after_count'], 9)

    def test_pipeline_exports_all_pairs_beyond_nine(self):
        for condition, bcount, acount in [('10', 13, 13), ('11', 8, 9), ('12', 2, 2)]:
            self.make_images(self.before, condition, bcount)
            self.make_images(self.after, condition, acount, prefix='result')
        output = self.root / 'output'

        def fake_measurement(image_path, stage, condition, ann_dir):
            # Only the image detector is replaced; inventory, V27 wrapper,
            # pairing, summaries, CSV and Excel generation all run normally.
            row = {
                'stage': stage, 'stage_zh': sem.V27_STAGE_ZH[stage],
                'condition': condition, 'condition_label': f'Condition {condition}',
                'image': image_path.name,
                'pattern': sem._classify_pattern_key_v19(image_path),
                'object_id': 1, 'x_width_nm': 10.0 + (stage == 'after'),
                'y_height_nm': 20.0 + (stage == 'after'),
                'tip_gap_y_nm': 30.0 + (stage == 'after'),
                'center_image_x_px': 50.0, 'center_image_y_px': 50.0,
            }
            return [row], {**row, 'status': 'OK'}, []

        with contextlib.redirect_stdout(io.StringIO()), \
                patch.object(sem, '_process_one_image_v19', side_effect=fake_measurement) as process, \
                patch.object(sem, '_v27_save_pair_plots'):
            outputs = sem.main_v27(self.before, self.after, output, conditions=['10', '11', '12'], via_pattern_mode='via40')
        self.assertEqual(process.call_count, 30)
        inventory = sem.pd.read_csv(outputs['inventory'])
        statuses = sem.pd.read_csv(outputs['image_status'])
        pairs = sem.pd.read_csv(outputs['image_pair_comparison'])
        objects = sem.pd.read_csv(outputs['object_comparison'])
        errors = sem.pd.read_csv(outputs['errors'])
        self.assertEqual(len(inventory), 30)
        self.assertEqual(len(statuses), 30)
        self.assertEqual(len(pairs), 30)
        self.assertEqual(len(objects), 30)
        self.assertEqual(pairs['pair_key'].nunique(), 15)
        self.assertTrue((pairs['delta_after_minus_before_nm'] == 1.0).all())
        self.assertTrue((objects['match_status'] == 'MATCHED').all())
        last_pair = pairs[pairs['pair_key'] == 'C10_R6_trench']
        self.assertEqual(len(last_pair), 2)
        self.assertEqual(set(last_pair['before_image']), {'image_13.tif'})
        self.assertEqual(set(last_pair['after_image']), {'result_13.tif'})
        self.assertEqual(errors['error_type'].tolist(), ['IMAGE_COUNT_MISMATCH'])
        settings = json.loads(outputs['settings'].read_text(encoding='utf-8'))
        self.assertIsNone(settings['required_tif_count_per_condition_stage'])
        self.assertEqual(settings['regions'], [2, 3, 4, 5, 6])
        from openpyxl import load_workbook
        workbook = load_workbook(outputs['excel'], read_only=True)
        self.addCleanup(workbook.close)
        self.assertEqual(workbook['inventory'].max_row, 31)
        self.assertEqual(workbook['image_pair_compare'].max_row, 31)


if __name__ == '__main__':
    unittest.main()
