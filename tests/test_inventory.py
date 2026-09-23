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
sys.path.insert(0, str(ROOT / 'SEM_0830_0831_compare_V29_bundle'))
import sem_before_after_compare_0831_V29_V25only as sem


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

    @staticmethod
    def fake_measurement(image_path, stage, condition, ann_dir):
        # Only the image detector is replaced; inventory, V29 wrapper,
        # pairing, summaries, CSV and Excel generation all run normally.
        row = {
            'stage': stage, 'stage_zh': sem.V29_STAGE_ZH[stage],
            'condition': condition, 'condition_label': f'Condition {condition}',
            'image': image_path.name,
            'pattern': sem._classify_pattern_key_v19(image_path),
            'object_id': 1, 'x_width_nm': 10.0 + (stage == 'after'),
            'y_height_nm': 20.0 + (stage == 'after'),
            'tip_gap_y_nm': 30.0 + (stage == 'after'),
            'center_image_x_px': 50.0, 'center_image_y_px': 50.0,
        }
        return [row], {**row, 'status': 'OK'}, []

    def test_discovery_uses_only_numeric_immediate_after_directories(self):
        expected = ['0', '2', '007', '7', '10', '103', '12345678901234567890']
        for name in reversed(expected):
            self.make_images(self.after, name, 0)
        for name in ('notes', '-1', '1.5', '12_suffix', '²'):
            (self.after / name / '99').mkdir(parents=True)
        (self.after / '12').touch()
        self.make_images(self.before, '999', 1)
        self.assertEqual(sem._v29_discover_conditions(self.after), expected)

    def test_after_without_numeric_directories_has_clear_error(self):
        self.after.mkdir()
        with self.assertRaisesRegex(ValueError, '没有以纯数字命名'):
            sem._v29_discover_conditions(self.after)
        (self.after / 'notes' / '1').mkdir(parents=True)
        (self.after / '2').touch()
        with self.assertRaisesRegex(ValueError, '没有以纯数字命名'):
            sem.main_v29(self.before, self.after, self.root / 'output')

    def test_missing_or_non_directory_after_has_clear_error(self):
        with self.assertRaisesRegex(FileNotFoundError, 'after 根目录不存在'):
            sem._v29_discover_conditions(self.after)
        self.after.touch()
        with self.assertRaisesRegex(NotADirectoryError, 'after 路径不是文件夹'):
            sem._v29_discover_conditions(self.after)

    def test_leading_zeros_require_exact_before_folder_names(self):
        self.make_images(self.after, '007', 1)
        self.make_images(self.after, '7', 1)
        self.make_images(self.before, '7', 1)
        inv, errors = sem._v29_build_inventory(self.before, self.after)
        self.assertEqual(inv['condition'].unique().tolist(), ['7'])
        self.assertEqual(errors['condition'].tolist(), ['007'])
        self.assertEqual(errors['stage'].tolist(), ['before'])
        self.assertEqual(errors['error_type'].tolist(), ['MISSING_CONDITION_FOLDER'])
        self.make_images(self.before, '007', 1)
        inv, errors = sem._v29_build_inventory(self.before, self.after)
        self.assertTrue(errors.empty)
        self.assertEqual(inv['condition'].unique().tolist(), ['007', '7'])
        self.assertEqual(inv['pair_key'].unique().tolist(), ['C007_R2_trench', 'C7_R2_trench'])

    def test_all_missing_before_folders_still_export_diagnostics(self):
        self.make_images(self.after, '25', 1)
        self.before.mkdir()
        output = self.root / 'output'
        with contextlib.redirect_stdout(io.StringIO()), patch.object(sem, '_process_one_image_v19') as process:
            with self.assertRaisesRegex(RuntimeError, '没有建立任何有效 inventory'):
                sem.main_v29(self.before, self.after, output)
            process.assert_not_called()
        errors = sem.pd.read_csv(output / 'processing_errors.csv', dtype={'condition': str})
        self.assertEqual(errors['condition'].tolist(), ['25'])
        self.assertEqual(errors['error_type'].tolist(), ['MISSING_CONDITION_FOLDER'])
        self.assertEqual(errors['stage'].tolist(), ['before'])

    def test_cli_uses_after_folders_and_exports_numeric_order(self):
        valid = ['0', '2', '007', '10', '103']
        for condition in reversed(valid):
            self.make_images(self.before, condition, 1, names=['source.png'])
            self.make_images(self.after, condition, 1, names=['result.jpg'])
        self.make_images(self.before, '999', 1)
        self.make_images(self.after, '25', 1)  # No corresponding BEFORE folder.
        self.make_images(self.before, '71', 1)
        self.make_images(self.after, '71', 2)  # Count mismatch.
        (self.after / 'notes').mkdir()
        (self.after / '19').touch()
        output = self.root / 'SEM_0830_0831_compare_V29'
        with contextlib.redirect_stdout(io.StringIO()), \
                patch.object(sem, '_process_one_image_v19', side_effect=self.fake_measurement) as process, \
                patch.object(sem, '_v29_save_pair_plots') as plots:
            result = sem._v29_cli(['--before', str(self.before), '--after', str(self.after), '--via-pattern', 'via40'])
        self.assertEqual(result, 0)
        self.assertEqual([call.args[2] for call in process.call_args_list], [c for c in valid for _ in range(2)])
        for filename in ('inventory_mapping.csv', 'image_status.csv', 'all_object_measurements.csv',
                         'image_metric_summary.csv', 'image_pair_comparison.csv',
                         'overall_comparison.csv', 'object_before_after_comparison.csv'):
            with self.subTest(output=filename):
                data = sem.pd.read_csv(output / filename, dtype={'condition': str})
                self.assertEqual(data['condition'].drop_duplicates().tolist(), valid)
        errors = sem.pd.read_csv(output / 'processing_errors.csv', dtype={'condition': str})
        self.assertEqual(errors['condition'].tolist(), ['25', '71'])
        self.assertEqual(errors['error_type'].tolist(), ['MISSING_CONDITION_FOLDER', 'IMAGE_COUNT_MISMATCH'])
        settings = json.loads((output / 'settings.json').read_text(encoding='utf-8'))
        self.assertEqual(settings['conditions'], ['0', '2', '007', '10', '25', '71', '103'])
        self.assertIn('after_root', settings['condition_discovery_rule'])
        from openpyxl import load_workbook
        workbook = load_workbook(output / 'SEM_0830_0831_before_after_V29_results.xlsx', read_only=True)
        self.addCleanup(workbook.close)
        rows = list(workbook['inventory'].values)
        column = rows[0].index('condition')
        self.assertEqual(list(dict.fromkeys(row[column] for row in rows[1:])), valid)

        # Exercise real chart ordering while suppressing image file writes.
        from matplotlib.figure import Figure
        delta_labels = []
        def capture_labels(figure, filename, *args, **kwargs):
            if Path(filename).name.startswith('ALL__'):
                delta_labels.append([label.get_text() for label in figure.axes[0].get_xticklabels()])
        with patch.object(Figure, 'savefig', autospec=True, side_effect=capture_labels):
            sem._v29_save_pair_plots(plots.call_args.args[0], output / 'plot_order_check')
        self.assertEqual(len(delta_labels), 2)
        for labels in delta_labels:
            self.assertEqual(labels, [f'C{condition}-R2' for condition in valid])

    def test_equal_counts_have_no_fixed_limit_or_multiple_of_three_requirement(self):
        for count in (0, 1, 2, 3, 8, 9, 10, 11, 12, 37):
            with self.subTest(count=count):
                condition = str(100 + count)
                self.make_images(self.before, condition, count)
                self.make_images(self.after, condition, count, prefix='different')
                inv, errors = sem._v29_build_inventory(self.before, self.after, [condition])
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
        inv, errors = sem._v29_build_inventory(self.before, self.after, ['10'])
        before = inv[inv['stage'] == 'before'].sort_values('sequence_index')
        self.assertEqual(list(zip(before['region'], before['pattern'])), [
            (2, 'trench'), (2, 'slot'), (2, 'via'),
            (3, 'trench'), (3, 'slot'), (3, 'via'),
            (4, 'trench'), (4, 'slot'), (4, 'via'),
            (5, 'trench'), (5, 'slot'),
        ])

    def test_natural_order_mixed_image_formats_and_different_filenames(self):
        before_names = ['scan_1.TIF', 'scan_2.png', 'scan_3.JPG', 'scan_10.jpeg', 'scan_11.tiff']
        after_names = ['result_4.jpg', 'result_8.TIFF', 'result_9.JPEG', 'result_20.PNG', 'result_21.tif']
        before_folder = self.make_images(self.before, '10', 3, names=before_names)
        self.make_images(self.after, '10', 3, names=after_names)
        (before_folder / 'notes.csv').touch()
        (before_folder / 'notes.txt').touch()
        nested = before_folder / 'nested'
        nested.mkdir()
        (nested / 'extra.tif').touch()
        inv, errors = sem._v29_build_inventory(self.before, self.after, ['10'])
        self.assertTrue(errors.empty)
        self.assertEqual(len(inv), 10)
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
        inv, errors = sem._v29_build_inventory(self.before, self.after, conditions + ['14'])
        self.assertEqual(set(inv['condition']), {'14'})
        self.assertEqual(len(inv), 24)
        self.assertEqual(errors['error_type'].tolist(), ['IMAGE_COUNT_MISMATCH'] * 4)
        self.assertEqual(list(zip(errors['before_count'], errors['after_count'])), [(10, 9), (9, 10), (0, 1), (1, 0)])

    def test_missing_folder_does_not_leave_one_sided_inventory(self):
        self.make_images(self.before, '10', 10)
        self.make_images(self.after, '11', 0)
        inv, errors = sem._v29_build_inventory(self.before, self.after, ['10', '11'])
        self.assertTrue(inv.empty)
        self.assertEqual(errors['error_type'].tolist(), ['MISSING_CONDITION_FOLDER'] * 2)

    def test_metadata_failure_does_not_change_the_image_count(self):
        folder = self.make_images(self.before, '10', 10)
        self.make_images(self.after, '10', 10)
        (folder / 'image_10.txt').unlink()
        inv, errors = sem._v29_build_inventory(self.before, self.after, ['10'])
        self.assertEqual(len(inv), 20)
        self.assertEqual(errors['error_type'].tolist(), ['PIXELSIZE_PREFLIGHT_ERROR'])
        self.assertEqual(errors.iloc[0]['image'], 'image_10.tif')

    def test_all_invalid_inputs_still_export_count_diagnostics(self):
        self.make_images(self.before, '10', 10)
        self.make_images(self.after, '10', 9)
        output = self.root / 'output'
        with contextlib.redirect_stdout(io.StringIO()), patch.object(sem, '_process_one_image_v19') as process:
            with self.assertRaisesRegex(RuntimeError, '照片数量是否一致'):
                sem.main_v29(self.before, self.after, output)
            process.assert_not_called()
        errors = sem.pd.read_csv(output / 'processing_errors.csv')
        self.assertEqual(errors.iloc[0]['error_type'], 'IMAGE_COUNT_MISMATCH')
        self.assertEqual(errors.iloc[0]['before_count'], 10)
        self.assertEqual(errors.iloc[0]['after_count'], 9)

    def test_png_and_jpg_only_folders_pair_without_tifs(self):
        self.make_images(self.before, '10', 12, names=[f'image_{i}.PNG' for i in range(1, 13)])
        self.make_images(self.after, '10', 12, names=[f'result_{i}.jpg' for i in range(1, 13)])
        inv, errors = sem._v29_build_inventory(self.before, self.after, ['10'])
        self.assertTrue(errors.empty)
        self.assertEqual(len(inv), 24)
        self.assertEqual(inv['pair_key'].nunique(), 12)
        self.assertTrue((inv.groupby('pair_key')['stage'].nunique() == 2).all())

    def test_extra_png_counts_toward_a_mismatch(self):
        folder = self.make_images(self.before, '10', 9)
        self.make_images(self.after, '10', 9)
        (folder / 'image_10.PNG').touch()
        (folder / 'image_10.txt').write_text('PixelSize=0.8', encoding='utf-8')
        inv, errors = sem._v29_build_inventory(self.before, self.after, ['10'])
        self.assertTrue(inv.empty)
        self.assertEqual(errors['error_type'].tolist(), ['IMAGE_COUNT_MISMATCH'])
        self.assertEqual(errors.iloc[0]['before_count'], 10)
        self.assertEqual(errors.iloc[0]['after_count'], 9)

    def test_same_stem_different_formats_sort_deterministically(self):
        names = ['image_1.png', 'image_1.jpg', 'image_2.tif']
        self.make_images(self.before, '10', 3, names=names)
        self.make_images(self.after, '10', 3, names=list(reversed(names)))
        inv, errors = sem._v29_build_inventory(self.before, self.after, ['10'])
        self.assertTrue(errors.empty)
        for stage in ('before', 'after'):
            ordered = inv[inv['stage'] == stage].sort_values('sequence_index')
            self.assertEqual(ordered['image'].tolist(), ['image_1.jpg', 'image_1.png', 'image_2.tif'])

    def test_png_still_requires_matching_pixel_size_txt(self):
        folder = self.make_images(self.before, '10', 1, names=['photo.PNG'])
        self.make_images(self.after, '10', 1, names=['result.jpg'])
        (folder / 'photo.txt').unlink()
        inv, errors = sem._v29_build_inventory(self.before, self.after, ['10'])
        self.assertEqual(len(inv), 2)
        self.assertEqual(errors['error_type'].tolist(), ['PIXELSIZE_PREFLIGHT_ERROR'])
        self.assertEqual(errors.iloc[0]['image'], 'photo.PNG')

    def test_real_image_decoding_and_metadata_for_all_formats(self):
        folder = self.root / '中文 图片'
        folder.mkdir()
        gray = sem.np.arange(32 * 48, dtype=sem.np.uint16).reshape(32, 48)
        gray8 = (gray % 256).astype(sem.np.uint8)
        color = sem.np.stack([gray8, sem.np.flipud(gray8), gray8], axis=2)
        alpha = sem.np.dstack([color, sem.np.full_like(gray8, 255)])
        cases = [('gray.tif', gray8), ('depth.TIFF', gray), ('gray.png', gray8),
                 ('color.PNG', color), ('alpha.png', alpha), ('color.JPG', color),
                 ('gray.JPEG', gray8)]
        for name, data in cases:
            with self.subTest(image=name):
                image = folder / name
                sem._imwrite_unicode_v19(image, data)
                image.with_suffix('.TXT').write_text('PixelSize=0.8', encoding='utf-8')
                decoded = sem._read_gray_image_unicode_v19(image)
                self.assertEqual(decoded.shape, (32, 48))
                self.assertEqual(decoded.dtype, sem.np.uint8)
                self.assertGreater(int(decoded.max()), int(decoded.min()))
                self.assertAlmostEqual(sem.read_pixel_size_nm(image).value_nm, 0.8)

    def test_pipeline_exports_all_pairs_beyond_nine(self):
        for condition, bcount, acount in [('10', 13, 13), ('11', 8, 9), ('12', 2, 2)]:
            suffixes = ('.tif', '.PNG', '.jpg', '.TIFF', '.jpeg')
            before_names = [f'image_{i}{suffixes[(i - 1) % 5]}' for i in range(1, bcount + 1)]
            after_names = [f'result_{i}{suffixes[i % 5]}' for i in range(1, acount + 1)]
            self.make_images(self.before, condition, bcount, names=before_names)
            self.make_images(self.after, condition, acount, names=after_names)
        output = self.root / 'output'

        with contextlib.redirect_stdout(io.StringIO()), \
                patch.object(sem, '_process_one_image_v19', side_effect=self.fake_measurement) as process, \
                patch.object(sem, '_v29_save_pair_plots'):
            outputs = sem.main_v29(self.before, self.after, output, via_pattern_mode='via40')
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
        self.assertEqual(set(last_pair['before_image']), {'image_13.jpg'})
        self.assertEqual(set(last_pair['after_image']), {'result_13.TIFF'})
        self.assertEqual(errors['error_type'].tolist(), ['IMAGE_COUNT_MISMATCH'])
        settings = json.loads(outputs['settings'].read_text(encoding='utf-8'))
        self.assertIsNone(settings['required_image_count_per_condition_stage'])
        self.assertEqual(settings['regions'], [2, 3, 4, 5, 6])
        self.assertEqual(settings['supported_image_extensions'], ['.jpeg', '.jpg', '.png', '.tif', '.tiff'])
        self.assertIn('V29', settings['script_version'])
        self.assertEqual(outputs['excel'].name, 'SEM_0830_0831_before_after_V29_results.xlsx')
        from openpyxl import load_workbook
        workbook = load_workbook(outputs['excel'], read_only=True)
        self.addCleanup(workbook.close)
        self.assertEqual(workbook['inventory'].max_row, 31)
        self.assertEqual(workbook['image_pair_compare'].max_row, 31)


if __name__ == '__main__':
    unittest.main()
