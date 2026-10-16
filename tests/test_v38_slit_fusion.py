"""V35/V37/V38 slit 回退对照；背景/邻位/模糊不能改变真实几何。"""
import contextlib
import importlib
import io
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
VERSIONS = {}
for version in (35, 37, 38):
    sys.path.insert(0, str(ROOT / f'SEM_0830_0831_compare_V{version}_bundle'))
    VERSIONS[version] = importlib.import_module(f'sem_before_after_compare_0831_V{version}')
sem = VERSIONS[38]
np, cv = sem.np, sem.cv2
from test_v38_edges import partial_feature


def slit_scene(width=16, pitch=120, texture=0., blur=.8, noise=.2, missing=()):
    yy, xx = np.indices((400, 2*pitch+2*max(60,width)))
    background = 150.+.05*xx+.03*yy+texture*np.sin(xx/29.)*np.cos(yy/37.)
    mask = np.zeros(xx.shape, np.uint8)
    for r, y in enumerate((55,145,235,325)):
        for c in range(3):
            if (r,c) in missing:
                continue
            x, hh = max(60,width)+c*pitch, (55,43,31,21)[r]//2
            cv.rectangle(mask, (x-width//2,y-hh), (x+width//2,y+hh), 1, -1)
    raw = background-80.*cv.GaussianBlur(mask.astype(np.float32),(0,0),blur)
    raw += np.random.default_rng(38).normal(0,noise,raw.shape)
    return np.clip(raw,0,255).astype(np.uint8)


def run_slots(raw, version=38, scale=1., **options):
    module = VERSIONS[version]
    pool, _ = module.extract_descriptors_with_rejections(module.segment_features(raw))
    with contextlib.redirect_stdout(io.StringIO()), module._v33_detection_context(options):
        return module._select_and_measure_bottom_slots_v19(
            raw, pool, scale, 'slit.png', 'single', '007', 'slot210')


class SlitFusionTests(unittest.TestCase):
    def test_curved_illumination_reproduces_v37_loss_and_preserves_v35_sites(self):
        raw = slit_scene(texture=4.)
        old, _, _, _, _ = run_slots(raw,35)
        lost, _, _, rejected, _ = run_slots(raw,37)
        rows, selected, bottom, _, meta = run_slots(raw)
        self.assertEqual((len(old),len(lost),len(rows),len(selected),len(bottom)), (3,0,3,12,3))
        self.assertTrue(any(r.get('edge_status') == 'UNRESOLVED_OUTER_BOUNDARY' for r in rejected))
        np.testing.assert_allclose([[r['x_width_nm'],r['y_height_nm']] for r in rows],
                                   [[17.,21.]]*3, atol=.5)
        self.assertEqual(meta['v38_slit_primary_accepted_count'],12)
        self.assertTrue(all(r['slit_locator_path'] == 'V35_primary' for r in rows))

    def test_stronger_curved_background_and_bright_polarity_keep_real_dimensions(self):
        raw = slit_scene(texture=8.)
        for polarity in ('dark','bright'):
            with self.subTest(polarity=polarity), patch.object(sem,'FEATURE_POLARITY',polarity):
                rows, _, _, _, _ = run_slots(raw if polarity=='dark' else 255-raw, scale=2.,
                                             slit_width_nm=34.,slit_length_nm=42.)
            self.assertEqual(len(rows),3)
            np.testing.assert_allclose([[r['x_width_nm'],r['y_height_nm']] for r in rows],
                                       [[34.,42.]]*3, atol=1.)

    def test_dense_neighbors_do_not_count_as_unexplained_outer_foreground(self):
        raw = slit_scene(width=64,pitch=80)
        old, _, _, _, _ = run_slots(raw,35)
        lost, _, _, _, _ = run_slots(raw,37)
        rows, selected, _, _, _ = run_slots(raw)
        self.assertEqual((len(old),len(lost),len(rows),len(selected)), (3,0,3,12))
        np.testing.assert_allclose([[r['x_width_nm'],r['y_height_nm']] for r in rows],
                                   [[65.,21.]]*3, atol=.6)

    def test_blur_tail_does_not_expand_the_half_gray_measurement(self):
        raw = slit_scene(blur=3.,noise=1.)
        old, _, _, _, _ = run_slots(raw,35)
        inflated, _, _, _, _ = run_slots(raw,37)
        rows, _, _, _, _ = run_slots(raw)
        self.assertEqual((len(old),len(inflated),len(rows)),(3,3,3))
        self.assertTrue(all(r['x_width_nm'] > 20. for r in inflated))
        actual = np.array([[r['x_width_nm'],r['y_height_nm']] for r in rows])
        np.testing.assert_allclose(actual, [[r['x_width_nm'],r['y_height_nm']] for r in old], atol=.1)
        np.testing.assert_allclose(actual, [[17.,21.]]*3, atol=1.5)
        self.assertTrue(all(not r['boundary_relocalized'] for r in rows))

    def test_slit_specific_path_still_recovers_a_weak_lower_half(self):
        for weak in (6.,12.):
            for polarity in ('dark','bright'):
                raw, seed = partial_feature('slot',weak=weak)
                with self.subTest(weak=weak,polarity=polarity), patch.object(sem,'FEATURE_POLARITY',polarity):
                    final, audit = sem._v32_refine_desc(raw if polarity=='dark' else 255-raw,
                                                       sem.describe_contour(seed),'slot')
                self.assertIsNotNone(final,audit)
                np.testing.assert_allclose(np.ptp(final.contour[:,0],axis=0),[17.,71.],atol=1.)
                self.assertTrue(audit['boundary_relocalized'])
                self.assertFalse(audit['edge_template_forced'])

    def test_actual_short_slit_is_not_extended(self):
        raw, seed = partial_feature('slot',actual_lower=False)
        final, audit = sem._v32_refine_desc(raw,sem.describe_contour(seed),'slot')
        self.assertIsNotNone(final,audit)
        self.assertAlmostEqual(np.ptp(final.contour[:,0,1]),35.,delta=1.)
        self.assertFalse(audit['boundary_relocalized'])

    def test_failed_narrow_seed_is_repaired_only_from_actual_image_pixels(self):
        seed = np.array([[[179,315]],[[181,315]],[[181,325]],[[181,335]],
                         [[179,335]],[[179,325]]],np.float32)
        desc = sem.describe_contour(seed)
        raw = slit_scene()
        for blank in (False,True):
            pairs, rejects = [], []
            measurement = np.full(raw.shape,180,np.float32) if blank else raw
            with self.subTest(blank=blank), patch.object(sem,'_V32_MEASUREMENT_GRAY',measurement):
                added, repaired = sem._v38_refine_slot_observations(
                    raw,[desc],[{'slot_row':4}],{'slot_horizontal_pitch_px':120.,'slot_row_pitch_px':90.},
                    pairs,rejects,'V35_primary')
            self.assertEqual((added,repaired),(0,0) if blank else (1,1))
            if not blank:
                np.testing.assert_allclose(np.ptp(pairs[0][0].contour[:,0],axis=0),[17.,21.],atol=.6)
            else:
                self.assertTrue(rejects)

    def test_one_complete_column_survives_when_other_columns_lack_array_support(self):
        raw = slit_scene(missing=((1,1),(2,1),(0,2),(1,2),(2,2)))
        old, _, _, _, _ = run_slots(raw,35)
        lost, _, _, _, _ = run_slots(raw,37)
        rows, selected, _, _, meta = run_slots(raw)
        self.assertEqual((len(old),len(lost),len(rows),len(selected)),(1,0,1,4))
        self.assertTrue(meta['v38_slit_complete_column_support_only'])
        self.assertEqual(rows[0]['slot_observed_rows'],'1,2,3,4')
        self.assertFalse(rows[0]['slot_is_partial_column'])
        np.testing.assert_allclose([rows[0]['x_width_nm'],rows[0]['y_height_nm']],[17.,21.],atol=.5)

    def test_whole_via_pipeline_preserves_v37_weak_and_elongated_objects(self):
        from test_v38_edges import shaded_vias
        from test_v38_morphology import via_scene
        for raw in (shaded_vias(),via_scene(35.)):
            results=[]
            for version in (37,38):
                module=VERSIONS[version]
                pool,_=module.extract_descriptors_with_rejections(module.segment_features(raw))
                with contextlib.redirect_stdout(io.StringIO()):
                    rows, selected, _, _=module._select_and_measure_vias_v19(
                        raw,pool,1.,'via.png','single','007','via60')
                # Match actual sites before comparing measurements. Tiny center
                # changes must not reorder separate rows by raw floating X.
                values=[(r['center_image_x_px'],r['center_image_y_px'],
                         r['x_width_nm'],r['y_height_nm']) for r in rows]
                results.append(sorted(values,key=lambda t:(round(t[1]/5),round(t[0]/5))))
                self.assertEqual(len(rows),20 if raw.shape==(360,430) else 18)
            np.testing.assert_allclose(results[0],results[1],atol=1.)
            pitch,margin=(75,65) if raw.shape==(360,430) else (85,60)
            expected=np.array([(margin+c*pitch,margin+r*pitch) for r in range(4) for c in range(5)
                               if raw.shape==(360,430) or (r,c) not in ((1,1),(2,2))])
            np.testing.assert_allclose(np.asarray(results[1])[:,:2],expected,atol=1.)


if __name__=='__main__':
    unittest.main()
