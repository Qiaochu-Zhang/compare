"""V37 接纳平滑拉长 via、拒绝异常边缘，允许 slit 阵列局部缺失。"""
import contextlib
import io
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'SEM_0830_0831_compare_V37_bundle'))
import sem_before_after_compare_0831_V37 as sem
from test_v30 import gradient_slots
from test_v32 import analytic_trenches

np, cv = sem.np, sem.cv2


def contour(kind='ellipse', angle=0.):
    t = np.linspace(0, 2*np.pi, 160, endpoint=False)
    radius = (1 + .35*np.cos(8*t) if kind == 'star' else
              1 + .25*np.cos(3*t) if kind == 'lobed' else
              1 + .20*np.cos(24*t) if kind == 'saw' else np.ones_like(t))
    points = np.column_stack([30*radius*np.cos(t), (10 if kind == 'ellipse' else 30)*radius*np.sin(t)])
    a = math.radians(angle)
    rotation = np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])
    return (points @ rotation.T + [60,60]).astype(np.float32).reshape(-1,1,2)


def via_scene(angle=0., defects=True, mixed=False):
    yy, xx = np.mgrid[:375, :460]
    mask = np.zeros(xx.shape, np.uint8)
    for r, y in enumerate((60, 145, 230, 315)):
        for c, x in enumerate((60, 145, 230, 315, 400)):
            if defects and (r,c) == (1,1):
                continue
            if defects and (r,c) == (2,2):
                points = contour('star').reshape(-1,2) - [60,60] + [x,y]
                cv.fillPoly(mask, [points.astype(np.int32)], 1)
            else:
                axes = (30,10)
                if mixed:
                    axes = (36,12) if (r,c) == (2,2) else (15,15)
                cv.ellipse(mask, (x,y), axes, angle, 0, 360, 1, -1)
    return (180+.02*xx+.01*yy-110*cv.GaussianBlur(mask.astype(np.float32),(0,0),.8)).astype(np.uint8)


def polluted_slots(missing=(), stain=False):
    image = gradient_slots(40.)
    yy,xx=np.indices(image.shape)
    background = 80.+.24*xx+.12*yy
    for r,c in missing:
        y=(55,145,235,325)[r]; x=(60,180,300)[c]
        image[y-35:y+36,x-20:x+21]=background[y-35:y+36,x-20:x+21].astype(np.uint8)
    if stain:
        points = contour('star').reshape(-1,2) * .65 + [68,130]
        cv.fillPoly(image,[points.astype(np.int32)],25)
    return image


def run_vias(gray, **options):
    primary,_=sem.extract_descriptors_with_rejections(sem.segment_features(gray))
    with sem._v33_detection_context(options), contextlib.redirect_stdout(io.StringIO()):
        return sem._select_and_measure_vias_v19(gray,primary,1.,'via.png','single','007','via60')


def run_slots(gray, **options):
    with sem._v33_detection_context(options), contextlib.redirect_stdout(io.StringIO()):
        return sem._select_and_measure_bottom_slots_v19(gray,[],1.,'slot.png','single','007','slot210')


class MorphologyTests(unittest.TestCase):
    def test_defaults_disable_circularity_and_peer_size_clipping_but_allow_explicit_controls(self):
        options=sem.V37DetectionOptions()
        self.assertEqual((options.via_min_circularity,options.via_max_size_deviation,options.slit_min_column_rows),(0.,0.,2))
        args=sem._v32_build_parser().parse_args(['--via-min-circularity','.8','--via-max-size-deviation','.3',
                                               '--slot-min-column-rows','4'])
        parsed=sem.V37DetectionOptions(**{k:getattr(args,k) for k in sem.asdict(options)})
        self.assertEqual((parsed.via_min_circularity,parsed.via_max_size_deviation,parsed.slit_min_column_rows),(.8,.3,4))
        for invalid in (1,5,2.5,True,'2'):
            with self.subTest(value=invalid), self.assertRaises(ValueError):
                sem.V37DetectionOptions(slit_min_column_rows=invalid)

    def test_ellipse_normalization_is_rotation_invariant_without_altering_contour(self):
        for angle in (0.,27.,73.):
            points=contour('ellipse',angle); original=points.copy()
            ok,metrics,reasons=sem._via_shape_qc_v24(points,40.,'single','via60')
            self.assertTrue(ok,reasons)
            self.assertLess(metrics['circularity'],.70)
            self.assertGreater(metrics['aspect'],2.9)
            self.assertLess(metrics['via_ellipse_residual_p95'],.01)
            np.testing.assert_array_equal(points,original)

    def test_jagged_and_lobed_shapes_fail_even_with_roundness_disabled_and_fallback(self):
        for kind in ('star','saw','lobed'):
            for level in (0,1,2):
                with self.subTest(kind=kind,level=level), patch.object(sem,'_V25_ZERO_FALLBACK_LEVEL',level,create=True):
                    ok,metrics,reasons=sem._via_shape_qc_v24(contour(kind),60.,'after','via40')
                self.assertFalse(ok)
                self.assertTrue(reasons)
                if kind=='lobed':
                    self.assertGreater(metrics['circularity'],.75)  # roundness alone would accept it
                    self.assertGreater(metrics['via_ellipse_residual_p95'],.18)

    def test_real_oval_arrays_survive_missing_site_and_star_without_duplicate_tips(self):
        for angle in (0.,35.):
            with self.subTest(angle=angle):
                rows,selected,rejects,bundle=run_vias(via_scene(angle))
                self.assertEqual(len(rows),18)
                centers=np.array([d.center for d in selected])
                self.assertEqual(len({(round(x/85),round(y/85)) for x,y in centers}),18)
                self.assertTrue(all(r['edge_signal_source']=='original_image_gray' for r in rows))
                self.assertTrue(all(not r['edge_template_forced'] for r in rows))
                self.assertLess(max(r['circularity'] for r in rows),.75)
                self.assertTrue(any(r['reason_code'] in ('via_irregular_boundary','original_image_edge_not_supported') for r in rejects))
                self.assertTrue(all(np.linalg.norm(d.center-[230,230])>30 for d in selected))
                self.assertTrue(all(np.linalg.norm(d.center-[145,145])>30 for d in selected))
                a=math.radians(angle)
                expected_x=2*math.sqrt(30**2*math.cos(a)**2+10**2*math.sin(a)**2)+1
                expected_y=2*math.sqrt(30**2*math.sin(a)**2+10**2*math.cos(a)**2)+1
                np.testing.assert_allclose([r['x_width_nm'] for r in rows],expected_x,atol=2.)
                np.testing.assert_allclose([r['y_height_nm'] for r in rows],expected_y,atol=2.)
                self.assertEqual(bundle['meta']['via_v24_final_count'],18)

    def test_smooth_larger_oval_is_not_removed_by_default_peer_size_filter(self):
        gray=via_scene(defects=False,mixed=True)
        rows,selected,_,_=run_vias(gray)
        self.assertEqual(len(rows),20)
        oval=next(r for r in rows if abs(r['center_image_x_px']-230)<2 and abs(r['center_image_y_px']-230)<2)
        self.assertGreater(oval['via_size_deviation_fraction'],.30)
        self.assertFalse(oval['via_size_filter_applied'])
        filtered,_,rejects,_=run_vias(gray,via_max_size_deviation=.30)
        self.assertEqual(len(filtered),19)
        self.assertTrue(any(r['reason_code']=='via_size_outlier' for r in rejects))

    def test_explicit_roundness_can_still_request_strict_filtering(self):
        rows,_,rejects,_=run_vias(via_scene(),via_min_circularity=.9)
        self.assertEqual(rows,[])
        self.assertTrue(any(r['reason_code']=='via_circularity_below_threshold' for r in rejects))

    def test_ellipse_locator_on_blank_original_cannot_create_measurements(self):
        gray=via_scene()
        with patch.object(sem,'_V32_MEASUREMENT_GRAY',np.full(gray.shape,180,np.float32)):
            rows,_,rejects,_=run_vias(gray)
        self.assertEqual(rows,[])
        self.assertTrue(any(r['reason_code']=='original_image_edge_not_supported' for r in rejects))

    def test_upper_holes_in_every_column_keep_three_measurable_bottom_slots(self):
        rows,selected,bottom,rejects,meta=run_slots(polluted_slots([(0,0),(1,1),(2,2)]))
        self.assertEqual((len(rows),len(selected),len(bottom)),(3,9,3))
        self.assertEqual(meta['v23_complete_slot_columns'],0)
        self.assertEqual(meta['v36_slot_partial_measured_count'],3)
        self.assertEqual(meta['v36_slot_missing_site_count'],3)
        self.assertEqual([r['slot_missing_rows'] for r in rows],['1','2','3'])
        self.assertEqual(sum(r['reason_code']=='slot_missing_or_obscured_site' for r in rejects),3)
        np.testing.assert_allclose([r['x_width_nm'] for r in rows],17.,atol=1.)
        np.testing.assert_allclose([r['y_height_nm'] for r in rows],21.,atol=1.)
        self.assertTrue(all(r['edge_support_fraction']>=.85 for r in rows))

    def test_two_missing_upper_members_allowed_and_strict_column_mode_available(self):
        gray=polluted_slots([(0,0),(1,0)])
        rows,_,_,_,_=run_slots(gray)
        self.assertEqual(len(rows),3)
        self.assertEqual(rows[0]['slot_observed_rows'],'3,4')
        strict,_,_,rejects,_=run_slots(gray,slit_min_column_rows=4)
        self.assertEqual(len(strict),2)
        self.assertTrue(any(r['reason_code']=='slot_insufficient_column_support' for r in rejects))

    def test_missing_bottom_skipped_and_entire_missing_bottom_never_promotes_upper_row(self):
        rows,selected,bottom,_,meta=run_slots(polluted_slots([(3,1)]))
        self.assertEqual((len(rows),len(selected),len(bottom)),(2,11,2))
        self.assertEqual(meta['v36_slot_missing_bottom_count'],1)
        np.testing.assert_allclose([r['center_image_y_px'] for r in rows],325.,atol=1.)
        rows,_,_,_,meta=run_slots(polluted_slots([(3,0),(3,1),(3,2)]))
        self.assertEqual(rows,[])
        self.assertFalse(meta['v36_slot_grid_supported'])

    def test_unrelated_stain_does_not_change_bottom_sizes_or_create_extra_columns(self):
        rows,_,_,_,_=run_slots(polluted_slots([(0,0)],stain=True))
        self.assertEqual(len(rows),3)
        np.testing.assert_allclose([r['x_width_nm'] for r in rows],17.,atol=1.)
        np.testing.assert_allclose([r['y_height_nm'] for r in rows],21.,atol=1.)

    def test_no_edges_or_only_isolated_bottom_row_do_not_establish_an_array(self):
        for gray in (np.full((400,360),140,np.uint8),
                     polluted_slots([(r,c) for r in range(3) for c in range(3)])):
            rows,_,_,_,_=run_slots(gray)
            self.assertEqual(rows,[])

    def test_single_cli_exports_morphology_metrics_missing_sites_and_same_stats(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); folder=root/'输入'/'007'; folder.mkdir(parents=True)
            for i,gray in enumerate((analytic_trenches(),polluted_slots([(0,0),(1,1),(2,2)]),via_scene()),1):
                sem._imwrite_unicode_v19(folder/f'{i}.png',gray)
                (folder/f'{i}.txt').write_text('PixelSize=2')
            output=root/'结果'
            result=subprocess.run([sys.executable,str(Path(sem.__file__)),'--input',str(folder.parent),
                                   '--output',str(output),'--image-timeout-seconds','90'],
                                  capture_output=True,text=True,timeout=240,
                                  env={**os.environ,'MPLBACKEND':'Agg','OPENBLAS_NUM_THREADS':'1','PYTHONDONTWRITEBYTECODE':'1'})
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            data=sem.pd.read_csv(output/'all_object_measurements.csv',dtype={'condition':str})
            self.assertEqual(len(data),30)
            self.assertEqual(set(data.condition),{'007'})
            self.assertEqual(set(data.stage),{'single'})
            slots=data[data.pattern=='slot']; vias=data[data.pattern=='via']
            self.assertEqual(len(vias),18)
            self.assertTrue(slots.slot_is_partial_column.astype(bool).all())
            self.assertTrue((vias.circularity<.75).all())
            self.assertTrue((vias.via_ellipse_residual_p95<.18).all())
            np.testing.assert_allclose(slots.x_width_nm,34.,atol=2.)
            np.testing.assert_allclose(vias.x_width_nm,122.,atol=3.)
            metrics=sem.pd.read_csv(output/'image_metric_summary.csv')
            self.assertEqual(metrics[metrics.pattern=='via'].n.tolist(),[18,18])
            status=sem.pd.read_csv(output/'image_status.csv')
            self.assertTrue(all(Path(p).is_file() for p in status.annotated_path))
            self.assertEqual(len(list((output/'plots').glob('*.png'))),18)
            rejected=sem.pd.read_csv(output/'rejected_objects.csv')
            self.assertEqual(sum(rejected.reason_code=='slot_missing_or_obscured_site'),3)
            settings=json.loads((output/'settings.json').read_text())
            self.assertEqual(settings['detection_options']['slit_min_column_rows'],2)
            self.assertEqual(settings['detection_options']['via_min_circularity'],0.)
            self.assertIn('V37',settings['script_version'])
            self.assertTrue((output/'SEM_0830_0831_single_V37_results.xlsx').is_file())


if __name__=='__main__':
    unittest.main()
