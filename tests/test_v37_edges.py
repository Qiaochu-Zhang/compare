"""V37 reproduces false internal caps and verifies raw outer-boundary recovery."""
import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'SEM_0830_0831_compare_V37_bundle'))
import sem_before_after_compare_0831_V37 as sem
from test_v30 import gradient_slots
np,cv=sem.np,sem.cv2


def partial_feature(pattern='via', weak=12., actual_lower=True, jagged=False):
    yy,xx=np.indices((180,180)); mask=np.zeros(xx.shape,np.uint8)
    if pattern=='slot':
        cv.rectangle(mask,(72,45),(88,115),1,-1)
    else:
        cv.ellipse(mask,(80,85),(22,38),0,0,360,1,-1)
    if jagged:
        t=np.linspace(0,2*np.pi,80,endpoint=False)
        radii=1.+.35*np.cos(8*t)
        points=np.column_stack([80+24*radii*np.cos(t),85+38*radii*np.sin(t)]).astype(np.int32)
        mask[:]=0;cv.fillPoly(mask,[points],1)
    if not actual_lower:mask[80:]=0
    amplitude=np.where(yy<80,90.,weak)
    image=(190.+.03*xx+.02*yy-cv.GaussianBlur((mask*amplitude).astype(np.float32),(0,0),.8)).astype(np.float32)
    initial=(mask*(yy<80)).astype(np.uint8)
    contours,_=cv.findContours(initial,cv.RETR_EXTERNAL,cv.CHAIN_APPROX_NONE)
    return image,contours[0].astype(np.float32)


def shaded_slots():
    yy,xx=np.indices((400,360));mask=np.zeros(xx.shape,np.float32)
    for row,y in enumerate((55,145,235,325)):
        height=(55,43,31,21)[row]
        for col,x in enumerate((60,180,300)):
            amp=np.where(yy<y,80.,10.)
            if (row,col)==(3,1):amp=np.full(xx.shape,7.)
            member=np.zeros(xx.shape,np.uint8)
            cv.rectangle(member,(x-8,y-height//2),(x+8,y+height//2),1,-1)
            mask+=member*amp
    gray=180.+.03*xx+.02*yy-cv.GaussianBlur(mask,(0,0),.8)
    gray+=np.random.default_rng(37).normal(0,.3,gray.shape)
    return np.clip(gray,0,255).astype(np.uint8)


def shaded_vias():
    yy,xx=np.indices((360,430));signal=np.zeros(xx.shape,np.float32)
    for row,y in enumerate((65,140,215,290)):
        for col,x in enumerate((65,140,215,290,365)):
            member=np.zeros(xx.shape,np.uint8)
            cv.circle(member,(x,y),20,1,-1)
            amp=np.where(yy<y,90.,12.)
            if (row,col)==(1,3):amp=np.full(xx.shape,7.)
            signal+=member*amp
    image=180.+.02*xx+.02*yy-cv.GaussianBlur(signal,(0,0),.8)
    return np.clip(image,0,255).astype(np.uint8)


def shaded_trenches(mode='steps'):
    yy,xx=np.mgrid[:400,:720]
    sigmoid=lambda z:1/(1+np.exp(np.clip(-z,-100,100)))
    signal=np.zeros(xx.shape,float)
    for i,x in enumerate(range(40,720,80)):
        amplitude=8. if mode=='weak_bars' and i%2 else 110.
        shape=sigmoid((xx-(x-10))/.8)*sigmoid(((x+10)-xx)/.8)
        top=sigmoid((170-yy)/.8);bottom=sigmoid((yy-230)/.8)
        ta=np.where(yy<140,amplitude,.10*amplitude) if mode=='steps' else amplitude
        ba=np.where(yy>265,amplitude,.10*amplitude) if mode=='steps' else (8. if mode=='weak_lower' else amplitude)
        signal+=shape*(top*ta+bottom*ba)
    return (170.+.03*xx+.02*yy-signal).astype(np.uint8)


class FullBoundaryTests(unittest.TestCase):
    def test_reproduces_internal_caps_then_recovers_true_full_slit_and_via(self):
        for pattern,expected in (('slot',(17.,71.)),('via',(45.,77.))):
            with self.subTest(pattern=pattern):
                gray,initial=partial_feature(pattern)
                old,old_audit=sem._v32_supported_contour_legacy(gray,initial)
                self.assertIsNotNone(old)
                self.assertEqual(old_audit['edge_support_fraction'],1.)
                self.assertLess(np.ptp(old[:,0,1]),40.)
                final,audit=sem._v32_supported_contour(gray,initial)
                self.assertIsNotNone(final,audit)
                np.testing.assert_allclose(np.ptp(final[:,0],axis=0),expected,atol=1.)
                self.assertTrue(audit['boundary_relocalized'])
                self.assertGreater(audit['boundary_seed_outer_contact_fraction'],.12)
                self.assertLessEqual(audit['boundary_outer_contact_fraction'],.12)
                self.assertEqual(audit['edge_signal_source'],'original_image_gray')
                self.assertFalse(audit['edge_template_forced'])

    def test_weak_bottom_and_bright_polarity_are_measured_on_raw_edges(self):
        gray,initial=partial_feature('via',weak=6.)
        for polarity in ('dark','bright'):
            with self.subTest(polarity=polarity),patch.object(sem,'FEATURE_POLARITY',polarity):
                image=gray if polarity=='dark' else 255.-gray
                final,audit=sem._v32_supported_contour(image,initial)
            self.assertIsNotNone(final,audit)
            self.assertAlmostEqual(np.ptp(final[:,0,1]),77.,delta=1.)

    def test_actual_short_object_is_not_extended_to_a_template(self):
        gray,initial=partial_feature('slot',actual_lower=False)
        final,audit=sem._v32_supported_contour(gray,initial)
        self.assertIsNotNone(final,audit)
        self.assertAlmostEqual(np.ptp(final[:,0,1]),35.,delta=1.)
        self.assertFalse(audit['boundary_relocalized'])

    def test_flat_or_sloped_original_does_not_accept_an_enhanced_locator(self):
        _,initial=partial_feature('via')
        yy,xx=np.indices((180,180))
        for raw in (np.full(xx.shape,180,np.float32),(140.+.1*xx+.1*yy).astype(np.float32)):
            final,audit=sem._v32_supported_contour(raw,initial)
            self.assertIsNone(final,audit)
        locator,_=partial_feature('via')
        desc=sem.describe_contour(initial)
        with patch.object(sem,'_V32_MEASUREMENT_GRAY',np.full(locator.shape,180,np.float32)):
            refined,_=sem._v32_refine_desc(locator,desc,'via')
        self.assertIsNone(refined)

    def test_irregular_real_outer_boundary_is_not_replaced_by_an_ellipse(self):
        gray,initial=partial_feature('via',jagged=True)
        final,audit=sem._v32_supported_contour(gray,initial)
        if final is not None:
            ok,_,_=sem._via_shape_qc_v24(final,60.,'single','via60')
            self.assertFalse(ok)
            self.assertFalse(audit['edge_template_forced'])
        else:
            self.assertIn(audit['edge_status'],('INSUFFICIENT_ORIGINAL_IMAGE_SUPPORT','SHAPE_CHANGE_EXCEEDS_LIMIT','UNRESOLVED_OUTER_BOUNDARY'))

    def test_shaded_slots_measure_full_bottom_objects_despite_partial_success(self):
        with contextlib.redirect_stdout(io.StringIO()):
            rows,selected,bottom,rejected,meta=sem._select_and_measure_bottom_slots_v19(
                shaded_slots(),[],1.,'slot.png','single','007','slot210')
        self.assertEqual((len(rows),len(bottom)),(3,3))
        np.testing.assert_allclose([r['x_width_nm'] for r in rows],17.,atol=1.)
        np.testing.assert_allclose([r['y_height_nm'] for r in rows],21.,atol=1.)
        np.testing.assert_allclose([r['center_image_y_px'] for r in rows],325.,atol=1.)
        self.assertTrue(all(not r['edge_template_forced'] for r in rows))

    def test_shaded_via_array_includes_whole_weak_hole_and_full_lower_halves(self):
        gray=shaded_vias()
        pool,_=sem.extract_descriptors_with_rejections(sem.segment_features(gray))
        with contextlib.redirect_stdout(io.StringIO()):
            rows,selected,rejected,bundle=sem._select_and_measure_vias_v19(gray,pool,1.,'via.png','single','007','via60')
        self.assertEqual(len(rows),20)
        self.assertGreater(bundle['meta']['v37_local_contrast_candidates'],0)
        np.testing.assert_allclose([r['x_width_nm'] for r in rows],41.,atol=2.)
        np.testing.assert_allclose([r['y_height_nm'] for r in rows],41.,atol=2.)
        self.assertEqual(len({(round(d.center[0]/75),round(d.center[1]/75)) for d in selected}),20)
        self.assertTrue(all(r['edge_signal_source']=='original_image_gray' for r in rows))

    def test_trenches_use_real_gap_across_internal_steps_and_weak_bodies(self):
        for mode in ('steps','weak_lower','weak_bars'):
            with self.subTest(mode=mode),contextlib.redirect_stdout(io.StringIO()):
                rows,pairs,rejected,meta=sem._select_and_measure_trenches_v19(
                    shaded_trenches(mode),[],1.,'trench.png','single','007','trench160')
            self.assertEqual(len(rows),9)
            np.testing.assert_allclose([r['x_width_nm'] for r in rows],20.,atol=1.)
            np.testing.assert_allclose([r['tip_gap_y_nm'] for r in rows],60.,atol=1.)
            np.testing.assert_allclose([p['top_tip_y_px'] for p in pairs],170.,atol=1.)
            np.testing.assert_allclose([p['bottom_tip_y_px'] for p in pairs],230.,atol=1.)
            self.assertEqual(meta['v37_trench_raw_pairs'],9)
            self.assertTrue(all(r['edge_method']=='V37_observed_outer_plateau50' for r in rows))

    def test_trench_polarity_hints_and_scale_preserve_real_dimensions(self):
        with sem._v33_detection_context(dict(trench_width_nm=40.,trench_gap_nm=120.)),patch.object(sem,'FEATURE_POLARITY','bright'):
            rows,pairs,_,_=sem._select_and_measure_trenches_v19(
                255-shaded_trenches(),[],2.,'trench.png','after','007','trench160')
        self.assertEqual(len(rows),9)
        np.testing.assert_allclose([r['x_width_nm'] for r in rows],40.,atol=2.)
        np.testing.assert_allclose([r['tip_gap_y_nm'] for r in rows],120.,atol=2.)

    def test_noise_and_single_sided_bars_do_not_create_trench_pairs(self):
        random=np.random.default_rng(37)
        yy,xx=np.indices((400,720))
        for image in ((160.+.03*xx+.02*yy+random.normal(0,.4,xx.shape)).astype(np.uint8),
                      np.full(xx.shape,160,np.uint8)):
            rows,_,_,_=sem._select_and_measure_trenches_v19(image,[],1.,'blank.png','single','007','trench160')
            self.assertEqual(rows,[])
        image=shaded_trenches();image[200:]=np.clip(170.+.03*xx[200:]+.02*yy[200:],0,255).astype(np.uint8)
        rows,_,_,_=sem._select_and_measure_trenches_v19(image,[],1.,'one_side.png','single','007','trench160')
        self.assertEqual(rows,[])

    def test_real_cli_exports_complete_shaded_patterns_and_annotations(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);folder=root/'原图'/'007';folder.mkdir(parents=True)
            for i,gray in enumerate((shaded_trenches(),shaded_slots(),shaded_vias()),1):
                sem._imwrite_unicode_v19(folder/f'{i}.png',gray)
                (folder/f'{i}.txt').write_text('PixelSize=2')
            output=root/'结果'
            completed=subprocess.run([sys.executable,str(Path(sem.__file__)),'--input',str(folder.parent),
                                      '--output',str(output),'--image-timeout-seconds','120'],
                                     capture_output=True,text=True,timeout=300,
                                     env={**os.environ,'MPLBACKEND':'Agg','OPENBLAS_NUM_THREADS':'1','PYTHONDONTWRITEBYTECODE':'1'})
            self.assertEqual(completed.returncode,0,completed.stdout+completed.stderr)
            data=sem.pd.read_csv(output/'all_object_measurements.csv',dtype={'condition':str})
            self.assertEqual(len(data),32)
            for pattern,width,height in [('slot',34.,42.),('via',82.,82.)]:
                group=data[data.pattern==pattern]
                np.testing.assert_allclose(group.x_width_nm,width,atol=4.)
                np.testing.assert_allclose(group.y_height_nm,height,atol=4.)
            trenches=data[data.pattern=='trench'];np.testing.assert_allclose(trenches.tip_gap_y_nm,120.,atol=2.)
            self.assertTrue((data.edge_signal_source=='original_image_gray').all())
            status=sem.pd.read_csv(output/'image_status.csv')
            self.assertTrue(all(Path(p).is_file() for p in status.annotated_path))
            self.assertEqual(len(list((output/'plots').glob('*.png'))),18)
            settings=json.loads((output/'settings.json').read_text())
            self.assertIn('V37',settings['script_version'])
            self.assertTrue((output/'SEM_0830_0831_single_V37_results.xlsx').is_file())


if __name__=='__main__':
    unittest.main()
