#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SEM before/after comparison V32: default via60, image-derived via search scale, relaxed relative size gates.
V31 original-gray edge validation, four-row slots and 300s image limit are retained.

Input: numeric-named immediate subfolders discovered under the after root.
Process them in numeric order and match exact folder names under the before root.
Extra before folders and nonnumeric after folders are ignored; leading zeros stay.
Supported images: .tif, .tiff, .png, .jpg, .jpeg (case-insensitive; may be mixed).
Each condition may contain any number of images, with equal counts before/after.
Images pair by natural filename order, even when filenames or formats differ.
Patterns repeat trench, slot, via; regions start at 2 and increase every 3 images.
Every image requires a same-stem TXT containing PixelSize=<number> in nm/pixel.

Run with --before, --after and --output to select the input and output roots.
Outputs include annotated images, object-level CSV, summary CSV/Excel,
rejected object logs, and before/after comparison plots.

Dependencies:
  pip install numpy pandas scipy opencv-python matplotlib openpyxl
"""

from __future__ import annotations

import math
import re
import sys
import time
import warnings
from dataclasses import dataclass
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Suppress OpenCV warning-level log messages by default. The vendor TIFF metadata
# warning is harmless and does not mean the image failed to load.
try:
    cv2.setLogLevel(2)  # ERROR level; suppress WARNING/INFO
except Exception:
    pass


# ============================================================
# 1) USER CONFIGURATION
# ============================================================
VERSION_NUMBER = 24

# ===== V19: normally you only need to edit these two paths =====
BEFORE_ROOT_DIR_V19 = Path(r"C:\Users\z00027644\PycharmProjects\sem_compare\before")
AFTER_ROOT_DIR_V19 = Path(r"C:\Users\z00027644\PycharmProjects\sem_compare\after")
# None -> create sem_before_after_outputV24 beside BEFORE_ROOT_DIR_V19.
OUTPUT_ROOT_DIR_V19: Optional[Path] = None

# Legacy V18 constants retained because the mature helper functions below refer to
# their names. The V19 main pipeline does not use ROOT_DIR/CONDITIONS directly.
ROOT_DIR = BEFORE_ROOT_DIR_V19
CONDITIONS = [str(i) for i in range(1, 9)]
# Folder 1-8 are eight independent conditions. Specific process names are not
# assigned yet, so all user-facing labels remain generic.
CONDITION_LABELS_ZH = {str(i): f"Condition {i}" for i in range(1, 9)}
OUTPUT_DIR_NAME = f"analysis_outputV{VERSION_NUMBER}"

# Feature polarity in SEM image: usually etched holes/trenches/slots are darker.
# Options: "dark" or "bright".
FEATURE_POLARITY = "dark"

# Pattern type is determined strictly by filename prefix (case-insensitive).
# ZEP150* -> opposing trench-tip pattern
# ZEP200* -> vertical slots with different lengths
# ZEP210* -> via / round-hole pattern
PATTERN_PREFIX_MAP = {
    "ZEP150": "trench",
    "ZEP200": "slot",
    "ZEP210": "via",
}

# Optional exact-stem override. This takes priority over prefix mapping.
MANUAL_PATTERN_BY_STEM: Dict[str, str] = {}

# Do not guess morphology when filename prefix is unknown.
AUTO_CLASSIFY_IF_UNKNOWN = False

# Print detailed PixelSize diagnostics for every TIF.
PIXELSIZE_DEBUG = False
SCRIPT_VERSION = "2026-08-26-V25-ZERO-RESULT-AUTO-FALLBACK"

# Global segmentation tuning.
GAUSSIAN_BLUR_KSIZE = 5      # must be odd, 0/1 disables blur
MORPH_KERNEL = 3
MORPH_OPEN_ITER = 0  # keep tiny candidates visible; <30 px are rejected/logged later
MORPH_CLOSE_ITER = 2
MIN_COMPONENT_AREA_PX = 30
MAX_COMPONENT_AREA_FRAC = 0.25
BORDER_MARGIN_PX = 3

# Via filtering: reject obviously non-circular holes.
VIA_MIN_CIRCULARITY = 0.68
VIA_MAX_AXIS_RATIO = 1.30    # max(width,height)/min(width,height)
VIA_MIN_SOLIDITY = 0.88
VIA_MIN_DIAMETER_PX = 8

# Trench filtering and pairing.
# ZEP150 geometry assumed here:
#   one dark vertical trench above + one dark vertical trench below,
#   with a bright tip-to-tip gap between them.
# Measurement definition:
#   trench width = arithmetic mean of the upper/lower edge-detected rectangle widths
#   tip-to-tip    = distance along image Y between the lower tip of the upper
#                   trench and the upper tip of the lower trench
TRENCH_MIN_VERTICAL_ASPECT = 1.8
TRENCH_MAX_ANGLE_FROM_Y_DEG = 25.0
TRENCH_MIN_LENGTH_PX = 15
TRENCH_MAX_X_OFFSET_FACTOR = 1.75      # |center_x(top)-center_x(bottom)| / mean width
TRENCH_MIN_X_OVERLAP_FRAC = 0.30      # X-overlap / smaller trench X-span
TRENCH_MAX_GAP_FACTOR = 20.0           # Y tip-gap / mean trench width
TRENCH_MAX_WIDTH_MISMATCH_FRAC = 0.80

# Robust measurement settings. These deliberately avoid single-pixel protrusions.
TRENCH_TIP_CENTER_X_FRAC = 0.50        # use central 50% of trench width to locate each tip
TRENCH_WIDTH_BODY_Y_FRAC = 0.60        # use central 60% of trench length to estimate X width
TRENCH_MIN_TIP_SAMPLE_COLUMNS = 3
TRENCH_MIN_WIDTH_SAMPLE_ROWS = 5

# Slot filtering.
# ZEP200 already tells us the image is a slot pattern, so morphology filtering is
# intentionally permissive. Near-circular / strongly trimmed objects are allowed.
# The global connected-component pixel-area threshold is the main hard filter.
SLOT_MIN_AREA_PX = 30

# Additional slot-only relative area filter.
# After the basic >=30 px^2 screening, calculate the mean filled-pixel area
# of all remaining slot candidates in that image. Reject a candidate if:
#     area < mean_area * SLOT_AREA_MIN_MEAN_FACTOR
# or  area > mean_area * SLOT_AREA_MAX_MEAN_FACTOR
# V32 widens these relative area limits to 1/20 and 20x.
SLOT_USE_RELATIVE_AREA_FILTER = True
SLOT_AREA_MIN_MEAN_FACTOR = 0.05
SLOT_AREA_MAX_MEAN_FACTOR = 20.0

SLOT_USE_ASPECT_FILTER = False
SLOT_MIN_ASPECT = 1.0
SLOT_USE_LENGTH_FILTER = False
SLOT_MIN_LENGTH_PX = 0
SLOT_USE_ANGLE_FILTER = False
SLOT_EXPECT_VERTICAL = True
SLOT_MAX_ANGLE_FROM_Y_DEG = 45.0
SLOT_REQUIRE_AT_LEAST_N = 4

# Slot regularity QC.
# Goal: accept reasonably regular rectangles/squares/circles/ellipses, including
# strongly trimmed / near-circular slots, while rejecting jagged, branched,
# highly concave or irregular connected components. Aspect ratio is NOT used here.
SLOT_USE_REGULARITY_FILTER = True
SLOT_MIN_REGULARITY_SOLIDITY = 0.80
SLOT_MIN_CONVEXITY = 0.70       # convex-hull perimeter / actual perimeter; lower = more jagged
SLOT_MIN_TEMPLATE_IOU = 0.58     # best IoU to rotated rectangle or fitted ellipse
SLOT_TEMPLATE_PAD_PX = 6

# ------------------------------------------------------------
# V13: FAST search + equal-pitch array assistance for VIA/SLOT
# ------------------------------------------------------------
# The expensive V12 strategy (many full-image masks + O(N^2) duplicate merging)
# is replaced by:
#   1) primary Otsu segmentation;
#   2) at most a few capped full-image fallback masks;
#   3) equal-pitch grid fitting;
#   4) expensive enhancement ONLY inside missing grid-site ROIs.
USE_ENHANCED_SEARCH_V12 = True   # kept name internally for compatibility
ARRAY_ASSIST_VIA = True
ARRAY_ASSIST_SLOT = True

# Suppress harmless OpenCV/libtiff warnings about vendor-specific TIFF tags
# (for example tags 40092/40093). Set False if you want to see all OpenCV warnings.
SUPPRESS_OPENCV_WARNINGS = True

# Full-image fallback search: intentionally small and capped.
ENHANCED_CLAHE_CLIP = 2.5
ENHANCED_CLAHE_TILE = 8
ENHANCED_ADAPTIVE_BLOCK = 61   # odd number
ENHANCED_ADAPTIVE_C = 3
ENHANCED_LOCAL_SIGMAS = (12.0,)  # V13: one scale only, not 6/12/22 + percentile variants
ENHANCED_MORPH_CLOSE_ITER = 1
FULL_IMAGE_FALLBACK_MIN_PRIMARY = 6
MAX_ENHANCED_CANDIDATES_PER_MASK = 220
MAX_CANDIDATE_POOL = 320
MERGE_SPATIAL_CELL_PX = 32.0

# Limit local recovery work even if an incorrectly fitted grid becomes huge.
MAX_LOCAL_RECOVERY_SITES = 80
LOCAL_RECOVERY_PROGRESS_EVERY = 10

# Candidate merging / array fitting.
ARRAY_CENTER_CLUSTER_TOL_PX = 8.0
ARRAY_MIN_GRID_SEEDS = 4
ARRAY_MIN_AXIS_POSITIONS = 2
ARRAY_MIN_AXIS_CLUSTER_SUPPORT = 2
ARRAY_MAX_EXPECTED_OBJECTS = 240
ARRAY_GRID_RESIDUAL_FRAC = 0.22
ARRAY_GRID_MIN_OCCUPANCY = 0.28
ARRAY_MATCH_PITCH_FRAC = 0.32
ARRAY_MATCH_SIZE_FACTOR = 0.90

# Local recovery around a missing expected lattice site.
ARRAY_LOCAL_ROI_PITCH_FRAC = 0.45
ARRAY_LOCAL_ROI_SIZE_FACTOR = 2.4
ARRAY_LOCAL_CENTER_TOL_FRAC = 0.45

# Relaxed shape requirements are only used when a candidate sits on a
# well-supported array site. They are deliberately looser than normal QC.
VIA_GRID_MIN_CIRCULARITY = 0.42
VIA_GRID_MAX_AXIS_RATIO = 1.75
VIA_GRID_MIN_SOLIDITY = 0.72

SLOT_GRID_MIN_SOLIDITY = 0.62
SLOT_GRID_MIN_CONVEXITY = 0.50
SLOT_GRID_MIN_TEMPLATE_IOU = 0.34

# Diagnostics.
ARRAY_DEBUG = True
ANNOTATE_ARRAY_EXPECTED_CENTERS = False

# Annotation: draw one representative measurement example in each image.
ANNOTATE_ONE_EXAMPLE_PER_IMAGE = True
ANNOTATION_FONT_SCALE = 0.55
ANNOTATION_THICKNESS = 1

# Also outline EVERY structure that actually participates in the measurement.
# Accepted/selected objects are GREEN; rejected objects remain RED.
ANNOTATE_SELECTED_OBJECTS = True
SELECTED_CONTOUR_THICKNESS = 1
SELECTED_LABEL_FONT_SCALE = 0.42

# Set True if you want to inspect binary masks as well.
SAVE_BINARY_MASK = False

# Rejection/QC visualization and logging.
ANNOTATE_REJECTED_OBJECTS = True
ANNOTATE_SEGMENTATION_REJECTS = True
MAX_REJECTED_ANNOTATIONS_PER_IMAGE = 60
REJECTED_CONTOUR_THICKNESS = 1
REJECTED_X_THICKNESS = 1
REJECTED_LABEL_FONT_SCALE = 0.42


# ------------------------------------------------------------
# V14: arbitrary 2D lattice fitting (rectangular / rotated / rhombic)
# ------------------------------------------------------------
# A lattice is represented by two non-collinear basis vectors a and b:
#     p(i,j) = origin + i*a + j*b
# The vectors do NOT need to be horizontal/vertical or perpendicular.
LATTICE_ASSIST_VIA = True
LATTICE_ASSIST_SLOT = True
LATTICE_MIN_SEEDS = 4
LATTICE_NEAREST_NEIGHBORS = 8
LATTICE_MAX_VECTOR_CANDIDATES = 18
LATTICE_VECTOR_ANGLE_TOL_DEG = 10.0
LATTICE_VECTOR_LENGTH_REL_TOL = 0.22
LATTICE_MIN_BASIS_ANGLE_DEG = 18.0
LATTICE_MAX_BASIS_ANGLE_DEG = 162.0
LATTICE_RESIDUAL_FRAC = 0.24
LATTICE_MIN_OCCUPANCY = 0.24
LATTICE_MAX_EXPECTED_OBJECTS = 260
LATTICE_MATCH_FRAC = 0.32
LATTICE_REFINE_ITERS = 2
LATTICE_FINAL_EXTEND_RINGS = 1
LATTICE_FINAL_MAX_SITES = 100
LATTICE_FINAL_PROGRESS_EVERY = 10
LATTICE_FINAL_REQUIRE_LOCAL_CONTRAST = True
LATTICE_FINAL_MIN_CONTRAST_GRAY = 1.2

# Final-pass acceptance is intentionally looser because the position itself is
# strongly constrained by an already-fitted lattice. The hard 30-pixel area
# threshold is still respected.
VIA_FINAL_MIN_CIRCULARITY = 0.24
VIA_FINAL_MAX_AXIS_RATIO = 2.30
VIA_FINAL_MIN_SOLIDITY = 0.52
SLOT_FINAL_MIN_SOLIDITY = 0.46
SLOT_FINAL_MIN_CONVEXITY = 0.34
SLOT_FINAL_MIN_TEMPLATE_IOU = 0.20
FINAL_CENTER_TOL_FRAC = 0.46
FINAL_SLOT_AREA_MIN_FACTOR = 0.10
FINAL_SLOT_AREA_MAX_FACTOR = 10.0

# Extra local threshold percentiles used only at predicted lattice sites.
FINAL_LOCAL_PERCENTILES_DARK = (18.0, 24.0, 30.0, 36.0, 42.0)
FINAL_LOCAL_PERCENTILES_BRIGHT = (58.0, 64.0, 70.0, 76.0, 82.0)

# ------------------------------------------------------------
# V15: center-to-center recursive array completion
# ------------------------------------------------------------
# This is an EXTRA final recall pass after the V14 affine-lattice recovery.
# Predictions are made from the CENTER of every already-selected object.
# Horizontal/vertical steps are estimated from center-to-center spacing.
# For the 45-degree directions, each X/Y component is:
#     horizontal_pitch / sqrt(2)  ~= horizontal_pitch / 1.414
# so the 45-degree center-to-center vector length is approximately the
# horizontal center pitch. Every newly recovered object is pushed back into
# the search queue and can generate more neighbors until the image boundary.
CENTER_RECURSIVE_SEARCH_VIA = True
CENTER_RECURSIVE_SEARCH_SLOT = True
CENTER_DIRECTION_ANGLE_TOL_DEG = 20.0
CENTER_PITCH_MIN_PX = 8.0
CENTER_PITCH_MAX_IMAGE_FRAC = 0.75
CENTER_DIAGONAL_DIVISOR = 1.41421356237
CENTER_RECURSIVE_MAX_ATTEMPTS = 420
CENTER_RECURSIVE_MAX_RECOVERED = 180
CENTER_RECURSIVE_PROGRESS_EVERY = 20
CENTER_RECURSIVE_SITE_TOL_FRAC = 0.24
CENTER_RECURSIVE_BOUNDARY_MARGIN_PX = 4
CENTER_RECURSIVE_LOCAL_CONTRAST_GRAY = 0.45
CENTER_RECURSIVE_CENTER_TOL_FRAC = 0.58

# Very loose shape QC is allowed ONLY at a center position predicted from an
# already accepted equal-spacing pattern. The 30-pixel hard area filter remains.
VIA_CENTER_MIN_CIRCULARITY = 0.14
VIA_CENTER_MAX_AXIS_RATIO = 3.20
VIA_CENTER_MIN_SOLIDITY = 0.38
SLOT_CENTER_MIN_SOLIDITY = 0.30
SLOT_CENTER_MIN_CONVEXITY = 0.20
SLOT_CENTER_MIN_TEMPLATE_IOU = 0.10

# Extra local percentile masks used only in V15 center-predicted ROIs.
CENTER_LOCAL_PERCENTILES_DARK = (12.0, 18.0, 24.0, 30.0, 36.0, 42.0, 48.0)
CENTER_LOCAL_PERCENTILES_BRIGHT = (52.0, 58.0, 64.0, 70.0, 76.0, 82.0, 88.0)


# ------------------------------------------------------------
# V16: fixed-pitch center array completion
# ------------------------------------------------------------
# Filename rule:
#   stem contains "-1" -> search horizontal, vertical, and four 45-degree directions
#   otherwise          -> search horizontal and vertical only
# The pitch is estimated ONCE from the initial accepted object CENTERS and is frozen.
# Recovered objects NEVER trigger re-estimation or pitch subdivision.
V16_DIAGONAL_FILENAME_TOKEN = "-1"
V16_DIRECTION_ANGLE_TOL_DEG = 16.0
V16_PITCH_MIN_PX = 8.0
V16_PITCH_MAX_IMAGE_FRAC = 0.75
V16_DIAGONAL_DIVISOR = 1.41421356237
V16_RECURSIVE_MAX_ATTEMPTS = 520
V16_RECURSIVE_MAX_RECOVERED = 220
V16_RECURSIVE_PROGRESS_EVERY = 20
V16_SITE_TOL_FRAC = 0.22
V16_OCCUPIED_TOL_FRAC = 0.18
V16_BOUNDARY_MARGIN_PX = 4
V16_LOCAL_CONTRAST_GRAY = 0.30
V16_LOCAL_CENTER_TOL_FRAC = 0.48
V16_LOCAL_ROI_BORDER_MARGIN_PX = 2

# Recovery may be blurrier, but shape must remain complete and regular.
# VIA: keep an intact near-round/elliptical closed object.
VIA_V16_RECOVER_MIN_CIRCULARITY = 0.50
VIA_V16_RECOVER_MAX_AXIS_RATIO = 1.60
VIA_V16_RECOVER_MIN_SOLIDITY = 0.82

# SLOT: rectangle/square/circle/ellipse are all allowed; jagged, branched,
# strongly concave, clipped or incomplete objects are rejected.
SLOT_V16_RECOVER_MIN_SOLIDITY = 0.76
SLOT_V16_RECOVER_MIN_CONVEXITY = 0.62
SLOT_V16_RECOVER_MIN_TEMPLATE_IOU = 0.48

# Additional expected-size guard. Slot's existing 0.1x~10x image-area rule is still kept.
V16_RECOVER_MIN_EXPECTED_AREA_FACTOR = 0.22
V16_RECOVER_MAX_EXPECTED_AREA_FACTOR = 4.5

# Local threshold masks used only at a fixed predicted center. These improve blur recall;
# morphology standards above are NOT relaxed.
V16_LOCAL_PERCENTILES_DARK = (16.0, 22.0, 28.0, 34.0, 40.0, 46.0)
V16_LOCAL_PERCENTILES_BRIGHT = (54.0, 60.0, 66.0, 72.0, 78.0, 84.0)


# ------------------------------------------------------------
# V17: SLOT-SPECIFIC 4-ROW ARRAY COMPLETION
# ------------------------------------------------------------
# ZEP200 layout prior supplied by user:
#   - exactly four short-slot rows, numbered from top to bottom 1..4;
#   - each row contains roughly ten equally spaced short slots;
#   - membership is based on the object CENTER being near the row/site, not on
#     strict graph connectivity to neighboring objects;
#   - blur may be tolerated at predicted sites, but recovered contours must remain
#     complete and reasonably regular (rectangle/square/circle/ellipse-like).
V17_SLOT_EXPECTED_ROWS = 4
V17_SLOT_APPROX_OBJECTS_PER_ROW = 10
V17_SLOT_MIN_ROW_PITCH_PX = 10.0
V17_SLOT_MAX_ROW_PITCH_FRAC = 0.45
V17_SLOT_ROW_MODEL_RESIDUAL_FRAC = 0.24
V17_SLOT_ROW_CENTER_TOL_FRAC = 0.38
V17_SLOT_SITE_X_TOL_FRAC = 0.40
V17_SLOT_SITE_Y_TOL_FRAC = 0.34
V17_SLOT_OCCUPIED_TOL_FRAC = 0.23
V17_SLOT_MAX_SITES_PER_ROW = 18
V17_SLOT_MIN_EXPECTED_SITES_PER_ROW = 5
V17_SLOT_MAX_EXPECTED_SITES_PER_ROW = 16
V17_SLOT_PROGRESS_EVERY = 10

# Initial slot candidates used for fitting can be somewhat looser than the old
# all-purpose V16 array membership rule. These thresholds still reject strong
# deformation / jagged / branched contours.
V17_SLOT_SEED_MIN_SOLIDITY = 0.70
V17_SLOT_SEED_MIN_CONVEXITY = 0.54
V17_SLOT_SEED_MIN_TEMPLATE_IOU = 0.38

# Recovery at an already-predicted 4-row site may be blurry, but shape integrity
# remains mandatory. Aggressive recovery changes image thresholding, NOT geometry QC.
V17_SLOT_RECOVER_MIN_SOLIDITY = 0.70
V17_SLOT_RECOVER_MIN_CONVEXITY = 0.54
V17_SLOT_RECOVER_MIN_TEMPLATE_IOU = 0.38
V17_SLOT_RECOVER_LOCAL_CONTRAST_GRAY = 0.10
V17_SLOT_RECOVER_CENTER_X_FRAC = 0.46
V17_SLOT_RECOVER_CENTER_Y_FRAC = 0.42
V17_SLOT_RECOVER_ROI_BORDER_MARGIN_PX = 1
V17_SLOT_LOCAL_PERCENTILES_DARK = (10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0)
V17_SLOT_LOCAL_PERCENTILES_BRIGHT = (50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0)

# If a row is weak, run the local masks again without the local-contrast pre-gate.
# The final geometry thresholds above are unchanged, so this raises recall for blur
# without accepting incomplete/deformed objects.
V17_SLOT_AGGRESSIVE_RECOVERY_FOR_MISSING_ROWS = True
V17_SLOT_ROW_WEAK_COUNT = 5


# ============================================================
# 2) DATA STRUCTURES
# ============================================================
@dataclass
class ShapeDesc:
    contour: np.ndarray
    area: float
    perimeter: float
    circularity: float
    solidity: float
    center: np.ndarray     # [x, y]
    u: np.ndarray          # major-axis unit vector [x, y]
    v: np.ndarray          # minor-axis unit vector [x, y]
    major_min: float
    major_max: float
    minor_min: float
    minor_max: float
    length_px: float
    width_px: float
    angle_from_x_deg: float

    @property
    def aspect(self) -> float:
        return self.length_px / max(self.width_px, 1e-9)


# ============================================================
# 3) IO AND BASIC IMAGE PROCESSING
# ============================================================
def read_gray_image(path: Path) -> np.ndarray:
    img = cv2.imread(str(path), cv2.IMREAD_UNCHANGED)
    if img is None:
        raise RuntimeError(f"Cannot read image: {path}")

    if img.ndim == 3:
        if img.shape[2] == 4:
            img = cv2.cvtColor(img, cv2.COLOR_BGRA2GRAY)
        else:
            img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    # Normalize any bit depth to uint8 for thresholding/display.
    img = img.astype(np.float32)
    lo, hi = np.percentile(img, [1.0, 99.0])
    if hi <= lo:
        lo, hi = float(img.min()), float(img.max())
    if hi <= lo:
        return np.zeros_like(img, dtype=np.uint8)
    img8 = np.clip((img - lo) * 255.0 / (hi - lo), 0, 255).astype(np.uint8)
    return img8


class PixelSizeTxtNotFound(FileNotFoundError):
    """The expected same-stem TXT does not exist."""


class PixelSizeFieldNotFound(ValueError):
    """TXT exists/readable, but no PixelSize field could be parsed."""


@dataclass
class PixelSizeInfo:
    value_nm: float
    txt_path: Path
    encoding: str
    matched_line: str


_PIXEL_RE = re.compile(
    r"Pix(?:el|cel)\s*[_ ]?\s*Size\s*=\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)",
    flags=re.IGNORECASE,
)


def _decode_txt_variants(raw: bytes) -> List[Tuple[str, str]]:
    """Return distinct decoded variants. UTF-16 is important for many SEM metadata files."""
    variants: List[Tuple[str, str]] = []
    seen_text = set()
    encodings = [
        "utf-8-sig",
        "utf-16",
        "utf-16-le",
        "utf-16-be",
        "gb18030",
        "cp936",
        "latin-1",
    ]
    for enc in encodings:
        try:
            text = raw.decode(enc)
        except Exception:
            continue
        if text not in seen_text:
            variants.append((enc, text))
            seen_text.add(text)

    # Last-resort form for ASCII metadata stored with interleaved NUL bytes.
    try:
        nul_stripped = raw.replace(b"\x00", b"").decode("latin-1", errors="ignore")
        if nul_stripped not in seen_text:
            variants.append(("nul-stripped-latin1", nul_stripped))
    except Exception:
        pass
    return variants


def _find_same_stem_txt(tif_path: Path) -> Tuple[Optional[Path], List[Path]]:
    """Find exactly the TXT whose filename stem equals the image stem, case-insensitively."""
    folder = tif_path.parent
    expected = tif_path.with_suffix(".txt")

    # List all files because Windows may hide extensions/case differences.
    all_files = [p for p in folder.iterdir() if p.is_file()]
    txt_like = [p for p in all_files if p.suffix.lower() == ".txt"]
    matches = [p for p in txt_like if p.stem.lower() == tif_path.stem.lower()]

    if matches:
        # Prefer exact expected path if it exists, otherwise first case-insensitive match.
        for m in matches:
            if m.name == expected.name:
                return m, txt_like
        return sorted(matches, key=lambda x: x.name.lower())[0], txt_like
    return None, txt_like


def read_pixel_size_nm(tif_path: Path) -> PixelSizeInfo:
    """
    Read PixelSize from the exact same-stem TXT only.

    Example:
      ZEP210-160-1.tif -> ZEP210-160-1.txt
      TXT line: PixelSize=0.8268229

    Errors are intentionally split into:
      1) PixelSizeTxtNotFound  -> same-stem TXT file was not found
      2) PixelSizeFieldNotFound -> TXT exists but PixelSize could not be parsed
    """
    txt_path, txt_files = _find_same_stem_txt(tif_path)
    expected = tif_path.with_suffix(".txt")

    if PIXELSIZE_DEBUG:
        print(f"\n    [PixelSize DEBUG]")
        print(f"      Image absolute path   : {tif_path.resolve()}")
        print(f"      Expected TXT path    : {expected.resolve()}")
        print(f"      Expected TXT exists? : {expected.exists()}")
        print(f"      Folder exists?       : {tif_path.parent.exists()}")
        print(f"      TXT files in folder  : {len(txt_files)}")
        if txt_files:
            for f in txt_files[:30]:
                print(f"        - {f.name}")
            if len(txt_files) > 30:
                print(f"        ... and {len(txt_files)-30} more")

    if txt_path is None:
        # Extra diagnostic: maybe Windows Explorer hid a second extension, e.g. xxx.txt.txt
        similarly_named = [
            p for p in tif_path.parent.iterdir()
            if p.is_file() and p.name.lower().startswith(tif_path.stem.lower())
        ]
        extra = ""
        if similarly_named:
            extra = "\n      Files starting with same stem: " + ", ".join(p.name for p in similarly_named)
        raise PixelSizeTxtNotFound(
            f"Same-stem TXT FILE NOT FOUND for {tif_path.name}.\n"
            f"      Expected: {expected}\n"
            f"      Image folder: {tif_path.parent}\n"
            f"      Found TXT files: {[p.name for p in txt_files]}"
            f"{extra}"
        )

    if PIXELSIZE_DEBUG:
        print(f"      Matched TXT path     : {txt_path.resolve()}")
        print(f"      Matched TXT size     : {txt_path.stat().st_size} bytes")

    try:
        raw = txt_path.read_bytes()
    except Exception as e:
        raise PixelSizeFieldNotFound(
            f"TXT FILE FOUND but could not be read: {txt_path}\n"
            f"      Read error: {type(e).__name__}: {e}"
        ) from e

    variants = _decode_txt_variants(raw)
    previews: List[str] = []
    for enc, text in variants:
        # Normalize BOM/NULs just in case, but preserve line structure for reporting.
        normalized = text.replace("\ufeff", "").replace("\x00", "")
        m = _PIXEL_RE.search(normalized)

        # Also collect candidate lines containing 'pixel' for diagnostics.
        candidate_lines = [ln.strip() for ln in normalized.splitlines() if "pixel" in ln.lower()]
        preview = " | ".join(candidate_lines[:5])
        if not preview:
            preview = normalized[:300].replace("\r", " ").replace("\n", " | ")
        previews.append(f"[{enc}] {preview}")

        if m:
            value = float(m.group(1))
            if not np.isfinite(value) or value <= 0:
                raise PixelSizeFieldNotFound(
                    f"TXT FILE FOUND and PixelSize field matched, but value is invalid: {m.group(1)!r}\n"
                    f"      TXT: {txt_path}\n"
                    f"      Encoding: {enc}"
                )

            matched_line = ""
            for ln in normalized.splitlines():
                if m.group(0) in ln:
                    matched_line = ln.strip()
                    break
            if not matched_line:
                matched_line = m.group(0)

            if PIXELSIZE_DEBUG:
                print(f"      TXT encoding used    : {enc}")
                print(f"      Matched metadata line: {matched_line}")
                print(f"      Parsed PixelSize     : {value:.10g} nm/pixel")

            return PixelSizeInfo(
                value_nm=value,
                txt_path=txt_path,
                encoding=enc,
                matched_line=matched_line,
            )

    # TXT exists, was readable, but no field matched in any decoding.
    diagnostic = "\n".join("      " + x for x in previews[:8])
    raise PixelSizeFieldNotFound(
        f"TXT FILE FOUND, but PixelSize FIELD NOT FOUND for {tif_path.name}.\n"
        f"      TXT: {txt_path}\n"
        f"      Expected format example: PixelSize=0.8268229\n"
        f"      Decoded previews:\n{diagnostic}"
    )


def segment_features(gray: np.ndarray) -> np.ndarray:
    img = gray.copy()
    if GAUSSIAN_BLUR_KSIZE and GAUSSIAN_BLUR_KSIZE >= 3:
        k = GAUSSIAN_BLUR_KSIZE if GAUSSIAN_BLUR_KSIZE % 2 == 1 else GAUSSIAN_BLUR_KSIZE + 1
        img = cv2.GaussianBlur(img, (k, k), 0)

    if FEATURE_POLARITY.lower() == "dark":
        _, mask = cv2.threshold(img, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
    elif FEATURE_POLARITY.lower() == "bright":
        _, mask = cv2.threshold(img, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    else:
        raise ValueError("FEATURE_POLARITY must be 'dark' or 'bright'")

    k = max(1, int(MORPH_KERNEL))
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k, k))
    if MORPH_OPEN_ITER > 0:
        mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel, iterations=MORPH_OPEN_ITER)
    if MORPH_CLOSE_ITER > 0:
        mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, kernel, iterations=MORPH_CLOSE_ITER)
    return mask


def describe_contour(c: np.ndarray) -> Optional[ShapeDesc]:
    pts = c.reshape(-1, 2).astype(np.float64)
    if len(pts) < 5:
        return None

    area = float(cv2.contourArea(c))
    perim = float(cv2.arcLength(c, True))
    if perim <= 0:
        return None
    circularity = float(4.0 * math.pi * area / (perim * perim))

    hull = cv2.convexHull(c)
    hull_area = float(cv2.contourArea(hull))
    solidity = area / hull_area if hull_area > 0 else 0.0

    center = pts.mean(axis=0)
    centered = pts - center
    cov = np.cov(centered.T)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    u = eigvecs[:, order[0]].astype(np.float64)
    u /= np.linalg.norm(u) + 1e-12
    # Fix sign for reproducibility.
    if u[0] < 0 or (abs(u[0]) < 1e-12 and u[1] < 0):
        u = -u
    v = np.array([-u[1], u[0]], dtype=np.float64)

    p_major = pts @ u
    p_minor = pts @ v
    major_min, major_max = float(p_major.min()), float(p_major.max())
    minor_min, minor_max = float(p_minor.min()), float(p_minor.max())
    length_px = major_max - major_min
    width_px = minor_max - minor_min
    if width_px > length_px:
        # Swap axes so u is always major axis.
        u, v = v, -u
        p_major = pts @ u
        p_minor = pts @ v
        major_min, major_max = float(p_major.min()), float(p_major.max())
        minor_min, minor_max = float(p_minor.min()), float(p_minor.max())
        length_px = major_max - major_min
        width_px = minor_max - minor_min

    angle = math.degrees(math.atan2(u[1], u[0]))
    # Reduce to [0, 90] angle from x-axis.
    angle_abs = abs(angle) % 180.0
    if angle_abs > 90.0:
        angle_abs = 180.0 - angle_abs

    return ShapeDesc(
        contour=c,
        area=area,
        perimeter=perim,
        circularity=circularity,
        solidity=solidity,
        center=center,
        u=u,
        v=v,
        major_min=major_min,
        major_max=major_max,
        minor_min=minor_min,
        minor_max=minor_max,
        length_px=float(length_px),
        width_px=float(width_px),
        angle_from_x_deg=float(angle_abs),
    )


def contour_pixel_count(contour: np.ndarray) -> int:
    """Count filled foreground pixels inside a contour (true pixel-area QC)."""
    if contour is None or len(contour) == 0:
        return 0
    x, y, w, h = cv2.boundingRect(contour.astype(np.int32))
    if w <= 0 or h <= 0:
        return 0
    local = np.zeros((h, w), dtype=np.uint8)
    c = contour.astype(np.int32).copy()
    c[:, 0, 0] -= x
    c[:, 0, 1] -= y
    cv2.drawContours(local, [c], -1, 255, thickness=cv2.FILLED)
    return int(cv2.countNonZero(local))


def _make_reject_from_contour(
    contour: np.ndarray,
    stage: str,
    reason_code: str,
    reason: str,
    desc: Optional[ShapeDesc] = None,
    **metrics,
) -> dict:
    """Internal rejection record. Private contour/descriptor fields are removed before CSV export."""
    if contour is not None and len(contour):
        pts = contour.reshape(-1, 2).astype(float)
        cx = float(np.mean(pts[:, 0]))
        cy = float(np.mean(pts[:, 1]))
        area = float(cv2.contourArea(contour))
    else:
        cx = cy = area = np.nan
    rec = {
        "stage": stage,
        "reason_code": reason_code,
        "reason": reason,
        "center_x_px": cx,
        "center_y_px": cy,
        "area_px2": area,
        "_contour": contour,
        "_desc": desc,
    }
    rec.update(metrics)
    return rec


def extract_descriptors_with_rejections(mask: np.ndarray) -> Tuple[List[ShapeDesc], List[dict]]:
    """
    Extract usable connected components and retain anything discarded by the
    global segmentation QC. This makes formerly silent rejections traceable.
    """
    h, w = mask.shape
    max_area = MAX_COMPONENT_AREA_FRAC * h * w
    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

    descs: List[ShapeDesc] = []
    rejected: List[dict] = []

    for c in contours:
        area = float(cv2.contourArea(c))
        pixel_area = contour_pixel_count(c)
        x, y, cw, ch = cv2.boundingRect(c)

        if pixel_area < MIN_COMPONENT_AREA_PX:
            rejected.append(_make_reject_from_contour(
                c, "segmentation_qc", "area_too_small",
                f"filled pixel area {pixel_area} px < minimum {MIN_COMPONENT_AREA_PX} px",
                bbox_x=x, bbox_y=y, bbox_w=cw, bbox_h=ch, pixel_area_px=pixel_area,
            ))
            continue

        if area > max_area:
            rejected.append(_make_reject_from_contour(
                c, "segmentation_qc", "area_too_large",
                f"component area {area:.1f} px^2 > maximum {max_area:.1f} px^2",
                bbox_x=x, bbox_y=y, bbox_w=cw, bbox_h=ch, pixel_area_px=pixel_area,
            ))
            continue

        touches_border = (
            x <= BORDER_MARGIN_PX
            or y <= BORDER_MARGIN_PX
            or x + cw >= w - BORDER_MARGIN_PX
            or y + ch >= h - BORDER_MARGIN_PX
        )
        if touches_border:
            rejected.append(_make_reject_from_contour(
                c, "segmentation_qc", "touches_image_border",
                f"component touches/crosses image border margin ({BORDER_MARGIN_PX} px)",
                bbox_x=x, bbox_y=y, bbox_w=cw, bbox_h=ch,
            ))
            continue

        d = describe_contour(c)
        if d is None:
            rejected.append(_make_reject_from_contour(
                c, "segmentation_qc", "descriptor_failed",
                "contour could not produce a stable shape descriptor",
                bbox_x=x, bbox_y=y, bbox_w=cw, bbox_h=ch,
            ))
            continue
        descs.append(d)

    return descs, rejected


def all_descriptors(mask: np.ndarray) -> List[ShapeDesc]:
    """Backward-compatible helper; main() uses extract_descriptors_with_rejections()."""
    descs, _ = extract_descriptors_with_rejections(mask)
    return descs


# ============================================================
# 4) PATTERN CLASSIFICATION
# ============================================================
def classify_by_name(tif_path: Path, condition_dir: Path) -> Optional[str]:
    stem = tif_path.stem
    if stem in MANUAL_PATTERN_BY_STEM:
        p = MANUAL_PATTERN_BY_STEM[stem].lower().strip()
        if p not in {"via", "trench", "slot"}:
            raise ValueError(f"Bad manual pattern type for {stem}: {p}")
        return p

    name_upper = tif_path.name.upper()
    for prefix, ptype in PATTERN_PREFIX_MAP.items():
        if name_upper.startswith(prefix.upper()):
            return ptype
    return None


def auto_classify(descs: List[ShapeDesc]) -> str:
    round_count = sum(
        1 for d in descs
        if d.circularity >= VIA_MIN_CIRCULARITY
        and d.aspect <= VIA_MAX_AXIS_RATIO
        and d.solidity >= VIA_MIN_SOLIDITY
        and min(d.length_px, d.width_px) >= VIA_MIN_DIAMETER_PX
    )
    elongated = [d for d in descs if d.aspect >= SLOT_MIN_ASPECT and d.length_px >= SLOT_MIN_LENGTH_PX]

    if round_count >= 1 and round_count >= len(elongated):
        return "via"

    # Slot images usually contain many repeated elongated objects (>=4, often much more).
    if len(elongated) >= SLOT_REQUIRE_AT_LEAST_N:
        vertical = [d for d in elongated if abs(90.0 - d.angle_from_x_deg) <= SLOT_MAX_ANGLE_FROM_Y_DEG]
        if not SLOT_EXPECT_VERTICAL or len(vertical) >= SLOT_REQUIRE_AT_LEAST_N:
            return "slot"

    if len(elongated) >= 2:
        return "trench"

    return "unknown"


# ============================================================
# 5) VIA MEASUREMENT
# ============================================================
def measure_via(
    descs: List[ShapeDesc], px_nm: float, image_name: str, condition: str
) -> Tuple[List[dict], Optional[ShapeDesc], List[dict]]:
    accepted: List[ShapeDesc] = []
    rejected: List[dict] = []
    rows: List[dict] = []

    for d in descs:
        axis_ratio = d.aspect
        reasons = []
        codes = []
        if d.circularity < VIA_MIN_CIRCULARITY:
            codes.append("via_low_circularity")
            reasons.append(f"circularity {d.circularity:.3f} < {VIA_MIN_CIRCULARITY:.3f}")
        if axis_ratio > VIA_MAX_AXIS_RATIO:
            codes.append("via_axis_ratio_too_high")
            reasons.append(f"axis ratio {axis_ratio:.3f} > {VIA_MAX_AXIS_RATIO:.3f}")
        if d.solidity < VIA_MIN_SOLIDITY:
            codes.append("via_low_solidity")
            reasons.append(f"solidity {d.solidity:.3f} < {VIA_MIN_SOLIDITY:.3f}")
        if min(d.length_px, d.width_px) < VIA_MIN_DIAMETER_PX:
            codes.append("via_too_small")
            reasons.append(
                f"minimum PCA diameter {min(d.length_px, d.width_px):.1f} px < {VIA_MIN_DIAMETER_PX} px"
            )

        if reasons:
            rejected.append(_make_reject_from_contour(
                d.contour,
                "via_shape_filter",
                "+".join(codes),
                "; ".join(reasons),
                desc=d,
                circularity=d.circularity,
                solidity=d.solidity,
                axis_ratio=axis_ratio,
                length_px=d.length_px,
                width_px=d.width_px,
            ))
            continue

        pts = d.contour.reshape(-1, 2).astype(np.float64)
        width_px_xy = float(pts[:, 0].max() - pts[:, 0].min())
        height_px_xy = float(pts[:, 1].max() - pts[:, 1].min())
        width_nm = width_px_xy * px_nm
        height_nm = height_px_xy * px_nm
        rows.append({
            "condition": condition,
            "image": image_name,
            "pattern": "via",
            "object_id": len(rows) + 1,
            "width_nm": width_nm,
            "height_nm": height_nm,
            "height_width_ratio": height_nm / max(width_nm, 1e-12),
            "circularity": d.circularity,
            "solidity": d.solidity,
            "axis_ratio": axis_ratio,
            "pixel_size_nm": px_nm,
        })
        accepted.append(d)

    representative = None
    if accepted:
        representative = max(accepted, key=lambda d: (d.circularity, -d.aspect))
    return rows, representative, rejected


# ============================================================
# 6) TRENCH MEASUREMENT
# ============================================================
def _contour_xy_bounds(desc: ShapeDesc) -> Tuple[float, float, float, float]:
    pts = desc.contour.reshape(-1, 2).astype(np.float64)
    return (
        float(pts[:, 0].min()), float(pts[:, 0].max()),
        float(pts[:, 1].min()), float(pts[:, 1].max()),
    )


def _filled_component_local(desc: ShapeDesc) -> Tuple[np.ndarray, int, int]:
    x, y, w, h = cv2.boundingRect(desc.contour)
    pad = 2
    local = np.zeros((h + 2 * pad, w + 2 * pad), dtype=np.uint8)
    c = desc.contour.copy().astype(np.int32)
    c[:, 0, 0] -= (x - pad)
    c[:, 0, 1] -= (y - pad)
    cv2.drawContours(local, [c], -1, 255, thickness=cv2.FILLED)
    return local, x - pad, y - pad


def _vertical_candidate_evaluate(desc: ShapeDesc) -> Tuple[Optional[dict], List[str], List[str]]:
    """Return geometry plus human-readable reasons if a component is not a vertical trench candidate."""
    x0, x1, y0, y1 = _contour_xy_bounds(desc)
    xspan = x1 - x0
    yspan = y1 - y0
    reasons: List[str] = []
    codes: List[str] = []

    if xspan <= 0 or yspan <= 0:
        return None, ["zero/invalid X or Y span"], ["trench_invalid_span"]

    vertical_aspect = yspan / max(xspan, 1e-9)
    angle_from_y = abs(90.0 - desc.angle_from_x_deg)
    if vertical_aspect < TRENCH_MIN_VERTICAL_ASPECT:
        codes.append("trench_aspect_too_low")
        reasons.append(f"vertical aspect {vertical_aspect:.2f} < {TRENCH_MIN_VERTICAL_ASPECT:.2f}")
    if yspan < TRENCH_MIN_LENGTH_PX:
        codes.append("trench_too_short")
        reasons.append(f"Y length {yspan:.1f} px < {TRENCH_MIN_LENGTH_PX} px")
    if angle_from_y > TRENCH_MAX_ANGLE_FROM_Y_DEG:
        codes.append("trench_not_vertical")
        reasons.append(f"angle from Y {angle_from_y:.1f} deg > {TRENCH_MAX_ANGLE_FROM_Y_DEG:.1f} deg")

    geom = {
        "desc": desc,
        "x0": x0, "x1": x1, "y0": y0, "y1": y1,
        "cx": 0.5 * (x0 + x1), "cy": 0.5 * (y0 + y1),
        "xspan": xspan, "yspan": yspan,
        "vertical_aspect": vertical_aspect,
        "angle_from_y_deg": angle_from_y,
    }
    return geom, reasons, codes


def _vertical_candidate_metrics(desc: ShapeDesc) -> Optional[dict]:
    geom, reasons, _ = _vertical_candidate_evaluate(desc)
    return geom if not reasons else None


def _robust_x_width(desc: ShapeDesc) -> Tuple[float, float, float, float]:
    local, ox, oy = _filled_component_local(desc)
    ys, xs = np.where(local > 0)
    if len(xs) == 0:
        return 0.0, float(desc.center[1]), float(desc.center[0]), float(desc.center[0])

    y_min, y_max = int(ys.min()), int(ys.max())
    height = max(1, y_max - y_min)
    body_frac = float(np.clip(TRENCH_WIDTH_BODY_Y_FRAC, 0.10, 1.0))
    margin = 0.5 * (1.0 - body_frac)
    ya = int(round(y_min + margin * height))
    yb = int(round(y_max - margin * height))
    if yb <= ya:
        ya, yb = y_min, y_max

    row_data = []
    for yy in range(ya, yb + 1):
        xx = np.where(local[yy] > 0)[0]
        if len(xx) >= 2:
            row_data.append((float(xx.max() - xx.min()), yy, float(xx.min()), float(xx.max())))

    if len(row_data) < TRENCH_MIN_WIDTH_SAMPLE_ROWS:
        row_data = []
        for yy in range(y_min, y_max + 1):
            xx = np.where(local[yy] > 0)[0]
            if len(xx) >= 2:
                row_data.append((float(xx.max() - xx.min()), yy, float(xx.min()), float(xx.max())))

    if not row_data:
        pts = desc.contour.reshape(-1, 2).astype(np.float64)
        width = float(pts[:, 0].max() - pts[:, 0].min())
        return width, float(desc.center[1]), float(pts[:, 0].min()), float(pts[:, 0].max())

    widths = np.array([r[0] for r in row_data], dtype=float)
    med = float(np.median(widths))
    r = min(row_data, key=lambda z: abs(z[0] - med))
    return med, float(r[1] + oy), float(r[2] + ox), float(r[3] + ox)


def _robust_tip_y(desc: ShapeDesc, facing: str, common_x: float, sample_width_px: float) -> float:
    local, ox, oy = _filled_component_local(desc)
    x0g, x1g, _, _ = _contour_xy_bounds(desc)
    half = max(1.0, 0.5 * sample_width_px)
    xa_g = max(x0g, common_x - half)
    xb_g = min(x1g, common_x + half)
    xa = int(math.ceil(xa_g - ox))
    xb = int(math.floor(xb_g - ox))

    tip_samples = []
    for xx in range(max(0, xa), min(local.shape[1] - 1, xb) + 1):
        yy = np.where(local[:, xx] > 0)[0]
        if len(yy) == 0:
            continue
        tip_samples.append(float((yy.max() if facing == "down" else yy.min()) + oy))

    if len(tip_samples) >= TRENCH_MIN_TIP_SAMPLE_COLUMNS:
        return float(np.median(tip_samples))
    pts = desc.contour.reshape(-1, 2).astype(np.float64)
    return float(np.percentile(pts[:, 1], 99.0 if facing == "down" else 1.0))


def trench_vertical_pair_evaluate(a: ShapeDesc, b: ShapeDesc) -> Tuple[Optional[dict], List[str]]:
    """Evaluate a possible upper/lower pair and explain why it fails."""
    ga = _vertical_candidate_metrics(a)
    gb = _vertical_candidate_metrics(b)
    if ga is None or gb is None:
        return None, ["one or both components are not valid vertical trench candidates"]

    top, bottom = (ga, gb) if ga["cy"] < gb["cy"] else (gb, ga)
    reasons: List[str] = []

    overlap = max(0.0, min(top["x1"], bottom["x1"]) - max(top["x0"], bottom["x0"]))
    overlap_frac = overlap / max(min(top["xspan"], bottom["xspan"]), 1e-9)
    if overlap_frac < TRENCH_MIN_X_OVERLAP_FRAC:
        reasons.append(f"X overlap {overlap_frac:.2f} < {TRENCH_MIN_X_OVERLAP_FRAC:.2f}")

    top_w, top_w_y, top_xl, top_xr = _robust_x_width(top["desc"])
    bot_w, bot_w_y, bot_xl, bot_xr = _robust_x_width(bottom["desc"])
    mean_width = 0.5 * (top_w + bot_w)
    if mean_width <= 0:
        reasons.append("robust trench width <= 0")
        return None, reasons

    x_offset = abs(top["cx"] - bottom["cx"])
    if x_offset > TRENCH_MAX_X_OFFSET_FACTOR * mean_width:
        reasons.append(
            f"X center offset {x_offset:.1f}px > {TRENCH_MAX_X_OFFSET_FACTOR:.2f} x mean width"
        )

    width_mismatch = abs(top_w - bot_w) / max(mean_width, 1e-9)
    if width_mismatch > TRENCH_MAX_WIDTH_MISMATCH_FRAC:
        reasons.append(
            f"upper/lower width mismatch {width_mismatch:.2f} > {TRENCH_MAX_WIDTH_MISMATCH_FRAC:.2f}"
        )

    common_x = 0.5 * (top["cx"] + bottom["cx"])
    tip_sample_width = max(3.0, TRENCH_TIP_CENTER_X_FRAC * mean_width)
    top_tip_y = _robust_tip_y(top["desc"], "down", common_x, tip_sample_width)
    bottom_tip_y = _robust_tip_y(bottom["desc"], "up", common_x, tip_sample_width)
    gap = bottom_tip_y - top_tip_y
    if gap <= 0:
        reasons.append(f"tip gap is non-positive ({gap:.1f}px): components overlap/cross in Y")
    elif gap > TRENCH_MAX_GAP_FACTOR * mean_width:
        reasons.append(
            f"tip gap {gap:.1f}px > {TRENCH_MAX_GAP_FACTOR:.1f} x mean trench width"
        )

    if reasons:
        return None, reasons

    score = (
        2.0 * x_offset / max(mean_width, 1e-9)
        + width_mismatch
        + 0.75 * (1.0 - overlap_frac)
        + 0.25 * gap / max(TRENCH_MAX_GAP_FACTOR * mean_width, 1e-9)
        + 0.25 * (top["angle_from_y_deg"] + bottom["angle_from_y_deg"]) /
          max(2.0 * TRENCH_MAX_ANGLE_FROM_Y_DEG, 1e-9)
    )

    return {
        "top": top["desc"], "bottom": bottom["desc"],
        "top_geom": top, "bottom_geom": bottom,
        "top_width_px": top_w, "bottom_width_px": bot_w, "width_px": mean_width,
        "top_width_y": top_w_y, "top_width_xl": top_xl, "top_width_xr": top_xr,
        "bottom_width_y": bot_w_y, "bottom_width_xl": bot_xl, "bottom_width_xr": bot_xr,
        "common_x": common_x, "top_tip_y": top_tip_y, "bottom_tip_y": bottom_tip_y,
        "gap_px": gap, "x_offset_px": x_offset, "x_overlap_frac": overlap_frac,
        "width_mismatch_frac": width_mismatch, "score": score,
    }, []


def trench_vertical_pair_metrics(a: ShapeDesc, b: ShapeDesc) -> Optional[dict]:
    m, _ = trench_vertical_pair_evaluate(a, b)
    return m


def find_trench_pairs_with_rejections(descs: List[ShapeDesc]) -> Tuple[List[dict], List[dict]]:
    valid_cand: List[ShapeDesc] = []
    rejected: List[dict] = []

    for d in descs:
        geom, reasons, codes = _vertical_candidate_evaluate(d)
        if reasons:
            rejected.append(_make_reject_from_contour(
                d.contour, "trench_candidate_filter", "+".join(codes), "; ".join(reasons), desc=d,
                vertical_aspect=(geom or {}).get("vertical_aspect", np.nan),
                angle_from_y_deg=(geom or {}).get("angle_from_y_deg", np.nan),
                y_length_px=(geom or {}).get("yspan", np.nan),
            ))
        else:
            valid_cand.append(d)

    potential = []
    compatible_partner_count = [0] * len(valid_cand)
    failure_examples: List[List[str]] = [[] for _ in valid_cand]
    for i in range(len(valid_cand)):
        for j in range(i + 1, len(valid_cand)):
            m, pair_reasons = trench_vertical_pair_evaluate(valid_cand[i], valid_cand[j])
            if m is not None:
                m["i"] = i
                m["j"] = j
                potential.append(m)
                compatible_partner_count[i] += 1
                compatible_partner_count[j] += 1
            else:
                short = "; ".join(pair_reasons[:2])
                if short:
                    failure_examples[i].append(short)
                    failure_examples[j].append(short)

    potential.sort(key=lambda x: x["score"])
    used = set()
    pairs = []
    for m in potential:
        if m["i"] in used or m["j"] in used:
            continue
        used.add(m["i"])
        used.add(m["j"])
        pairs.append(m)
    pairs.sort(key=lambda p: p["common_x"])

    for i, d in enumerate(valid_cand):
        if i in used:
            continue
        if compatible_partner_count[i] > 0:
            code = "trench_pair_conflict"
            reason = (
                f"valid vertical trench candidate had {compatible_partner_count[i]} compatible partner(s), "
                "but another lower-score pair was selected first"
            )
        else:
            code = "trench_no_compatible_partner"
            example = failure_examples[i][0] if failure_examples[i] else "no opposite vertical candidate"
            reason = f"valid vertical trench candidate could not form an upper/lower pair; example failure: {example}"
        rejected.append(_make_reject_from_contour(
            d.contour, "trench_pairing", code, reason, desc=d,
            compatible_partner_count=compatible_partner_count[i],
        ))

    return pairs, rejected


def find_trench_pairs(descs: List[ShapeDesc]) -> List[dict]:
    pairs, _ = find_trench_pairs_with_rejections(descs)
    return pairs


def measure_trench(
    descs: List[ShapeDesc], px_nm: float, image_name: str, condition: str
) -> Tuple[List[dict], Optional[dict], List[dict]]:
    pairs, rejected = find_trench_pairs_with_rejections(descs)
    rows = []
    for k, p in enumerate(pairs, 1):
        rows.append({
            "condition": condition,
            "image": image_name,
            "pattern": "trench",
            "object_id": k,
            "trench_width_nm": p["width_px"] * px_nm,
            "top_trench_width_nm": p["top_width_px"] * px_nm,
            "bottom_trench_width_nm": p["bottom_width_px"] * px_nm,
            "tip_to_tip_nm": p["gap_px"] * px_nm,
            "pair_x_offset_nm": p["x_offset_px"] * px_nm,
            "pair_x_overlap_frac": p["x_overlap_frac"],
            "width_mismatch_frac": p["width_mismatch_frac"],
            "pixel_size_nm": px_nm,
        })
    representative = min(pairs, key=lambda p: p["score"]) if pairs else None
    return rows, representative, rejected


# ============================================================
# 7) SLOT MEASUREMENT
# ============================================================
def cluster_1d_into_four_rows(y_values: np.ndarray) -> np.ndarray:
    y = np.asarray(y_values, dtype=np.float32).reshape(-1, 1)
    if len(y) < 4:
        raise ValueError("Need at least 4 slot objects for row clustering")
    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 200, 1e-3)
    _, labels, centers = cv2.kmeans(y, 4, None, criteria, 20, cv2.KMEANS_PP_CENTERS)
    centers = centers.ravel()
    order = np.argsort(centers)
    remap = {int(old): int(new + 1) for new, old in enumerate(order)}
    return np.array([remap[int(v)] for v in labels.ravel()], dtype=int)


def _binary_iou(a: np.ndarray, b: np.ndarray) -> float:
    aa = a > 0
    bb = b > 0
    inter = int(np.count_nonzero(aa & bb))
    union = int(np.count_nonzero(aa | bb))
    return float(inter / union) if union > 0 else 0.0


def slot_regularity_metrics(desc: ShapeDesc) -> Dict[str, float]:
    """
    Shape-regularity metrics for ZEP200 slot candidates.

    These metrics deliberately do not care whether the object is long, square,
    circular, or elliptical. They only ask whether it is a reasonably clean,
    compact geometric object instead of a jagged/branched/irregular blob.

    Metrics:
      solidity   = contour area / convex-hull area
      convexity  = convex-hull perimeter / contour perimeter
                   (falls as the contour gets saw-toothed / wavy)
      rect_iou   = pixel-mask IoU with its best rotated rectangle
      ellipse_iou= pixel-mask IoU with fitted ellipse (if fit is possible)
      template_iou = max(rect_iou, ellipse_iou)
    """
    c = desc.contour.astype(np.int32)
    if c is None or len(c) < 3:
        return {
            "regularity_solidity": 0.0,
            "regularity_convexity": 0.0,
            "regularity_rect_iou": 0.0,
            "regularity_ellipse_iou": 0.0,
            "regularity_template_iou": 0.0,
        }

    hull = cv2.convexHull(c)
    hull_perim = float(cv2.arcLength(hull, True))
    contour_perim = float(cv2.arcLength(c, True))
    convexity = hull_perim / contour_perim if contour_perim > 0 else 0.0

    x, y, w, h = cv2.boundingRect(c)
    pad = int(max(2, SLOT_TEMPLATE_PAD_PX))
    x0, y0 = x - pad, y - pad
    canvas_h = max(1, h + 2 * pad)
    canvas_w = max(1, w + 2 * pad)

    shape_mask = np.zeros((canvas_h, canvas_w), dtype=np.uint8)
    cl = c.copy()
    cl[:, 0, 0] -= x0
    cl[:, 0, 1] -= y0
    cv2.drawContours(shape_mask, [cl], -1, 255, thickness=cv2.FILLED)

    # Best rotated rectangle template.
    rect_mask = np.zeros_like(shape_mask)
    rect = cv2.minAreaRect(c)
    box = cv2.boxPoints(rect)
    box[:, 0] -= x0
    box[:, 1] -= y0
    box_i = np.round(box).astype(np.int32)
    cv2.fillPoly(rect_mask, [box_i], 255)
    rect_iou = _binary_iou(shape_mask, rect_mask)

    # Best fitted ellipse template. This naturally accepts circles and ellipses.
    ellipse_iou = 0.0
    if len(c) >= 5:
        try:
            (ecx, ecy), (ew, eh), eang = cv2.fitEllipse(c)
            if ew > 0 and eh > 0:
                ellipse_mask = np.zeros_like(shape_mask)
                center_local = (int(round(ecx - x0)), int(round(ecy - y0)))
                axes = (max(1, int(round(ew / 2.0))), max(1, int(round(eh / 2.0))))
                cv2.ellipse(
                    ellipse_mask, center_local, axes, float(eang),
                    0.0, 360.0, 255, thickness=cv2.FILLED, lineType=cv2.LINE_8
                )
                ellipse_iou = _binary_iou(shape_mask, ellipse_mask)
        except cv2.error:
            ellipse_iou = 0.0

    return {
        "regularity_solidity": float(desc.solidity),
        "regularity_convexity": float(convexity),
        "regularity_rect_iou": float(rect_iou),
        "regularity_ellipse_iou": float(ellipse_iou),
        "regularity_template_iou": float(max(rect_iou, ellipse_iou)),
    }


def measure_slot(
    descs: List[ShapeDesc], px_nm: float, image_name: str, condition: str
) -> Tuple[List[dict], List[ShapeDesc], Dict[int, ShapeDesc], List[dict]]:
    """
    Permissive ZEP200 slot measurement.

    Important design choice:
    - Filename prefix already tells us this is a slot image.
    - Therefore aspect ratio / near-circularity / orientation are NOT hard filters by default.
    - The base hard object filter is filled pixel area >= SLOT_MIN_AREA_PX (normally 30 px),
      plus the global border/oversize segmentation QC. Then an image-local relative-area QC
      rejects slot candidates below 0.1x or above 10x the mean candidate area.
    - We still calculate aspect and angle for diagnostics/output, but they do not reject
      the object unless the optional *_USE_* switches are enabled.
    """
    cand: List[ShapeDesc] = []
    rejected: List[dict] = []

    # Keep shape metrics for accepted objects so they can also be exported later.
    regularity_by_id: Dict[int, Dict[str, float]] = {}

    for d in descs:
        reasons = []
        codes = []
        pixel_area = contour_pixel_count(d.contour)
        angle_from_y = abs(90.0 - d.angle_from_x_deg)
        reg = slot_regularity_metrics(d)
        regularity_by_id[id(d)] = reg

        if pixel_area < SLOT_MIN_AREA_PX:
            codes.append("slot_area_too_small")
            reasons.append(f"filled pixel area {pixel_area} px < {SLOT_MIN_AREA_PX} px")
        if SLOT_USE_ASPECT_FILTER and d.aspect < SLOT_MIN_ASPECT:
            codes.append("slot_aspect_too_low")
            reasons.append(f"PCA aspect {d.aspect:.2f} < {SLOT_MIN_ASPECT:.2f}")
        if SLOT_USE_LENGTH_FILTER and d.length_px < SLOT_MIN_LENGTH_PX:
            codes.append("slot_too_short")
            reasons.append(f"PCA length {d.length_px:.1f}px < {SLOT_MIN_LENGTH_PX}px")
        if SLOT_USE_ANGLE_FILTER and SLOT_EXPECT_VERTICAL and angle_from_y > SLOT_MAX_ANGLE_FROM_Y_DEG:
            codes.append("slot_not_vertical")
            reasons.append(f"angle from Y {angle_from_y:.1f} deg > {SLOT_MAX_ANGLE_FROM_Y_DEG:.1f} deg")

        # Regularity is independent of aspect ratio: rectangle/square/circle/ellipse can all pass.
        if SLOT_USE_REGULARITY_FILTER:
            if reg["regularity_solidity"] < SLOT_MIN_REGULARITY_SOLIDITY:
                codes.append("slot_low_regularity_solidity")
                reasons.append(
                    f"solidity {reg['regularity_solidity']:.3f} < {SLOT_MIN_REGULARITY_SOLIDITY:.3f} "
                    "(large concavity/branching)"
                )
            if reg["regularity_convexity"] < SLOT_MIN_CONVEXITY:
                codes.append("slot_low_convexity")
                reasons.append(
                    f"convexity {reg['regularity_convexity']:.3f} < {SLOT_MIN_CONVEXITY:.3f} "
                    "(boundary too jagged/wavy)"
                )
            if reg["regularity_template_iou"] < SLOT_MIN_TEMPLATE_IOU:
                codes.append("slot_irregular_template_match")
                reasons.append(
                    f"best rectangle/ellipse IoU {reg['regularity_template_iou']:.3f} "
                    f"< {SLOT_MIN_TEMPLATE_IOU:.3f} "
                    f"(rect IoU={reg['regularity_rect_iou']:.3f}, "
                    f"ellipse IoU={reg['regularity_ellipse_iou']:.3f})"
                )

        if reasons:
            rejected.append(_make_reject_from_contour(
                d.contour, "slot_shape_filter", "+".join(codes), "; ".join(reasons), desc=d,
                pixel_area_px=pixel_area, aspect=d.aspect, length_px=d.length_px,
                angle_from_y_deg=angle_from_y,
                **reg,
            ))
            continue
        cand.append(d)

    # --------------------------------------------------------
    # Slot relative-area QC (image-local)
    # --------------------------------------------------------
    # Use the candidates that already passed the absolute 30-pixel area and
    # optional geometry checks to calculate this image's mean target area.
    # Then remove objects whose filled-pixel area is outside [0.1*mean, 10*mean].
    # This is intentionally slot-only and does not affect via/trench logic.
    slot_mean_area_px = np.nan
    slot_reference_median_area_px = np.nan
    if cand:
        base_areas = np.array([contour_pixel_count(d.contour) for d in cand], dtype=float)

        # A raw arithmetic mean can itself be pulled badly by a huge false-positive object.
        # To estimate the "average target area" robustly, first use the median only as a
        # broad seed to define the main target population (0.1x to 10x median), then take
        # the ARITHMETIC MEAN of that main population. The final accept/reject decision
        # below is still exactly based on 0.1x / 10x of this mean target area.
        slot_reference_median_area_px = float(np.median(base_areas))
        seed_low = slot_reference_median_area_px * SLOT_AREA_MIN_MEAN_FACTOR
        seed_high = slot_reference_median_area_px * SLOT_AREA_MAX_MEAN_FACTOR
        reference_pool = base_areas[(base_areas >= seed_low) & (base_areas <= seed_high)]
        if len(reference_pool) == 0:
            reference_pool = base_areas
        slot_mean_area_px = float(np.mean(reference_pool))

    if SLOT_USE_RELATIVE_AREA_FILTER and cand and np.isfinite(slot_mean_area_px) and slot_mean_area_px > 0:
        area_low_px = slot_mean_area_px * SLOT_AREA_MIN_MEAN_FACTOR
        area_high_px = slot_mean_area_px * SLOT_AREA_MAX_MEAN_FACTOR
        kept_after_area_qc: List[ShapeDesc] = []

        for d in cand:
            pixel_area = contour_pixel_count(d.contour)
            area_ratio = float(pixel_area / slot_mean_area_px)

            if pixel_area < area_low_px:
                rejected.append(_make_reject_from_contour(
                    d.contour,
                    "slot_relative_area_filter",
                    "slot_area_below_0.1x_mean",
                    (f"filled pixel area {pixel_area} px is {area_ratio:.4f}x image slot mean "
                     f"{slot_mean_area_px:.2f} px; below lower limit {area_low_px:.2f} px "
                     f"({SLOT_AREA_MIN_MEAN_FACTOR:.2f}x mean)"),
                    desc=d,
                    pixel_area_px=pixel_area,
                    slot_mean_area_px=slot_mean_area_px,
                    slot_reference_median_area_px=slot_reference_median_area_px,
                    area_ratio_to_mean=area_ratio,
                    lower_area_limit_px=area_low_px,
                    upper_area_limit_px=area_high_px,
                ))
                continue

            if pixel_area > area_high_px:
                rejected.append(_make_reject_from_contour(
                    d.contour,
                    "slot_relative_area_filter",
                    "slot_area_above_10x_mean",
                    (f"filled pixel area {pixel_area} px is {area_ratio:.4f}x image slot mean "
                     f"{slot_mean_area_px:.2f} px; above upper limit {area_high_px:.2f} px "
                     f"({SLOT_AREA_MAX_MEAN_FACTOR:.1f}x mean)"),
                    desc=d,
                    pixel_area_px=pixel_area,
                    slot_mean_area_px=slot_mean_area_px,
                    slot_reference_median_area_px=slot_reference_median_area_px,
                    area_ratio_to_mean=area_ratio,
                    lower_area_limit_px=area_low_px,
                    upper_area_limit_px=area_high_px,
                ))
                continue

            kept_after_area_qc.append(d)

        cand = kept_after_area_qc

    if len(cand) < SLOT_REQUIRE_AT_LEAST_N:
        for d in cand:
            rejected.append(_make_reject_from_contour(
                d.contour, "slot_group_filter", "slot_too_few_valid_objects",
                f"individually valid slot, but only {len(cand)} valid slot(s) remain after QC; need >= {SLOT_REQUIRE_AT_LEAST_N}",
                desc=d,
                slot_mean_area_px=slot_mean_area_px,
            ))
        return [], [], {}, rejected

    ys = np.array([d.center[1] for d in cand], dtype=np.float32)
    row_labels = cluster_1d_into_four_rows(ys)

    rows = []
    for i, (d, row_id) in enumerate(zip(cand, row_labels), 1):
        pts = d.contour.reshape(-1, 2).astype(np.float64)
        width_px_xy = float(pts[:, 0].max() - pts[:, 0].min())
        length_px_xy = float(pts[:, 1].max() - pts[:, 1].min())
        rows.append({
            "condition": condition,
            "image": image_name,
            "pattern": "slot",
            "object_id": i,
            "slot_row": int(row_id),
            "slot_width_nm": width_px_xy * px_nm,
            "slot_length_nm": length_px_xy * px_nm,
            "slot_aspect_ratio": length_px_xy / max(width_px_xy, 1e-12),
            "slot_angle_from_y_deg": abs(90.0 - d.angle_from_x_deg),
            "slot_pixel_area_px": contour_pixel_count(d.contour),
            "slot_mean_candidate_area_px": slot_mean_area_px,
            "slot_reference_median_area_px": slot_reference_median_area_px,
            "slot_area_ratio_to_mean": (
                contour_pixel_count(d.contour) / slot_mean_area_px
                if np.isfinite(slot_mean_area_px) and slot_mean_area_px > 0 else np.nan
            ),
            "slot_regularity_solidity": regularity_by_id.get(id(d), {}).get("regularity_solidity", np.nan),
            "slot_regularity_convexity": regularity_by_id.get(id(d), {}).get("regularity_convexity", np.nan),
            "slot_rect_iou": regularity_by_id.get(id(d), {}).get("regularity_rect_iou", np.nan),
            "slot_ellipse_iou": regularity_by_id.get(id(d), {}).get("regularity_ellipse_iou", np.nan),
            "slot_best_template_iou": regularity_by_id.get(id(d), {}).get("regularity_template_iou", np.nan),
            "pixel_size_nm": px_nm,
        })

    x_center = float(np.median([d.center[0] for d in cand]))
    reps: Dict[int, ShapeDesc] = {}
    for row_id in [1, 2, 3, 4]:
        idxs = [i for i, lab in enumerate(row_labels) if int(lab) == row_id]
        if idxs:
            reps[row_id] = min((cand[i] for i in idxs), key=lambda d: abs(d.center[0] - x_center))

    return rows, cand, reps, rejected



# ============================================================
# 7B) V12 ENHANCED SEARCH + ARRAY-ASSISTED VIA/SLOT RECOVERY
# ============================================================
def _mark_source(desc: ShapeDesc, source: str) -> ShapeDesc:
    setattr(desc, "_source", source)
    return desc


def _source_of(desc: ShapeDesc) -> str:
    return str(getattr(desc, "_source", "primary"))


def _bbox_from_desc(desc: ShapeDesc) -> Tuple[int, int, int, int]:
    return cv2.boundingRect(desc.contour.astype(np.int32))


def _bbox_iou_desc(a: ShapeDesc, b: ShapeDesc) -> float:
    ax, ay, aw, ah = _bbox_from_desc(a)
    bx, by, bw, bh = _bbox_from_desc(b)
    x1 = max(ax, bx)
    y1 = max(ay, by)
    x2 = min(ax + aw, bx + bw)
    y2 = min(ay + ah, by + bh)
    iw = max(0, x2 - x1)
    ih = max(0, y2 - y1)
    inter = iw * ih
    union = aw * ah + bw * bh - inter
    return float(inter / union) if union > 0 else 0.0


def _candidate_quality_v12(d: ShapeDesc, pattern: str) -> float:
    if pattern == "via":
        return float(1.5 * d.circularity + 1.0 * d.solidity - 0.15 * max(0.0, d.aspect - 1.0))
    if pattern == "slot":
        reg = slot_regularity_metrics(d)
        return float(
            0.8 * reg["regularity_solidity"]
            + 0.6 * reg["regularity_convexity"]
            + 1.2 * reg["regularity_template_iou"]
        )
    return float(d.solidity)


def _postprocess_binary_v12(mask: np.ndarray, close_iter: int = 1) -> np.ndarray:
    out = mask.astype(np.uint8)
    k = max(1, int(MORPH_KERNEL))
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k, k))
    if close_iter > 0:
        out = cv2.morphologyEx(out, cv2.MORPH_CLOSE, kernel, iterations=close_iter)
    return out


def enhanced_masks_v12(gray: np.ndarray) -> Dict[str, np.ndarray]:
    """V13 lightweight enhancement masks.

    IMPORTANT: this function is used both for the full image and for small grid-site
    ROIs. V12 generated many masks at several Gaussian scales and percentile levels.
    V13 deliberately keeps only three complementary masks to avoid candidate explosion:
      - CLAHE + Otsu
      - adaptive Gaussian threshold
      - one local-background subtraction scale
    """
    g = gray.astype(np.uint8)
    clahe = cv2.createCLAHE(
        clipLimit=float(ENHANCED_CLAHE_CLIP),
        tileGridSize=(int(ENHANCED_CLAHE_TILE), int(ENHANCED_CLAHE_TILE)),
    )
    ge = clahe.apply(g)
    masks: Dict[str, np.ndarray] = {}

    if FEATURE_POLARITY.lower() == "dark":
        _, m = cv2.threshold(ge, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
    else:
        _, m = cv2.threshold(ge, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    masks["enh_clahe_otsu"] = _postprocess_binary_v12(m, ENHANCED_MORPH_CLOSE_ITER)

    block = int(ENHANCED_ADAPTIVE_BLOCK)
    if block % 2 == 0:
        block += 1
    block = max(11, block)
    # Adaptive block cannot exceed ROI dimensions meaningfully. OpenCV accepts a
    # large block, but capping it improves small-ROI behavior and speed.
    min_dim = int(min(ge.shape[:2]))
    max_block = max(11, min_dim - 1 if min_dim % 2 == 0 else min_dim)
    if max_block % 2 == 0:
        max_block -= 1
    block = max(3, min(block, max_block))
    if block >= 3:
        thresh_type = cv2.THRESH_BINARY_INV if FEATURE_POLARITY.lower() == "dark" else cv2.THRESH_BINARY
        m = cv2.adaptiveThreshold(
            ge, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
            thresh_type, block, float(ENHANCED_ADAPTIVE_C)
        )
        masks["enh_adaptive"] = _postprocess_binary_v12(m, ENHANCED_MORPH_CLOSE_ITER)

    gf = ge.astype(np.float32)
    for sigma in ENHANCED_LOCAL_SIGMAS:
        # Skip absurdly large sigmas relative to a tiny local ROI.
        if min(ge.shape[:2]) < 10:
            break
        sig = min(float(sigma), max(2.0, 0.20 * min(ge.shape[:2])))
        bg = cv2.GaussianBlur(gf, (0, 0), sigmaX=sig, sigmaY=sig)
        response = (bg - gf) if FEATURE_POLARITY.lower() == "dark" else (gf - bg)
        response = np.maximum(response, 0.0)
        if float(response.max()) <= 1e-6:
            continue
        r8 = cv2.normalize(response, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
        _, m = cv2.threshold(r8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
        masks[f"enh_local_sigma{sig:g}"] = _postprocess_binary_v12(m, ENHANCED_MORPH_CLOSE_ITER)
        break

    return masks

def _fast_pattern_prefilter_v13(d: ShapeDesc, pattern: str) -> bool:
    """Very permissive prefilter used ONLY to cap fallback-search noise.

    Final via/slot QC remains unchanged later. This stage removes obvious threshold
    fragments before they can make merging or grid fitting expensive.
    """
    pa = contour_pixel_count(d.contour)
    if pa < MIN_COMPONENT_AREA_PX:
        return False
    if pattern == "via":
        return bool(
            d.solidity >= 0.50
            and d.circularity >= 0.18
            and d.aspect <= 2.4
            and min(d.length_px, d.width_px) >= 4.0
        )
    if pattern == "slot":
        # Slot may be rectangular, square, circular or elliptical after trimming.
        # Keep this intentionally loose; regularity/area QC is applied later.
        return bool(d.solidity >= 0.48 and d.length_px >= 4.0 and d.width_px >= 3.0)
    return True


def extract_loose_descriptors_v12(
    mask: np.ndarray,
    source: str,
    pattern: str = "unknown",
    max_keep: int = MAX_ENHANCED_CANDIDATES_PER_MASK,
) -> List[ShapeDesc]:
    """Fast candidate extraction for fallback search, with hard candidate cap."""
    h, w = mask.shape
    max_area = MAX_COMPONENT_AREA_FRAC * h * w
    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    out: List[ShapeDesc] = []
    for c in contours:
        # contourArea is much cheaper than drawing every tiny contour first.
        area = float(cv2.contourArea(c))
        if area <= 0 or area > max_area:
            continue
        x, y, cw, ch = cv2.boundingRect(c)
        if cw * ch < MIN_COMPONENT_AREA_PX:
            continue
        if (
            x <= BORDER_MARGIN_PX or y <= BORDER_MARGIN_PX
            or x + cw >= w - BORDER_MARGIN_PX
            or y + ch >= h - BORDER_MARGIN_PX
        ):
            continue
        d = describe_contour(c)
        if d is None:
            continue
        if not _fast_pattern_prefilter_v13(d, pattern):
            continue
        out.append(_mark_source(d, source))

    # Prefer larger / more regular candidates if a noisy threshold produces hundreds.
    if len(out) > max_keep:
        out.sort(
            key=lambda d: (
                _candidate_quality_v12(d, pattern),
                math.log1p(max(1, contour_pixel_count(d.contour))),
            ),
            reverse=True,
        )
        out = out[:max_keep]
    return out

def merge_candidate_descs_v12(
    primary: List[ShapeDesc], extra: List[ShapeDesc], pattern: str
) -> Tuple[List[ShapeDesc], int]:
    """V13 spatial-hash duplicate merge.

    V12 compared every new candidate with every existing candidate, which can become
    O(N^2) when a blurry adaptive threshold creates many fragments. V13 only compares
    candidates in the same/neighboring spatial buckets.
    """
    merged: List[ShapeDesc] = []
    cell = max(8.0, float(MERGE_SPATIAL_CELL_PX))
    buckets: Dict[Tuple[int, int], List[int]] = {}

    def key_for(d: ShapeDesc) -> Tuple[int, int]:
        return (int(math.floor(float(d.center[0]) / cell)), int(math.floor(float(d.center[1]) / cell)))

    def add_to_bucket(idx: int, d: ShapeDesc) -> None:
        buckets.setdefault(key_for(d), []).append(idx)

    for d0 in primary:
        d = _mark_source(d0, "primary")
        idx = len(merged)
        merged.append(d)
        add_to_bucket(idx, d)
        if len(merged) >= MAX_CANDIDATE_POOL:
            break

    added = 0
    for d in extra:
        if len(merged) >= MAX_CANDIDATE_POOL:
            break
        kx, ky = key_for(d)
        neighbor_indices: List[int] = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                neighbor_indices.extend(buckets.get((kx + dx, ky + dy), []))

        best_idx = -1
        for i in neighbor_indices:
            e = merged[i]
            iou = _bbox_iou_desc(d, e)
            dist = float(np.linalg.norm(d.center - e.center))
            sx = max(2.0, min(d.length_px, d.width_px, e.length_px, e.width_px))
            if iou >= 0.38 or dist <= max(3.0, 0.45 * sx):
                best_idx = i
                break

        if best_idx < 0:
            idx = len(merged)
            merged.append(d)
            add_to_bucket(idx, d)
            added += 1
        else:
            old = merged[best_idx]
            if _candidate_quality_v12(d, pattern) > _candidate_quality_v12(old, pattern) + 0.05:
                # Keep bucket entry of old position too; the replacement should be close
                # because it was classified as a duplicate. Also add its new bucket.
                merged[best_idx] = d
                add_to_bucket(best_idx, d)

    return merged, added

def build_candidate_pool_v12(
    gray: np.ndarray,
    primary_descs: List[ShapeDesc],
    pattern: str,
    image_name: str = "",
) -> Tuple[List[ShapeDesc], Dict[str, int]]:
    """V13 fast full-image search.

    Start from primary Otsu. Only a few fallback masks are evaluated, one at a time,
    with per-mask and total candidate caps. Heavy enhancement is reserved for missing
    lattice sites later.
    """
    for d in primary_descs:
        _mark_source(d, "primary")
    pool = list(primary_descs[:MAX_CANDIDATE_POOL])
    meta = {
        "primary_candidates": len(primary_descs),
        "enhanced_candidates_raw": 0,
        "enhanced_candidates_added": 0,
        "enhanced_masks_used": 0,
        "candidate_pool_capped": 0,
    }
    if not USE_ENHANCED_SEARCH_V12 or pattern not in {"via", "slot"}:
        return pool, meta

    masks = enhanced_masks_v12(gray)
    # Always try CLAHE/Otsu. Only use more permissive masks if the current pool is
    # still too sparse to fit a grid reliably.
    preferred = ["enh_clahe_otsu", "enh_adaptive"]
    preferred += [k for k in masks.keys() if k.startswith("enh_local_")]

    for name in preferred:
        if name not in masks:
            continue
        if len(pool) >= MAX_CANDIDATE_POOL:
            meta["candidate_pool_capped"] = 1
            break
        if name != "enh_clahe_otsu" and len(pool) >= FULL_IMAGE_FALLBACK_MIN_PRIMARY:
            # Once enough global candidates exist, stop broad full-image searching.
            break

        ds = extract_loose_descriptors_v12(
            masks[name], source=name, pattern=pattern,
            max_keep=MAX_ENHANCED_CANDIDATES_PER_MASK,
        )
        meta["enhanced_masks_used"] += 1
        meta["enhanced_candidates_raw"] += len(ds)
        pool, n_added = merge_candidate_descs_v12(pool, ds, pattern)
        meta["enhanced_candidates_added"] += n_added
        if image_name:
            print(
                f"      fallback mask {name}: raw={len(ds)}, added={n_added}, pool={len(pool)}",
                flush=True,
            )

    return pool, meta

def _cluster_axis_values_v12(values: List[float], tol: float) -> Tuple[np.ndarray, np.ndarray]:
    if not values:
        return np.array([], dtype=float), np.array([], dtype=int)
    vals = sorted(float(v) for v in values)
    groups: List[List[float]] = [[vals[0]]]
    for v in vals[1:]:
        med = float(np.median(groups[-1]))
        if abs(v - med) <= tol:
            groups[-1].append(v)
        else:
            groups.append([v])
    centers = np.array([float(np.median(g)) for g in groups], dtype=float)
    counts = np.array([len(g) for g in groups], dtype=int)
    return centers, counts


def _estimate_regular_axis_v12(
    obs: np.ndarray,
    weights: np.ndarray,
    min_pitch: float,
    image_limit: float,
) -> Tuple[np.ndarray, float, float]:
    """
    Fit an arithmetic progression to observed row/column centers.
    Returns (grid_positions, pitch, occupancy). Missing *internal* sites are filled.
    """
    obs = np.asarray(obs, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if len(obs) == 0:
        return np.array([], dtype=float), np.nan, 0.0
    order = np.argsort(obs)
    obs = obs[order]
    weights = weights[order]
    if len(obs) == 1:
        return obs.copy(), np.nan, 1.0
    if len(obs) == 2:
        pitch = float(obs[1] - obs[0])
        if pitch < min_pitch:
            return obs.copy(), pitch, 1.0
        return obs.copy(), pitch, 1.0

    # Generate plausible fundamental pitches from pair differences / integer multiples.
    cands: List[float] = []
    for i in range(len(obs)):
        for j in range(i + 1, len(obs)):
            diff = float(obs[j] - obs[i])
            if diff <= 0:
                continue
            for k in range(1, 7):
                p = diff / k
                if p >= min_pitch and p <= image_limit:
                    cands.append(p)
    if not cands:
        diffs = np.diff(obs)
        pitch = float(np.median(diffs))
        return obs.copy(), pitch, 1.0

    # Avoid evaluating thousands of nearly identical values.
    cands = sorted(set(round(p * 2.0) / 2.0 for p in cands if p > 0))
    best = None
    total_w = float(np.sum(weights)) if float(np.sum(weights)) > 0 else float(len(obs))

    for pitch in cands:
        tol = max(2.0, ARRAY_GRID_RESIDUAL_FRAC * pitch)
        for phase0 in obs:
            n = np.rint((obs - phase0) / pitch).astype(int)
            phase = float(np.average(obs - n * pitch, weights=weights))
            pred = phase + n * pitch
            resid = np.abs(obs - pred)
            good = resid <= tol
            if int(np.count_nonzero(good)) < 2:
                continue
            ng = n[good]
            unique_n = np.unique(ng)
            nmin, nmax = int(unique_n.min()), int(unique_n.max())
            grid_count = max(1, nmax - nmin + 1)
            occupancy = float(len(unique_n) / grid_count)
            weighted_support = float(np.sum(weights[good]) / total_w)
            rms = float(np.sqrt(np.mean(resid[good] ** 2)))
            # True pitch tends to maximize occupancy; half-pitch solutions have many empty sites.
            score = 2.0 * weighted_support + 1.8 * occupancy - 0.03 * grid_count - 0.02 * rms
            if best is None or score > best[0]:
                best = (score, phase, pitch, nmin, nmax, occupancy)

    if best is None:
        pitch = float(np.median(np.diff(obs)))
        return obs.copy(), pitch, 1.0

    _, phase, pitch, nmin, nmax, occupancy = best
    positions = phase + np.arange(nmin, nmax + 1, dtype=float) * pitch
    positions = positions[(positions >= 0) & (positions <= image_limit - 1)]
    return positions.astype(float), float(pitch), float(occupancy)


def _trim_axis_to_count_v12(
    positions: np.ndarray,
    count: int,
    observations: List[float],
) -> np.ndarray:
    positions = np.asarray(positions, dtype=float)
    if len(positions) <= count:
        return positions
    obs = np.asarray(observations, dtype=float)
    best = None
    for start in range(0, len(positions) - count + 1):
        win = positions[start:start + count]
        if obs.size:
            dist = np.min(np.abs(obs[:, None] - win[None, :]), axis=1)
            score = -float(np.median(dist)) - 0.1 * float(np.mean(dist))
        else:
            score = 0.0
        if best is None or score > best[0]:
            best = (score, win)
    return best[1].copy() if best is not None else positions[:count]



def _extend_axis_to_count_v12(
    positions: np.ndarray,
    pitch: float,
    count: int,
    image_limit: float,
    observations: List[float],
) -> np.ndarray:
    """Extend a fitted consecutive lattice to a known count (slot Y has exactly 4 rows)."""
    positions = np.asarray(positions, dtype=float)
    if len(positions) >= count or not np.isfinite(pitch) or pitch <= 0 or len(positions) == 0:
        return positions
    missing = count - len(positions)
    obs = np.asarray(observations, dtype=float)
    candidates = []
    for n_before in range(missing + 1):
        n_after = missing - n_before
        start = positions[0] - n_before * pitch
        cand = start + np.arange(count, dtype=float) * pitch
        if cand[0] < 0 or cand[-1] > image_limit - 1:
            continue
        if obs.size:
            dist = np.min(np.abs(obs[:, None] - cand[None, :]), axis=1)
            # Favor candidates explaining all observations, while not letting one
            # far-away outlier dominate the decision.
            score = float(np.median(dist) + 0.20 * np.mean(np.minimum(dist, pitch)))
        else:
            score = 0.0
        candidates.append((score, cand))
    if not candidates:
        return positions
    candidates.sort(key=lambda x: x[0])
    return candidates[0][1]


def _grid_match_tolerance_v12(pitch: float, size_ref: float) -> float:
    vals = [3.0, ARRAY_MATCH_SIZE_FACTOR * max(1.0, size_ref)]
    if np.isfinite(pitch) and pitch > 0:
        vals.append(ARRAY_MATCH_PITCH_FRAC * pitch)
    return float(max(vals))


def _via_strict_metrics_v12(d: ShapeDesc) -> Tuple[bool, List[str], List[str]]:
    codes: List[str] = []
    reasons: List[str] = []
    if d.circularity < VIA_MIN_CIRCULARITY:
        codes.append("via_low_circularity")
        reasons.append(f"circularity {d.circularity:.3f} < {VIA_MIN_CIRCULARITY:.3f}")
    if d.aspect > VIA_MAX_AXIS_RATIO:
        codes.append("via_axis_ratio_too_high")
        reasons.append(f"axis ratio {d.aspect:.3f} > {VIA_MAX_AXIS_RATIO:.3f}")
    if d.solidity < VIA_MIN_SOLIDITY:
        codes.append("via_low_solidity")
        reasons.append(f"solidity {d.solidity:.3f} < {VIA_MIN_SOLIDITY:.3f}")
    if min(d.length_px, d.width_px) < VIA_MIN_DIAMETER_PX:
        codes.append("via_too_small")
        reasons.append(f"minimum diameter {min(d.length_px, d.width_px):.1f}px < {VIA_MIN_DIAMETER_PX}px")
    return len(codes) == 0, codes, reasons


def _via_grid_relaxed_ok_v12(d: ShapeDesc) -> bool:
    return bool(
        contour_pixel_count(d.contour) >= MIN_COMPONENT_AREA_PX
        and d.circularity >= VIA_GRID_MIN_CIRCULARITY
        and d.aspect <= VIA_GRID_MAX_AXIS_RATIO
        and d.solidity >= VIA_GRID_MIN_SOLIDITY
        and min(d.length_px, d.width_px) >= max(5.0, 0.65 * VIA_MIN_DIAMETER_PX)
    )


def _slot_regularity_ok_v12(d: ShapeDesc, relaxed: bool = False) -> Tuple[bool, Dict[str, float], List[str], List[str]]:
    reg = slot_regularity_metrics(d)
    codes: List[str] = []
    reasons: List[str] = []
    if relaxed:
        min_sol = SLOT_GRID_MIN_SOLIDITY
        min_conv = SLOT_GRID_MIN_CONVEXITY
        min_iou = SLOT_GRID_MIN_TEMPLATE_IOU
    else:
        min_sol = SLOT_MIN_REGULARITY_SOLIDITY
        min_conv = SLOT_MIN_CONVEXITY
        min_iou = SLOT_MIN_TEMPLATE_IOU
    if reg["regularity_solidity"] < min_sol:
        codes.append("slot_low_regularity_solidity")
        reasons.append(f"solidity {reg['regularity_solidity']:.3f} < {min_sol:.3f}")
    if reg["regularity_convexity"] < min_conv:
        codes.append("slot_low_convexity")
        reasons.append(f"convexity {reg['regularity_convexity']:.3f} < {min_conv:.3f}")
    if reg["regularity_template_iou"] < min_iou:
        codes.append("slot_irregular_template_match")
        reasons.append(
            f"best rectangle/ellipse IoU {reg['regularity_template_iou']:.3f} < {min_iou:.3f}"
        )
    return len(codes) == 0, reg, codes, reasons


def _local_candidate_score_v12(
    d: ShapeDesc,
    expected_xy: Tuple[float, float],
    expected_area: float,
    pattern: str,
) -> float:
    dist = float(np.linalg.norm(d.center - np.asarray(expected_xy, dtype=float)))
    size_scale = max(4.0, math.sqrt(max(1.0, expected_area)))
    area = float(contour_pixel_count(d.contour))
    area_term = abs(math.log(max(area, 1.0) / max(expected_area, 1.0)))
    q = _candidate_quality_v12(d, pattern)
    return 1.8 * dist / size_scale + 0.55 * area_term - 0.35 * q


def recover_local_object_v12(
    gray: np.ndarray,
    expected_xy: Tuple[float, float],
    expected_w: float,
    expected_h: float,
    pattern: str,
    x_pitch: float,
    y_pitch: float,
    expected_area: float,
) -> Optional[ShapeDesc]:
    """Search only near one missing lattice site using all enhanced threshold modes."""
    H, W = gray.shape
    cx, cy = expected_xy
    half_w = ARRAY_LOCAL_ROI_SIZE_FACTOR * max(4.0, expected_w) / 2.0
    half_h = ARRAY_LOCAL_ROI_SIZE_FACTOR * max(4.0, expected_h) / 2.0
    if np.isfinite(x_pitch) and x_pitch > 0:
        half_w = min(max(half_w, expected_w), ARRAY_LOCAL_ROI_PITCH_FRAC * x_pitch)
    if np.isfinite(y_pitch) and y_pitch > 0:
        half_h = min(max(half_h, expected_h), ARRAY_LOCAL_ROI_PITCH_FRAC * y_pitch)
        if pattern == "slot":
            # The longest of the four slot rows can nearly fill the vertical pitch.
            # Give local recovery enough Y context while still avoiding the neighbor row.
            half_h = max(half_h, 0.44 * y_pitch)
            half_h = min(half_h, 0.49 * y_pitch)
    half_w = max(10.0, half_w)
    half_h = max(10.0, half_h)

    x0 = max(0, int(math.floor(cx - half_w)))
    x1 = min(W, int(math.ceil(cx + half_w + 1)))
    y0 = max(0, int(math.floor(cy - half_h)))
    y1 = min(H, int(math.ceil(cy + half_h + 1)))
    if x1 - x0 < 8 or y1 - y0 < 8:
        return None
    roi = gray[y0:y1, x0:x1]

    masks = {"local_primary": segment_features(roi)}
    masks.update(enhanced_masks_v12(roi))

    # Additional site-guided threshold: compare the expected-center patch with the
    # ROI border background. This is especially useful for a globally blurry image
    # where the slot is still visible by eye but histogram thresholding finds nothing.
    rh0, rw0 = roi.shape
    band = max(2, min(6, min(rh0, rw0) // 8))
    border_vals = np.concatenate([
        roi[:band, :].ravel(), roi[-band:, :].ravel(),
        roi[:, :band].ravel(), roi[:, -band:].ravel(),
    ]).astype(np.float32)
    cxr = int(round(cx - x0))
    cyr = int(round(cy - y0))
    phw = max(3, int(round(min(expected_w, rw0 / 2) * 0.55)))
    phh = max(3, int(round(min(expected_h, rh0 / 2) * 0.40)))
    px0, px1 = max(0, cxr - phw), min(rw0, cxr + phw + 1)
    py0, py1 = max(0, cyr - phh), min(rh0, cyr + phh + 1)
    core = roi[py0:py1, px0:px1].astype(np.float32)
    if border_vals.size and core.size:
        bg_med = float(np.median(border_vals))
        if FEATURE_POLARITY.lower() == "dark":
            target = float(np.percentile(core, 28.0))
            if bg_med - target >= 1.5:
                thr = 0.52 * bg_med + 0.48 * target
                msite = (roi.astype(np.float32) < thr).astype(np.uint8) * 255
                masks["local_site_guided"] = _postprocess_binary_v12(msite, 1)
        else:
            target = float(np.percentile(core, 72.0))
            if target - bg_med >= 1.5:
                thr = 0.52 * bg_med + 0.48 * target
                msite = (roi.astype(np.float32) > thr).astype(np.uint8) * 255
                masks["local_site_guided"] = _postprocess_binary_v12(msite, 1)

    candidates: List[ShapeDesc] = []
    for source, mask in masks.items():
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        for c in contours:
            if contour_pixel_count(c) < MIN_COMPONENT_AREA_PX:
                continue
            bx, by, bw, bh = cv2.boundingRect(c)
            if bx <= 1 or by <= 1 or bx + bw >= roi.shape[1] - 1 or by + bh >= roi.shape[0] - 1:
                continue
            cg = c.astype(np.int32).copy()
            cg[:, 0, 0] += x0
            cg[:, 0, 1] += y0
            d = describe_contour(cg)
            if d is None:
                continue
            _mark_source(d, "grid_recovered")
            center_tol = ARRAY_LOCAL_CENTER_TOL_FRAC * max(half_w, half_h)
            if float(np.linalg.norm(d.center - np.asarray(expected_xy))) > center_tol:
                continue
            if pattern == "via":
                if not _via_grid_relaxed_ok_v12(d):
                    continue
            else:
                ok, _, _, _ = _slot_regularity_ok_v12(d, relaxed=True)
                if not ok:
                    continue
            candidates.append(d)

    if not candidates:
        return None
    return min(candidates, key=lambda d: _local_candidate_score_v12(d, expected_xy, expected_area, pattern))


def infer_via_grid_v12(seeds: List[ShapeDesc], image_shape: Tuple[int, int]) -> Dict[str, object]:
    H, W = image_shape
    meta: Dict[str, object] = {
        "used": False, "x_positions": np.array([]), "y_positions": np.array([]),
        "x_pitch": np.nan, "y_pitch": np.nan, "occupancy": 0.0,
    }
    if len(seeds) < ARRAY_MIN_GRID_SEEDS:
        return meta
    med_w = float(np.median([max(1.0, _contour_xy_bounds(d)[1] - _contour_xy_bounds(d)[0]) for d in seeds]))
    med_h = float(np.median([max(1.0, _contour_xy_bounds(d)[3] - _contour_xy_bounds(d)[2]) for d in seeds]))
    xtol = max(3.0, min(ARRAY_CENTER_CLUSTER_TOL_PX, 0.65 * med_w))
    ytol = max(3.0, min(ARRAY_CENTER_CLUSTER_TOL_PX, 0.65 * med_h))
    xc, xcount = _cluster_axis_values_v12([d.center[0] for d in seeds], xtol)
    yc, ycount = _cluster_axis_values_v12([d.center[1] for d in seeds], ytol)
    # Repeated-array rows/columns normally contain multiple objects. Single-support
    # coordinate clusters are often false positives from secondary thresholding.
    xmask = xcount >= ARRAY_MIN_AXIS_CLUSTER_SUPPORT
    ymask = ycount >= ARRAY_MIN_AXIS_CLUSTER_SUPPORT
    if int(np.count_nonzero(xmask)) >= ARRAY_MIN_AXIS_POSITIONS:
        xc, xcount = xc[xmask], xcount[xmask]
    if int(np.count_nonzero(ymask)) >= ARRAY_MIN_AXIS_POSITIONS:
        yc, ycount = yc[ymask], ycount[ymask]
    if len(xc) < ARRAY_MIN_AXIS_POSITIONS or len(yc) < ARRAY_MIN_AXIS_POSITIONS:
        return meta

    xp, x_pitch, x_occ = _estimate_regular_axis_v12(xc, xcount, max(2.0 * xtol, 1.2 * med_w), W)
    yp, y_pitch, y_occ = _estimate_regular_axis_v12(yc, ycount, max(2.0 * ytol, 1.2 * med_h), H)
    if len(xp) < 2 or len(yp) < 2:
        return meta
    expected = len(xp) * len(yp)
    if expected > ARRAY_MAX_EXPECTED_OBJECTS:
        return meta
    # Occupancy is based on DISTINCT lattice sites, not raw candidate count;
    # multi-threshold duplicates therefore cannot inflate confidence above 100%.
    tolx_site = max(xtol, 0.30 * x_pitch if np.isfinite(x_pitch) else xtol)
    toly_site = max(ytol, 0.30 * y_pitch if np.isfinite(y_pitch) else ytol)
    supported_sites = set()
    for d in seeds:
        ix, dx = _nearest_expected_index_v12(float(d.center[0]), xp)
        iy, dy = _nearest_expected_index_v12(float(d.center[1]), yp)
        if ix >= 0 and iy >= 0 and dx <= tolx_site and dy <= toly_site:
            supported_sites.add((ix, iy))
    occ = len(supported_sites) / max(1, expected)
    if occ < ARRAY_GRID_MIN_OCCUPANCY:
        return meta
    meta.update({
        "used": True, "x_positions": xp, "y_positions": yp,
        "x_pitch": x_pitch, "y_pitch": y_pitch,
        "occupancy": occ, "supported_sites": len(supported_sites),
        "median_w": med_w, "median_h": med_h,
    })
    return meta


def infer_slot_grid_v12(seeds: List[ShapeDesc], image_shape: Tuple[int, int]) -> Dict[str, object]:
    H, W = image_shape
    meta: Dict[str, object] = {
        "used": False, "x_positions": np.array([]), "y_positions": np.array([]),
        "x_pitch": np.nan, "y_pitch": np.nan, "occupancy": 0.0,
    }
    if len(seeds) < ARRAY_MIN_GRID_SEEDS:
        return meta

    med_w = float(np.median([max(1.0, _contour_xy_bounds(d)[1] - _contour_xy_bounds(d)[0]) for d in seeds]))
    med_h = float(np.median([max(1.0, _contour_xy_bounds(d)[3] - _contour_xy_bounds(d)[2]) for d in seeds]))
    # Centers of repeated objects in the same row/column should align much more tightly
    # than their physical length/width, so use a fixed small tolerance rather than height.
    xtol = max(4.0, ARRAY_CENTER_CLUSTER_TOL_PX)
    ytol = max(4.0, ARRAY_CENTER_CLUSTER_TOL_PX)
    xc, xcount = _cluster_axis_values_v12([d.center[0] for d in seeds], xtol)
    yc, ycount = _cluster_axis_values_v12([d.center[1] for d in seeds], ytol)
    xmask = xcount >= ARRAY_MIN_AXIS_CLUSTER_SUPPORT
    ymask = ycount >= ARRAY_MIN_AXIS_CLUSTER_SUPPORT
    if int(np.count_nonzero(xmask)) >= 2:
        xc, xcount = xc[xmask], xcount[xmask]
    if int(np.count_nonzero(ymask)) >= 2:
        yc, ycount = yc[ymask], ycount[ymask]
    if len(xc) < 2 or len(yc) < 2:
        return meta

    xp, x_pitch, x_occ = _estimate_regular_axis_v12(xc, xcount, max(2.0 * xtol, 1.5 * med_w), W)
    yp, y_pitch, y_occ = _estimate_regular_axis_v12(yc, ycount, max(2.0 * ytol, 0.35 * med_h), H)
    if len(xp) < 2 or len(yp) < 2:
        return meta

    # This design has four Y levels. If extra spurious Y levels exist, retain the
    # consecutive four-level window best supported by observations.
    yp = _trim_axis_to_count_v12(yp, 4, [d.center[1] for d in seeds])
    if len(yp) < 4:
        yp = _extend_axis_to_count_v12(
            yp, y_pitch, 4, H, [d.center[1] for d in seeds]
        )
    if len(yp) != 4:
        return meta
    expected = len(xp) * 4
    if expected > ARRAY_MAX_EXPECTED_OBJECTS:
        return meta
    tolx_site = max(xtol, 0.30 * x_pitch if np.isfinite(x_pitch) else xtol)
    toly_site = max(ytol, 0.30 * y_pitch if np.isfinite(y_pitch) else ytol)
    supported_sites = set()
    for d in seeds:
        ix, dx = _nearest_expected_index_v12(float(d.center[0]), xp)
        iy, dy = _nearest_expected_index_v12(float(d.center[1]), yp)
        if ix >= 0 and iy >= 0 and dx <= tolx_site and dy <= toly_site:
            supported_sites.add((ix, iy))
    occ = len(supported_sites) / max(1, expected)
    if occ < ARRAY_GRID_MIN_OCCUPANCY:
        return meta
    meta.update({
        "used": True, "x_positions": xp, "y_positions": yp,
        "x_pitch": x_pitch, "y_pitch": y_pitch,
        "occupancy": occ, "supported_sites": len(supported_sites),
        "median_w": med_w, "median_h": med_h,
    })
    return meta


def _nearest_expected_index_v12(value: float, positions: np.ndarray) -> Tuple[int, float]:
    if len(positions) == 0:
        return -1, np.inf
    idx = int(np.argmin(np.abs(positions - value)))
    return idx, float(abs(positions[idx] - value))


def _select_best_near_site_v12(
    pool: List[ShapeDesc],
    used_ids: set,
    expected_xy: Tuple[float, float],
    tol_x: float,
    tol_y: float,
    pattern: str,
    expected_area: float,
    slot_area_limits: Optional[Tuple[float, float]] = None,
) -> Optional[ShapeDesc]:
    candidates = []
    ex, ey = expected_xy
    for d in pool:
        if id(d) in used_ids:
            continue
        if abs(float(d.center[0]) - ex) > tol_x or abs(float(d.center[1]) - ey) > tol_y:
            continue
        pa = contour_pixel_count(d.contour)
        if pattern == "via":
            if not _via_grid_relaxed_ok_v12(d):
                continue
        else:
            ok, _, _, _ = _slot_regularity_ok_v12(d, relaxed=True)
            if not ok:
                continue
            if pa < SLOT_MIN_AREA_PX:
                continue
            if slot_area_limits is not None:
                lo, hi = slot_area_limits
                if pa < lo or pa > hi:
                    continue
        score = _local_candidate_score_v12(d, expected_xy, expected_area, pattern)
        # Primary/strict detections win ties over enhanced detections.
        if _source_of(d) == "primary":
            score -= 0.10
        candidates.append((score, d))
    if not candidates:
        return None
    candidates.sort(key=lambda t: t[0])
    return candidates[0][1]


def _desc_close_to_any_selected_v12(d: ShapeDesc, selected: List[ShapeDesc]) -> bool:
    for s in selected:
        if _bbox_iou_desc(d, s) >= 0.35:
            return True
        dist = float(np.linalg.norm(d.center - s.center))
        ref = max(3.0, 0.5 * min(d.length_px, d.width_px, s.length_px, s.width_px))
        if dist <= ref:
            return True
    return False



def _grid_model_score_v12(grid: Dict[str, object], n_seeds: int) -> float:
    if not grid.get("used", False):
        return -1e9
    expected = int(len(grid.get("x_positions", [])) * len(grid.get("y_positions", [])))
    occ = float(grid.get("occupancy", 0.0))
    supported = int(grid.get("supported_sites", round(occ * expected)))
    # Reward distinct-site occupancy and broad repeated-array support. Raw duplicate
    # candidate count is intentionally not rewarded.
    return occ + 0.010 * min(expected, 60) + 0.008 * min(supported, 60)



# ============================================================
# V14: GENERAL 2D LATTICE (RECTANGULAR / ROTATED / RHOMBIC)
# ============================================================
def _line_angle_diff_rad_v14(a: np.ndarray, b: np.ndarray) -> float:
    """Undirected angle between two vectors, in [0, pi/2]."""
    aa = math.atan2(float(a[1]), float(a[0]))
    bb = math.atan2(float(b[1]), float(b[0]))
    d = abs(aa - bb) % math.pi
    return min(d, math.pi - d)


def _canonical_line_vector_v14(v: np.ndarray) -> np.ndarray:
    """Give an undirected displacement a stable sign without assuming x/y axes."""
    v = np.asarray(v, dtype=float).copy()
    if v[0] < 0 or (abs(v[0]) < 1e-12 and v[1] < 0):
        v *= -1.0
    return v


def _lattice_vector_candidates_v14(seeds: List[ShapeDesc]) -> List[np.ndarray]:
    """
    Build a small set of repeated nearest-neighbour displacement vectors.
    Works for horizontal/vertical, rotated rectangular and rhombic lattices.
    """
    if len(seeds) < 2:
        return []
    P = np.asarray([d.center for d in seeds], dtype=float)
    min_dim = float(np.median([
        max(1.0, min(_contour_xy_bounds(d)[1] - _contour_xy_bounds(d)[0],
                     _contour_xy_bounds(d)[3] - _contour_xy_bounds(d)[2]))
        for d in seeds
    ]))
    min_len = max(5.0, 0.65 * min_dim)

    raw: List[np.ndarray] = []
    k = min(LATTICE_NEAREST_NEIGHBORS, len(P) - 1)
    for i in range(len(P)):
        dif = P - P[i]
        dist = np.linalg.norm(dif, axis=1)
        order = np.argsort(dist)
        used = 0
        for j in order:
            if j == i or dist[j] < min_len:
                continue
            raw.append(_canonical_line_vector_v14(dif[j]))
            used += 1
            if used >= k:
                break

    # Greedy cluster by undirected angle + relative length. Vectors in a cluster
    # are sign-aligned before averaging, so ~1 deg and ~179 deg line directions merge.
    clusters: List[List[np.ndarray]] = []
    ang_tol = math.radians(LATTICE_VECTOR_ANGLE_TOL_DEG)
    for v0 in sorted(raw, key=lambda q: float(np.linalg.norm(q))):
        lv = float(np.linalg.norm(v0))
        if lv <= 0:
            continue
        placed = False
        for g in clusters:
            ref = np.mean(np.vstack(g), axis=0)
            lr = float(np.linalg.norm(ref))
            if lr <= 0:
                continue
            if _line_angle_diff_rad_v14(v0, ref) <= ang_tol and abs(lv - lr) / max(lv, lr) <= LATTICE_VECTOR_LENGTH_REL_TOL:
                vv = v0.copy()
                if float(np.dot(vv, ref)) < 0:
                    vv *= -1.0
                g.append(vv)
                placed = True
                break
        if not placed:
            clusters.append([v0.copy()])

    scored: List[Tuple[float, np.ndarray]] = []
    for g in clusters:
        G = np.vstack(g)
        ref = np.mean(G, axis=0)
        lens = np.linalg.norm(G, axis=1)
        med_len = float(np.median(lens))
        if np.linalg.norm(ref) < 1e-9:
            continue
        ref = ref / np.linalg.norm(ref) * med_len
        support = len(g)
        # Favor repeated short vectors; long multi-pitch vectors can still survive
        # but are ranked behind a true nearest-neighbour basis.
        score = float(support) / max(1.0, math.sqrt(med_len))
        scored.append((score, ref))

    scored.sort(key=lambda z: z[0], reverse=True)
    out: List[np.ndarray] = []
    for _, v in scored:
        if all(
            _line_angle_diff_rad_v14(v, q) > math.radians(4.0)
            or abs(np.linalg.norm(v) - np.linalg.norm(q)) / max(np.linalg.norm(v), np.linalg.norm(q)) > 0.10
            for q in out
        ):
            out.append(v)
        if len(out) >= LATTICE_MAX_VECTOR_CANDIDATES:
            break
    return out


def _fit_affine_lattice_basis_v14(
    centers: np.ndarray,
    a0: np.ndarray,
    b0: np.ndarray,
    max_expected: int,
) -> Optional[Dict[str, object]]:
    """Fit p = origin + i*a + j*b using integer lattice coordinates."""
    P = np.asarray(centers, dtype=float)
    a0 = np.asarray(a0, dtype=float)
    b0 = np.asarray(b0, dtype=float)
    la, lb = float(np.linalg.norm(a0)), float(np.linalg.norm(b0))
    if la <= 1 or lb <= 1:
        return None
    angle = math.degrees(math.acos(np.clip(float(np.dot(a0, b0) / (la * lb)), -1.0, 1.0)))
    if angle < LATTICE_MIN_BASIS_ANGLE_DEG or angle > LATTICE_MAX_BASIS_ANGLE_DEG:
        return None

    best = None
    origin_candidates = P[:min(len(P), 8)]
    for p0 in origin_candidates:
        origin = p0.copy()
        a = a0.copy()
        b = b0.copy()
        ok = True
        for _ in range(max(1, LATTICE_REFINE_ITERS)):
            B = np.column_stack([a, b])
            det = float(np.linalg.det(B))
            if abs(det) < 1e-5:
                ok = False
                break
            invB = np.linalg.inv(B)
            coeff = (P - origin) @ invB.T
            ij = np.rint(coeff).astype(int)
            pred = origin + ij @ B.T
            resid = np.linalg.norm(P - pred, axis=1)
            tol = max(3.0, LATTICE_RESIDUAL_FRAC * min(np.linalg.norm(a), np.linalg.norm(b)))
            good = resid <= tol
            if int(np.count_nonzero(good)) < LATTICE_MIN_SEEDS:
                ok = False
                break
            # Refine affine origin/basis by least squares with fixed integer indices.
            X = np.column_stack([
                np.ones(int(np.count_nonzero(good))),
                ij[good, 0].astype(float),
                ij[good, 1].astype(float),
            ])
            coef_x, *_ = np.linalg.lstsq(X, P[good, 0], rcond=None)
            coef_y, *_ = np.linalg.lstsq(X, P[good, 1], rcond=None)
            origin = np.array([coef_x[0], coef_y[0]], dtype=float)
            a = np.array([coef_x[1], coef_y[1]], dtype=float)
            b = np.array([coef_x[2], coef_y[2]], dtype=float)
        if not ok:
            continue

        B = np.column_stack([a, b])
        if abs(float(np.linalg.det(B))) < 1e-5:
            continue
        invB = np.linalg.inv(B)
        coeff = (P - origin) @ invB.T
        ij = np.rint(coeff).astype(int)
        pred = origin + ij @ B.T
        resid = np.linalg.norm(P - pred, axis=1)
        pitch_min = min(float(np.linalg.norm(a)), float(np.linalg.norm(b)))
        tol = max(3.0, LATTICE_RESIDUAL_FRAC * pitch_min)
        good = resid <= tol
        if int(np.count_nonzero(good)) < LATTICE_MIN_SEEDS:
            continue
        uniq = sorted(set((int(i), int(j)) for i, j in ij[good]))
        if len(uniq) < LATTICE_MIN_SEEDS:
            continue
        ii = [q[0] for q in uniq]
        jj = [q[1] for q in uniq]
        imin, imax = min(ii), max(ii)
        jmin, jmax = min(jj), max(jj)
        expected = (imax - imin + 1) * (jmax - jmin + 1)
        if expected <= 0 or expected > max_expected:
            continue
        occupancy = len(uniq) / expected
        support_frac = float(np.count_nonzero(good)) / max(1, len(P))
        rms = float(np.sqrt(np.mean(resid[good] ** 2)))
        # Fundamental lattice bases explain many points with high occupancy and low residual.
        score = 3.2 * support_frac + 2.3 * occupancy + 0.012 * min(len(uniq), 80) - 0.020 * rms - 0.0015 * expected
        rec = {
            "score": score,
            "origin": origin,
            "a": a,
            "b": b,
            "good_mask": good,
            "ij": ij,
            "residuals": resid,
            "i_min": int(imin), "i_max": int(imax),
            "j_min": int(jmin), "j_max": int(jmax),
            "occupancy": float(occupancy),
            "supported_sites": int(len(uniq)),
            "expected_sites": int(expected),
            "rms": rms,
        }
        if best is None or score > float(best["score"]):
            best = rec
    return best


def _orient_slot_lattice_v14(model: Dict[str, object], centers: np.ndarray, image_shape: Tuple[int, int]) -> Dict[str, object]:
    """Choose one lattice coordinate as the four-level slot-row direction."""
    m = dict(model)
    ij = np.asarray(m["ij"], dtype=int)
    good = np.asarray(m["good_mask"], dtype=bool)
    a = np.asarray(m["a"], dtype=float)
    b = np.asarray(m["b"], dtype=float)

    vals0 = np.unique(ij[good, 0])
    vals1 = np.unique(ij[good, 1])
    # Prefer an axis whose observed level count is near four; verticality is a tie-breaker.
    def axis_score(vals: np.ndarray, vec: np.ndarray) -> float:
        verticality = abs(float(vec[1])) / max(1e-9, float(np.linalg.norm(vec)))
        return -1.5 * abs(len(vals) - 4) + 0.7 * verticality
    row_axis = 0 if axis_score(vals0, a) > axis_score(vals1, b) else 1

    if row_axis == 0:
        # New i = old j (group direction), new j = old i (slot-row direction).
        m["a"], m["b"] = b.copy(), a.copy()
        m["ij"] = ij[:, [1, 0]].copy()
        m["i_min"], m["i_max"] = int(m["j_min"]), int(m["j_max"])
        m["j_min"], m["j_max"] = int(model["i_min"]), int(model["i_max"])
    else:
        m["a"], m["b"] = a.copy(), b.copy()
        m["ij"] = ij.copy()

    # Make +j go downward in the image so j_min..j_max maps naturally to S1..S4.
    if float(np.asarray(m["b"])[1]) < 0:
        m["b"] = -np.asarray(m["b"], dtype=float)
        mij = np.asarray(m["ij"], dtype=int)
        mij[:, 1] *= -1
        m["ij"] = mij
        old_min, old_max = int(m["j_min"]), int(m["j_max"])
        m["j_min"], m["j_max"] = -old_max, -old_min

    mij = np.asarray(m["ij"], dtype=int)
    good = np.asarray(m["good_mask"], dtype=bool)
    jvals = mij[good, 1]
    if len(jvals) == 0:
        m["used"] = False
        return m

    # Find the best four consecutive row indices. If fewer than four are observed,
    # extend only enough to make four; all final sites still require image evidence.
    counts = Counter(int(v) for v in jvals)
    jmin_obs, jmax_obs = int(np.min(jvals)), int(np.max(jvals))
    starts = range(jmin_obs - 3, jmax_obs + 1)
    best = None
    H, W = image_shape
    origin = np.asarray(m["origin"], dtype=float)
    avec = np.asarray(m["a"], dtype=float)
    bvec = np.asarray(m["b"], dtype=float)
    for st in starts:
        en = st + 3
        support = sum(counts.get(j, 0) for j in range(st, en + 1))
        outside = sum(v for j, v in counts.items() if j < st or j > en)
        # Prefer four levels that explain most observations and stay mostly in image.
        probe_i = int(round((int(m["i_min"]) + int(m["i_max"])) / 2))
        ys = [(origin + probe_i * avec + j * bvec)[1] for j in range(st, en + 1)]
        in_image = sum(0 <= y < H for y in ys)
        score = 2.0 * support - 2.5 * outside + 0.25 * in_image
        if best is None or score > best[0]:
            best = (score, st, en)
    if best is None:
        m["used"] = False
        return m
    _, jst, jen = best
    m["j_min"], m["j_max"] = int(jst), int(jen)

    # Recompute i bounds using seeds that fall in the chosen 4-row window.
    in_rows = good & (mij[:, 1] >= jst) & (mij[:, 1] <= jen)
    if int(np.count_nonzero(in_rows)) >= LATTICE_MIN_SEEDS:
        m["i_min"] = int(np.min(mij[in_rows, 0]))
        m["i_max"] = int(np.max(mij[in_rows, 0]))
    return m


def infer_lattice_v14(seeds: List[ShapeDesc], image_shape: Tuple[int, int], pattern: str) -> Dict[str, object]:
    meta: Dict[str, object] = {"used": False, "lattice_type": "none"}
    if len(seeds) < LATTICE_MIN_SEEDS:
        return meta
    centers = np.asarray([d.center for d in seeds], dtype=float)
    vecs = _lattice_vector_candidates_v14(seeds)
    if len(vecs) < 2:
        return meta

    best = None
    for i in range(len(vecs)):
        for j in range(i + 1, len(vecs)):
            a, b = vecs[i], vecs[j]
            la, lb = np.linalg.norm(a), np.linalg.norm(b)
            if la <= 1 or lb <= 1:
                continue
            ang = math.degrees(math.acos(np.clip(float(np.dot(a, b) / (la * lb)), -1.0, 1.0)))
            if ang < LATTICE_MIN_BASIS_ANGLE_DEG or ang > LATTICE_MAX_BASIS_ANGLE_DEG:
                continue
            fit = _fit_affine_lattice_basis_v14(centers, a, b, LATTICE_MAX_EXPECTED_OBJECTS)
            if fit is None:
                continue
            # For slots, favor a basis that naturally yields a four-level coordinate axis.
            if pattern == "slot":
                ij = np.asarray(fit["ij"], dtype=int)
                good = np.asarray(fit["good_mask"], dtype=bool)
                n0 = len(np.unique(ij[good, 0]))
                n1 = len(np.unique(ij[good, 1]))
                fit["score"] = float(fit["score"]) + 0.45 * max(0.0, 1.0 - min(abs(n0 - 4), abs(n1 - 4)) / 4.0)
            if best is None or float(fit["score"]) > float(best["score"]):
                best = fit
    if best is None:
        return meta

    if pattern == "slot":
        best = _orient_slot_lattice_v14(best, centers, image_shape)
        if not best.get("used", True):
            return meta

    # Keep only plausible models after orientation.
    expected = (int(best["i_max"]) - int(best["i_min"]) + 1) * (int(best["j_max"]) - int(best["j_min"]) + 1)
    if expected <= 0 or expected > LATTICE_MAX_EXPECTED_OBJECTS:
        return meta
    # Occupancy is recalculated for chosen bounds.
    mij = np.asarray(best["ij"], dtype=int)
    good = np.asarray(best["good_mask"], dtype=bool)
    sites = set()
    for q in mij[good]:
        ii, jj = int(q[0]), int(q[1])
        if int(best["i_min"]) <= ii <= int(best["i_max"]) and int(best["j_min"]) <= jj <= int(best["j_max"]):
            sites.add((ii, jj))
    occ = len(sites) / max(1, expected)
    if occ < LATTICE_MIN_OCCUPANCY:
        return meta

    a = np.asarray(best["a"], dtype=float)
    b = np.asarray(best["b"], dtype=float)
    la, lb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    angle = math.degrees(math.acos(np.clip(float(np.dot(a, b) / max(1e-9, la * lb)), -1.0, 1.0)))
    med_w = float(np.median([max(1.0, _contour_xy_bounds(d)[1] - _contour_xy_bounds(d)[0]) for d in seeds]))
    med_h = float(np.median([max(1.0, _contour_xy_bounds(d)[3] - _contour_xy_bounds(d)[2]) for d in seeds]))
    best.update({
        "used": True,
        "lattice_type": "affine_2d",
        "occupancy": float(occ),
        "supported_sites": int(len(sites)),
        "expected_sites": int(expected),
        "basis_a_len": la,
        "basis_b_len": lb,
        "basis_angle_deg": float(angle),
        "median_w": med_w,
        "median_h": med_h,
    })
    return best


def _lattice_point_v14(model: Dict[str, object], i: int, j: int) -> np.ndarray:
    return (
        np.asarray(model["origin"], dtype=float)
        + float(i) * np.asarray(model["a"], dtype=float)
        + float(j) * np.asarray(model["b"], dtype=float)
    )


def lattice_sites_v14(
    model: Dict[str, object],
    image_shape: Tuple[int, int],
    expand_rings: int = 0,
    pattern: str = "via",
) -> List[Dict[str, object]]:
    if not model.get("used", False):
        return []
    H, W = image_shape
    imin, imax = int(model["i_min"]), int(model["i_max"])
    jmin, jmax = int(model["j_min"]), int(model["j_max"])
    if expand_rings > 0:
        imin -= expand_rings
        imax += expand_rings
        # Slot has exactly four row types. Extend groups/columns, not the 4-row axis.
        if pattern != "slot":
            jmin -= expand_rings
            jmax += expand_rings
    sites: List[Dict[str, object]] = []
    for j in range(jmin, jmax + 1):
        for i in range(imin, imax + 1):
            p = _lattice_point_v14(model, i, j)
            x, y = float(p[0]), float(p[1])
            if 2 <= x <= W - 3 and 2 <= y <= H - 3:
                sites.append({
                    "i": int(i), "j": int(j), "x": x, "y": y,
                    "row_id": int(j - int(model["j_min"]) + 1) if pattern == "slot" else None,
                    "is_internal": bool(
                        int(model["i_min"]) <= i <= int(model["i_max"])
                        and int(model["j_min"]) <= j <= int(model["j_max"])
                    ),
                })
    return sites


def _nearest_lattice_site_v14(model: Dict[str, object], point: np.ndarray) -> Tuple[int, int, float]:
    B = np.column_stack([np.asarray(model["a"], dtype=float), np.asarray(model["b"], dtype=float)])
    if abs(float(np.linalg.det(B))) < 1e-8:
        return 0, 0, np.inf
    q = np.linalg.inv(B) @ (np.asarray(point, dtype=float) - np.asarray(model["origin"], dtype=float))
    ij = np.rint(q).astype(int)
    pred = _lattice_point_v14(model, int(ij[0]), int(ij[1]))
    return int(ij[0]), int(ij[1]), float(np.linalg.norm(np.asarray(point, dtype=float) - pred))


def _lattice_match_tol_v14(model: Dict[str, object], size_ref: float) -> float:
    pitch = min(float(model.get("basis_a_len", np.nan)), float(model.get("basis_b_len", np.nan)))
    vals = [4.0, 0.75 * max(2.0, size_ref)]
    if np.isfinite(pitch) and pitch > 0:
        vals.append(LATTICE_MATCH_FRAC * pitch)
    return float(max(vals))


def _select_best_near_lattice_site_v14(
    pool: List[ShapeDesc],
    used_ids: set,
    expected_xy: Tuple[float, float],
    tol: float,
    pattern: str,
    expected_area: float,
    slot_area_limits: Optional[Tuple[float, float]] = None,
    final_relaxed: bool = False,
) -> Optional[ShapeDesc]:
    ex = np.asarray(expected_xy, dtype=float)
    cand: List[Tuple[float, ShapeDesc]] = []
    for d in pool:
        if id(d) in used_ids:
            continue
        dist = float(np.linalg.norm(d.center - ex))
        if dist > tol:
            continue
        pa = contour_pixel_count(d.contour)
        if pa < MIN_COMPONENT_AREA_PX:
            continue
        if pattern == "via":
            if final_relaxed:
                if not (
                    d.circularity >= VIA_FINAL_MIN_CIRCULARITY
                    and d.aspect <= VIA_FINAL_MAX_AXIS_RATIO
                    and d.solidity >= VIA_FINAL_MIN_SOLIDITY
                ):
                    continue
            elif not _via_grid_relaxed_ok_v12(d):
                continue
        else:
            if slot_area_limits is not None:
                lo, hi = slot_area_limits
                if pa < lo or pa > hi:
                    continue
            if final_relaxed:
                reg = slot_regularity_metrics(d)
                if not (
                    reg["regularity_solidity"] >= SLOT_FINAL_MIN_SOLIDITY
                    and reg["regularity_convexity"] >= SLOT_FINAL_MIN_CONVEXITY
                    and reg["regularity_template_iou"] >= SLOT_FINAL_MIN_TEMPLATE_IOU
                ):
                    continue
            else:
                ok, _, _, _ = _slot_regularity_ok_v12(d, relaxed=True)
                if not ok:
                    continue
        score = _local_candidate_score_v12(d, expected_xy, expected_area, pattern)
        score += 0.25 * dist / max(3.0, tol)
        if _source_of(d) == "primary":
            score -= 0.12
        cand.append((score, d))
    if not cand:
        return None
    cand.sort(key=lambda z: z[0])
    return cand[0][1]


def _local_contrast_at_site_v14(gray: np.ndarray, xy: Tuple[float, float], w: float, h: float) -> float:
    H, W = gray.shape
    cx, cy = xy
    hw = max(4, int(round(0.65 * max(4.0, w))))
    hh = max(4, int(round(0.65 * max(4.0, h))))
    x0, x1 = max(0, int(cx - hw)), min(W, int(cx + hw + 1))
    y0, y1 = max(0, int(cy - hh)), min(H, int(cy + hh + 1))
    roi = gray[y0:y1, x0:x1].astype(np.float32)
    if roi.size < 25:
        return 0.0
    band = max(1, min(4, min(roi.shape) // 6))
    border = np.concatenate([roi[:band].ravel(), roi[-band:].ravel(), roi[:, :band].ravel(), roi[:, -band:].ravel()])
    ch, cw = roi.shape
    core = roi[max(0, ch//2-hh//3):min(ch, ch//2+hh//3+1), max(0, cw//2-hw//3):min(cw, cw//2+hw//3+1)]
    if border.size == 0 or core.size == 0:
        return 0.0
    bg = float(np.median(border))
    if FEATURE_POLARITY.lower() == "dark":
        target = float(np.percentile(core, 35))
        return bg - target
    target = float(np.percentile(core, 65))
    return target - bg


def recover_local_object_v14(
    gray: np.ndarray,
    expected_xy: Tuple[float, float],
    expected_w: float,
    expected_h: float,
    pattern: str,
    nearest_pitch: float,
    expected_area: float,
    slot_area_limits: Optional[Tuple[float, float]] = None,
    final_relaxed: bool = False,
) -> Optional[ShapeDesc]:
    """Site-guided local recovery. In final_relaxed mode shape thresholds are looser."""
    H, W = gray.shape
    cx, cy = expected_xy
    half_w = max(10.0, ARRAY_LOCAL_ROI_SIZE_FACTOR * max(4.0, expected_w) / 2.0)
    half_h = max(10.0, ARRAY_LOCAL_ROI_SIZE_FACTOR * max(4.0, expected_h) / 2.0)
    if np.isfinite(nearest_pitch) and nearest_pitch > 0:
        half_w = min(max(half_w, expected_w), 0.48 * nearest_pitch)
        half_h = min(max(half_h, expected_h), 0.48 * nearest_pitch)
    x0 = max(0, int(math.floor(cx - half_w)))
    x1 = min(W, int(math.ceil(cx + half_w + 1)))
    y0 = max(0, int(math.floor(cy - half_h)))
    y1 = min(H, int(math.ceil(cy + half_h + 1)))
    if x1 - x0 < 8 or y1 - y0 < 8:
        return None

    if final_relaxed and LATTICE_FINAL_REQUIRE_LOCAL_CONTRAST:
        contrast = _local_contrast_at_site_v14(gray, expected_xy, expected_w, expected_h)
        if contrast < LATTICE_FINAL_MIN_CONTRAST_GRAY:
            return None

    roi = gray[y0:y1, x0:x1]
    masks: Dict[str, np.ndarray] = {"local_primary": segment_features(roi)}
    masks.update(enhanced_masks_v12(roi))

    if final_relaxed:
        # Percentile masks are cheap inside a tiny ROI and recover very blurred structures.
        percs = FINAL_LOCAL_PERCENTILES_DARK if FEATURE_POLARITY.lower() == "dark" else FINAL_LOCAL_PERCENTILES_BRIGHT
        for q in percs:
            thr = float(np.percentile(roi, q))
            if FEATURE_POLARITY.lower() == "dark":
                mm = (roi.astype(np.float32) <= thr).astype(np.uint8) * 255
            else:
                mm = (roi.astype(np.float32) >= thr).astype(np.uint8) * 255
            masks[f"final_pct_{q:g}"] = _postprocess_binary_v12(mm, 1)

    candidates: List[ShapeDesc] = []
    tol = FINAL_CENTER_TOL_FRAC * max(half_w, half_h) if final_relaxed else ARRAY_LOCAL_CENTER_TOL_FRAC * max(half_w, half_h)
    for source, mm in masks.items():
        contours, _ = cv2.findContours(mm, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        for c in contours:
            pa0 = contour_pixel_count(c)
            if pa0 < MIN_COMPONENT_AREA_PX:
                continue
            bx, by, bw, bh = cv2.boundingRect(c)
            # A contour touching the local crop is likely a clipped background region.
            if bx <= 0 or by <= 0 or bx + bw >= roi.shape[1] or by + bh >= roi.shape[0]:
                continue
            cg = c.astype(np.int32).copy()
            cg[:, 0, 0] += x0
            cg[:, 0, 1] += y0
            d = describe_contour(cg)
            if d is None:
                continue
            if float(np.linalg.norm(d.center - np.asarray(expected_xy, dtype=float))) > tol:
                continue
            pa = contour_pixel_count(d.contour)
            if pattern == "slot" and slot_area_limits is not None:
                lo, hi = slot_area_limits
                if pa < lo or pa > hi:
                    continue
            if pattern == "via":
                if final_relaxed:
                    if not (
                        d.circularity >= VIA_FINAL_MIN_CIRCULARITY
                        and d.aspect <= VIA_FINAL_MAX_AXIS_RATIO
                        and d.solidity >= VIA_FINAL_MIN_SOLIDITY
                    ):
                        continue
                elif not _via_grid_relaxed_ok_v12(d):
                    continue
            else:
                if final_relaxed:
                    reg = slot_regularity_metrics(d)
                    if not (
                        reg["regularity_solidity"] >= SLOT_FINAL_MIN_SOLIDITY
                        and reg["regularity_convexity"] >= SLOT_FINAL_MIN_CONVEXITY
                        and reg["regularity_template_iou"] >= SLOT_FINAL_MIN_TEMPLATE_IOU
                    ):
                        continue
                else:
                    ok, _, _, _ = _slot_regularity_ok_v12(d, relaxed=True)
                    if not ok:
                        continue
            _mark_source(d, "lattice_final_recovered" if final_relaxed else "lattice_recovered")
            candidates.append(d)
    if not candidates:
        return None
    return min(candidates, key=lambda d: _local_candidate_score_v12(d, expected_xy, expected_area, pattern))


def _refine_lattice_from_selected_v14(model: Dict[str, object], selected: List[ShapeDesc]) -> Dict[str, object]:
    """Refit origin/a/b from selected centers while keeping their nearest integer lattice IDs."""
    if not model.get("used", False) or len(selected) < LATTICE_MIN_SEEDS:
        return model
    rows = []
    pts = []
    for d in selected:
        i, j, dist = _nearest_lattice_site_v14(model, d.center)
        if not np.isfinite(dist):
            continue
        rows.append((i, j))
        pts.append(d.center)
    if len(rows) < LATTICE_MIN_SEEDS:
        return model
    ij = np.asarray(rows, dtype=float)
    P = np.asarray(pts, dtype=float)
    X = np.column_stack([np.ones(len(ij)), ij[:, 0], ij[:, 1]])
    cx, *_ = np.linalg.lstsq(X, P[:, 0], rcond=None)
    cy, *_ = np.linalg.lstsq(X, P[:, 1], rcond=None)
    m = dict(model)
    m["origin"] = np.array([cx[0], cy[0]], dtype=float)
    m["a"] = np.array([cx[1], cy[1]], dtype=float)
    m["b"] = np.array([cx[2], cy[2]], dtype=float)
    m["basis_a_len"] = float(np.linalg.norm(m["a"]))
    m["basis_b_len"] = float(np.linalg.norm(m["b"]))
    la, lb = m["basis_a_len"], m["basis_b_len"]
    m["basis_angle_deg"] = math.degrees(math.acos(np.clip(float(np.dot(m["a"], m["b"]) / max(1e-9, la * lb)), -1.0, 1.0)))
    return m


def _lattice_model_score_v14(model: Dict[str, object]) -> float:
    if not model.get("used", False):
        return -1e9
    return float(model.get("occupancy", 0.0)) + 0.012 * min(int(model.get("supported_sites", 0)), 80) - 0.01 * float(model.get("rms", 0.0))


def _slot_row_dims_v14(
    model: Dict[str, object],
    seeds: List[ShapeDesc],
    mean_area: float,
) -> Dict[int, Tuple[float, float, float]]:
    med_w = float(model.get("median_w", 10.0))
    med_h = float(model.get("median_h", 10.0))
    if np.isfinite(mean_area) and mean_area > 0:
        area_ref = float(mean_area)
    elif seeds:
        area_ref = float(np.median([contour_pixel_count(d.contour) for d in seeds]))
    else:
        area_ref = max(30.0, med_w * med_h * 0.6)
    out: Dict[int, Tuple[float, float, float]] = {}
    for j in range(int(model["j_min"]), int(model["j_max"]) + 1):
        rid = j - int(model["j_min"]) + 1
        near = []
        for d in seeds:
            i0, j0, dist = _nearest_lattice_site_v14(model, d.center)
            if j0 == j and dist <= _lattice_match_tol_v14(model, max(med_w, 0.25 * med_h)):
                near.append(d)
        if len(near) >= 2:
            ws = [max(1.0, _contour_xy_bounds(d)[1] - _contour_xy_bounds(d)[0]) for d in near]
            hs = [max(1.0, _contour_xy_bounds(d)[3] - _contour_xy_bounds(d)[2]) for d in near]
            aa = [contour_pixel_count(d.contour) for d in near]
            out[rid] = (float(np.median(ws)), float(np.median(hs)), float(np.median(aa)))
        else:
            out[rid] = (med_w, med_h, area_ref)
    return out



# ============================================================
# V15: CENTER-BASED RECURSIVE EQUAL-SPACING SEARCH
# ============================================================
def _model_axis_pitch_v15(model: Dict[str, object], axis: str) -> float:
    """Find the shortest lattice-vector combination close to horizontal/vertical."""
    if not model.get("used", False):
        return np.nan
    a = np.asarray(model.get("a", [np.nan, np.nan]), dtype=float)
    b = np.asarray(model.get("b", [np.nan, np.nan]), dtype=float)
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        return np.nan
    tol = math.radians(CENTER_DIRECTION_ANGLE_TOL_DEG)
    vals = []
    for i in range(-3, 4):
        for j in range(-3, 4):
            if i == 0 and j == 0:
                continue
            v = i * a + j * b
            L = float(np.linalg.norm(v))
            if L < CENTER_PITCH_MIN_PX:
                continue
            ang = abs(math.atan2(float(v[1]), float(v[0]))) % math.pi
            if axis == "x":
                d = min(ang, math.pi - ang)
                comp = abs(float(v[0]))
            else:
                d = abs(ang - math.pi / 2.0)
                comp = abs(float(v[1]))
            if d <= tol and comp >= CENTER_PITCH_MIN_PX:
                vals.append((L, comp))
    if not vals:
        return np.nan
    vals.sort(key=lambda z: z[0])
    return float(vals[0][1])


def _axis_pitch_from_centers_v15(selected: List[ShapeDesc], axis: str, image_shape: Tuple[int, int]) -> float:
    """Nearest repeated horizontal/vertical CENTER spacing; one nearest value per center."""
    if len(selected) < 2:
        return np.nan
    P = np.asarray([d.center for d in selected], dtype=float)
    tan_tol = math.tan(math.radians(CENTER_DIRECTION_ANGLE_TOL_DEG))
    nearest = []
    for i in range(len(P)):
        best = np.inf
        for j in range(len(P)):
            if i == j:
                continue
            dx = abs(float(P[j, 0] - P[i, 0]))
            dy = abs(float(P[j, 1] - P[i, 1]))
            if axis == "x":
                if dx < CENTER_PITCH_MIN_PX or dy > tan_tol * max(dx, 1e-9):
                    continue
                val = dx
                max_allowed = image_shape[1] * CENTER_PITCH_MAX_IMAGE_FRAC
            else:
                if dy < CENTER_PITCH_MIN_PX or dx > tan_tol * max(dy, 1e-9):
                    continue
                val = dy
                max_allowed = image_shape[0] * CENTER_PITCH_MAX_IMAGE_FRAC
            if val <= max_allowed:
                best = min(best, val)
        if np.isfinite(best):
            nearest.append(best)
    if not nearest:
        return np.nan
    arr = np.asarray(nearest, dtype=float)
    med = float(np.median(arr))
    # Remove obvious multiples/outliers and re-estimate from the central group.
    keep = arr[(arr >= 0.65 * med) & (arr <= 1.55 * med)]
    return float(np.median(keep)) if len(keep) else med


def estimate_center_pitches_v15(
    selected: List[ShapeDesc], model: Dict[str, object], image_shape: Tuple[int, int]
) -> Tuple[float, float]:
    """Estimate horizontal and vertical center pitch; affine model is only a fallback."""
    hp = _axis_pitch_from_centers_v15(selected, "x", image_shape)
    vp = _axis_pitch_from_centers_v15(selected, "y", image_shape)
    if not np.isfinite(hp):
        hp = _model_axis_pitch_v15(model, "x")
    if not np.isfinite(vp):
        vp = _model_axis_pitch_v15(model, "y")
    # Last fallback: use the short affine basis length. This keeps the recursion usable
    # on strongly rotated/rhombic arrays where no pair is nearly horizontal/vertical.
    basis = [float(model.get("basis_a_len", np.nan)), float(model.get("basis_b_len", np.nan))]
    basis = [x for x in basis if np.isfinite(x) and x >= CENTER_PITCH_MIN_PX]
    fallback = min(basis) if basis else np.nan
    if not np.isfinite(hp):
        hp = fallback
    if not np.isfinite(vp):
        vp = fallback
    return float(hp), float(vp)


def _center_ultra_shape_ok_v15(d: ShapeDesc, pattern: str, slot_area_limits: Optional[Tuple[float, float]]) -> bool:
    pa = contour_pixel_count(d.contour)
    if pa < MIN_COMPONENT_AREA_PX:
        return False
    if pattern == "via":
        return (
            d.circularity >= VIA_CENTER_MIN_CIRCULARITY
            and d.aspect <= VIA_CENTER_MAX_AXIS_RATIO
            and d.solidity >= VIA_CENTER_MIN_SOLIDITY
        )
    if slot_area_limits is not None:
        lo, hi = slot_area_limits
        if pa < max(SLOT_MIN_AREA_PX, lo) or pa > hi:
            return False
    reg = slot_regularity_metrics(d)
    return (
        reg["regularity_solidity"] >= SLOT_CENTER_MIN_SOLIDITY
        and reg["regularity_convexity"] >= SLOT_CENTER_MIN_CONVEXITY
        and reg["regularity_template_iou"] >= SLOT_CENTER_MIN_TEMPLATE_IOU
    )


def _select_best_near_center_v15(
    pool: List[ShapeDesc],
    used_ids: set,
    expected_xy: Tuple[float, float],
    tol: float,
    pattern: str,
    expected_area: float,
    slot_area_limits: Optional[Tuple[float, float]],
) -> Optional[ShapeDesc]:
    ex = np.asarray(expected_xy, dtype=float)
    cand = []
    for d in pool:
        if id(d) in used_ids:
            continue
        dist = float(np.linalg.norm(d.center - ex))
        if dist > tol:
            continue
        if not _center_ultra_shape_ok_v15(d, pattern, slot_area_limits):
            continue
        score = _local_candidate_score_v12(d, expected_xy, expected_area, pattern)
        score += 0.35 * dist / max(3.0, tol)
        if _source_of(d) == "primary":
            score -= 0.10
        cand.append((score, d))
    if not cand:
        return None
    cand.sort(key=lambda q: q[0])
    return cand[0][1]


def recover_local_object_v15_center(
    gray: np.ndarray,
    expected_xy: Tuple[float, float],
    expected_w: float,
    expected_h: float,
    pattern: str,
    nearest_pitch: float,
    expected_area: float,
    slot_area_limits: Optional[Tuple[float, float]],
) -> Optional[ShapeDesc]:
    """Ultra-relaxed recovery used ONLY at recursively predicted center locations."""
    H, W = gray.shape
    cx, cy = expected_xy
    half_w = max(11.0, 1.35 * ARRAY_LOCAL_ROI_SIZE_FACTOR * max(4.0, expected_w) / 2.0)
    half_h = max(11.0, 1.35 * ARRAY_LOCAL_ROI_SIZE_FACTOR * max(4.0, expected_h) / 2.0)
    if np.isfinite(nearest_pitch) and nearest_pitch > 0:
        cap = 0.48 * nearest_pitch
        half_w = min(max(half_w, expected_w), cap)
        half_h = min(max(half_h, expected_h), cap)
    x0 = max(0, int(math.floor(cx - half_w)))
    x1 = min(W, int(math.ceil(cx + half_w + 1)))
    y0 = max(0, int(math.floor(cy - half_h)))
    y1 = min(H, int(math.ceil(cy + half_h + 1)))
    if x1 - x0 < 8 or y1 - y0 < 8:
        return None

    contrast = _local_contrast_at_site_v14(gray, expected_xy, expected_w, expected_h)
    if contrast < CENTER_RECURSIVE_LOCAL_CONTRAST_GRAY:
        return None

    roi = gray[y0:y1, x0:x1]
    masks: Dict[str, np.ndarray] = {"center_primary": segment_features(roi)}
    masks.update(enhanced_masks_v12(roi))
    percs = CENTER_LOCAL_PERCENTILES_DARK if FEATURE_POLARITY.lower() == "dark" else CENTER_LOCAL_PERCENTILES_BRIGHT
    for q in percs:
        thr = float(np.percentile(roi, q))
        if FEATURE_POLARITY.lower() == "dark":
            mm = (roi.astype(np.float32) <= thr).astype(np.uint8) * 255
        else:
            mm = (roi.astype(np.float32) >= thr).astype(np.uint8) * 255
        masks[f"center_pct_{q:g}"] = _postprocess_binary_v12(mm, 1)

    candidates = []
    tol = CENTER_RECURSIVE_CENTER_TOL_FRAC * max(half_w, half_h)
    for _, mm in masks.items():
        contours, _ = cv2.findContours(mm, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        for c in contours:
            if contour_pixel_count(c) < MIN_COMPONENT_AREA_PX:
                continue
            bx, by, bw, bh = cv2.boundingRect(c)
            if bx <= 0 or by <= 0 or bx + bw >= roi.shape[1] or by + bh >= roi.shape[0]:
                continue
            cg = c.astype(np.int32).copy()
            cg[:, 0, 0] += x0
            cg[:, 0, 1] += y0
            d = describe_contour(cg)
            if d is None:
                continue
            if float(np.linalg.norm(d.center - np.asarray(expected_xy, dtype=float))) > tol:
                continue
            if not _center_ultra_shape_ok_v15(d, pattern, slot_area_limits):
                continue
            _mark_source(d, "center_recursive_recovered")
            candidates.append(d)
    if not candidates:
        return None
    return min(candidates, key=lambda d: _local_candidate_score_v12(d, expected_xy, expected_area, pattern))


def _slot_row_id_from_model_v15(model: Dict[str, object], center: np.ndarray) -> Optional[int]:
    if not model.get("used", False):
        return None
    B = np.column_stack([np.asarray(model["a"], dtype=float), np.asarray(model["b"], dtype=float)])
    if abs(float(np.linalg.det(B))) < 1e-8:
        return None
    q = np.linalg.inv(B) @ (np.asarray(center, dtype=float) - np.asarray(model["origin"], dtype=float))
    j = int(round(float(q[1])))
    rid = j - int(model["j_min"]) + 1
    if 1 <= rid <= 4:
        return rid
    # For a slightly off-model recovered center, assign the nearest of the four
    # row-center coordinates instead of discarding an otherwise valid object.
    js = list(range(int(model["j_min"]), int(model["j_max"]) + 1))
    if not js:
        return None
    jbest = min(js, key=lambda jj: abs(float(q[1]) - jj))
    return int(jbest - int(model["j_min"]) + 1)


def assign_slot_rows_v15(model: Dict[str, object], selected: List[ShapeDesc]) -> List[int]:
    if not selected:
        return []
    if model.get("used", False):
        out = [_slot_row_id_from_model_v15(model, d.center) for d in selected]
        if all(r is not None and 1 <= int(r) <= 4 for r in out):
            return [int(r) for r in out]
    if len(selected) >= 4:
        labs = cluster_1d_into_four_rows(np.asarray([d.center[1] for d in selected], dtype=np.float32))
        return [int(x) for x in labs]
    return []


def recursive_center_recovery_v15(
    gray: np.ndarray,
    pool: List[ShapeDesc],
    selected: List[ShapeDesc],
    pattern: str,
    model: Dict[str, object],
    image_name: str,
    slot_area_limits: Optional[Tuple[float, float]] = None,
) -> Tuple[List[ShapeDesc], Dict[str, object]]:
    """
    Recursively expand from ACCEPTED OBJECT CENTERS in 8 directions.

    H/V predictions use the measured center-to-center H/V pitch. For 45-degree
    predictions, dx=dy=horizontal_pitch/sqrt(2), which implements the requested
    horizontal_pitch/1.414 center offset along both image axes.
    """
    meta = {
        "center_recursive_used": False,
        "center_horizontal_pitch_px": np.nan,
        "center_vertical_pitch_px": np.nan,
        "center_diagonal_component_px": np.nan,
        "center_recursive_attempts": 0,
        "center_recursive_recovered": 0,
        "center_recursive_boundary_hits": 0,
    }
    if len(selected) < 2:
        return selected, meta
    hp, vp = estimate_center_pitches_v15(selected, model, gray.shape)
    if not np.isfinite(hp) and not np.isfinite(vp):
        return selected, meta
    if not np.isfinite(hp):
        hp = vp
    if not np.isfinite(vp):
        vp = hp
    diag = hp / CENTER_DIAGONAL_DIVISOR if np.isfinite(hp) else np.nan
    if hp < CENTER_PITCH_MIN_PX or vp < CENTER_PITCH_MIN_PX:
        return selected, meta

    meta.update({
        "center_recursive_used": True,
        "center_horizontal_pitch_px": float(hp),
        "center_vertical_pitch_px": float(vp),
        "center_diagonal_component_px": float(diag),
    })
    H, W = gray.shape
    nearest_pitch = float(min(hp, vp, hp if not np.isfinite(diag) else max(diag, 1.0)))
    med_w = float(np.median([max(1.0, _contour_xy_bounds(d)[1] - _contour_xy_bounds(d)[0]) for d in selected]))
    med_h = float(np.median([max(1.0, _contour_xy_bounds(d)[3] - _contour_xy_bounds(d)[2]) for d in selected]))
    med_area = float(np.median([contour_pixel_count(d.contour) for d in selected]))
    tol = max(4.0, CENTER_RECURSIVE_SITE_TOL_FRAC * min(hp, vp))
    occ_tol = max(4.0, 0.55 * tol)
    visit_cell = max(3.0, 0.55 * tol)

    directions = [
        np.array([ hp, 0.0]), np.array([-hp, 0.0]),
        np.array([0.0,  vp]), np.array([0.0, -vp]),
    ]
    if np.isfinite(diag) and diag >= 2.0:
        directions.extend([
            np.array([ diag,  diag]), np.array([ diag, -diag]),
            np.array([-diag,  diag]), np.array([-diag, -diag]),
        ])

    from collections import deque
    queue = deque(selected)
    used_ids = {id(d) for d in selected}
    visited = set()
    attempts = 0
    recovered = 0
    boundary_hits = 0

    def occupied(xy: np.ndarray) -> bool:
        return any(float(np.linalg.norm(d.center - xy)) <= occ_tol for d in selected)

    while queue and attempts < CENTER_RECURSIVE_MAX_ATTEMPTS and recovered < CENTER_RECURSIVE_MAX_RECOVERED:
        base = queue.popleft()
        for dv in directions:
            expected = np.asarray(base.center, dtype=float) + dv
            x, y = float(expected[0]), float(expected[1])
            margin_x = max(CENTER_RECURSIVE_BOUNDARY_MARGIN_PX, 0.55 * med_w)
            margin_y = max(CENTER_RECURSIVE_BOUNDARY_MARGIN_PX, 0.55 * med_h)
            if x <= margin_x or x >= W - 1 - margin_x or y <= margin_y or y >= H - 1 - margin_y:
                boundary_hits += 1
                continue
            key = (int(round(x / visit_cell)), int(round(y / visit_cell)))
            if key in visited:
                continue
            visited.add(key)
            if occupied(expected):
                continue
            attempts += 1
            if attempts == 1 or attempts % CENTER_RECURSIVE_PROGRESS_EVERY == 0:
                print(
                    f"      [{image_name}] V15 center recursive search {attempts}/{CENTER_RECURSIVE_MAX_ATTEMPTS}; recovered={recovered} ...",
                    flush=True,
                )

            # Prefer row-specific slot size when the affine model can identify S1..S4.
            ew, eh, ea = med_w, med_h, med_area
            if pattern == "slot" and model.get("used", False):
                rid = _slot_row_id_from_model_v15(model, expected)
                rowdims = _slot_row_dims_v14(model, selected, med_area)
                if rid in rowdims:
                    ew, eh, ea = rowdims[rid]

            d = _select_best_near_center_v15(
                pool, used_ids, (x, y), 1.20 * tol, pattern, ea, slot_area_limits
            )
            if d is None:
                d = recover_local_object_v15_center(
                    gray, (x, y), ew, eh, pattern, nearest_pitch, ea, slot_area_limits
                )
                if d is not None:
                    pool.append(d)
            if d is None or _desc_close_to_any_selected_v12(d, selected):
                continue
            if _source_of(d) != "center_recursive_recovered":
                _mark_source(d, "center_recursive_selected")
            selected.append(d)
            used_ids.add(id(d))
            queue.append(d)  # CRITICAL: newly found centers continue searching outward.
            recovered += 1

    meta.update({
        "center_recursive_attempts": int(attempts),
        "center_recursive_recovered": int(recovered),
        "center_recursive_boundary_hits": int(boundary_hits),
    })
    return selected, meta



# ============================================================
# V16: FIXED CENTER-PITCH ARRAY SEARCH
# ============================================================
def _image_allows_diagonal_v16(image_name: str) -> bool:
    """Per user rule, only filenames whose stem contains '-1' use 45-degree search."""
    return V16_DIAGONAL_FILENAME_TOKEN in Path(image_name).stem


def _robust_one_pitch_v16(values: List[float]) -> float:
    """Return one robust fundamental spacing; never subdivide/re-estimate later."""
    arr = np.asarray([v for v in values if np.isfinite(v) and v >= V16_PITCH_MIN_PX], dtype=float)
    if len(arr) == 0:
        return np.nan
    # Use nearest-neighbor values only (caller supplies these), then choose the densest
    # relative cluster. This avoids a few 2x gaps changing the fixed pitch.
    arr.sort()
    clusters = []
    for v in arr:
        placed = False
        for c in clusters:
            med = float(np.median(c))
            if abs(v - med) <= max(2.0, 0.16 * med):
                c.append(float(v)); placed = True; break
        if not placed:
            clusters.append([float(v)])
    clusters.sort(key=lambda c: (-len(c), float(np.median(c))))
    return float(np.median(clusters[0]))


def _fixed_axis_pitch_from_centers_v16(
    selected: List[ShapeDesc], axis: str, image_shape: Tuple[int, int]
) -> float:
    """Estimate exactly one H or V center pitch from initial accepted centers."""
    if len(selected) < 2:
        return np.nan
    P = np.asarray([d.center for d in selected], dtype=float)
    tan_tol = math.tan(math.radians(V16_DIRECTION_ANGLE_TOL_DEG))
    nearest = []
    max_allowed = (image_shape[1] if axis == 'x' else image_shape[0]) * V16_PITCH_MAX_IMAGE_FRAC
    for i in range(len(P)):
        vals = []
        for j in range(len(P)):
            if i == j:
                continue
            dx = abs(float(P[j, 0] - P[i, 0]))
            dy = abs(float(P[j, 1] - P[i, 1]))
            if axis == 'x':
                if dx >= V16_PITCH_MIN_PX and dx <= max_allowed and dy <= tan_tol * max(dx, 1e-9):
                    vals.append(dx)
            else:
                if dy >= V16_PITCH_MIN_PX and dy <= max_allowed and dx <= tan_tol * max(dy, 1e-9):
                    vals.append(dy)
        if vals:
            nearest.append(min(vals))
    return _robust_one_pitch_v16(nearest)


def _nearest_center_distance_v16(selected: List[ShapeDesc], image_shape: Tuple[int, int]) -> float:
    if len(selected) < 2:
        return np.nan
    P = np.asarray([d.center for d in selected], dtype=float)
    nearest = []
    max_allowed = max(image_shape) * V16_PITCH_MAX_IMAGE_FRAC
    for i in range(len(P)):
        vals = []
        for j in range(len(P)):
            if i == j:
                continue
            d = float(np.linalg.norm(P[j] - P[i]))
            if V16_PITCH_MIN_PX <= d <= max_allowed:
                vals.append(d)
        if vals:
            nearest.append(min(vals))
    return _robust_one_pitch_v16(nearest)


def estimate_fixed_pitches_v16(
    initial_selected: List[ShapeDesc], image_shape: Tuple[int, int]
) -> Tuple[float, float]:
    """Estimate H/V pitch once. Missing axis uses the other axis or nearest-center fallback."""
    hp = _fixed_axis_pitch_from_centers_v16(initial_selected, 'x', image_shape)
    vp = _fixed_axis_pitch_from_centers_v16(initial_selected, 'y', image_shape)
    fallback = _nearest_center_distance_v16(initial_selected, image_shape)
    if not np.isfinite(hp):
        hp = vp if np.isfinite(vp) else fallback
    if not np.isfinite(vp):
        vp = hp if np.isfinite(hp) else fallback
    return float(hp), float(vp)


def _complete_regular_shape_ok_v16(
    d: ShapeDesc,
    pattern: str,
    expected_area: float,
    slot_area_limits: Optional[Tuple[float, float]],
) -> bool:
    """Blur is allowed, deformation/incomplete boundaries are not."""
    pa = contour_pixel_count(d.contour)
    if pa < MIN_COMPONENT_AREA_PX:
        return False
    # Via objects are nominally same-size, so expected-area consistency is useful.
    # Short-slot rows intentionally have four different Y lengths, so do NOT impose
    # this narrower expected-area guard on slot; its existing 0.1x~10x image rule remains.
    if pattern == 'via' and np.isfinite(expected_area) and expected_area > 0:
        if pa < V16_RECOVER_MIN_EXPECTED_AREA_FACTOR * expected_area:
            return False
        if pa > V16_RECOVER_MAX_EXPECTED_AREA_FACTOR * expected_area:
            return False
    if pattern == 'via':
        return (
            d.circularity >= VIA_V16_RECOVER_MIN_CIRCULARITY
            and d.aspect <= VIA_V16_RECOVER_MAX_AXIS_RATIO
            and d.solidity >= VIA_V16_RECOVER_MIN_SOLIDITY
        )
    if slot_area_limits is not None:
        lo, hi = slot_area_limits
        if pa < max(SLOT_MIN_AREA_PX, lo) or pa > hi:
            return False
    reg = slot_regularity_metrics(d)
    return (
        reg['regularity_solidity'] >= SLOT_V16_RECOVER_MIN_SOLIDITY
        and reg['regularity_convexity'] >= SLOT_V16_RECOVER_MIN_CONVEXITY
        and reg['regularity_template_iou'] >= SLOT_V16_RECOVER_MIN_TEMPLATE_IOU
    )


def _select_pool_candidate_fixed_site_v16(
    pool: List[ShapeDesc],
    selected: List[ShapeDesc],
    expected_xy: Tuple[float, float],
    tol: float,
    pattern: str,
    expected_area: float,
    slot_area_limits: Optional[Tuple[float, float]],
) -> Optional[ShapeDesc]:
    ex = np.asarray(expected_xy, dtype=float)
    cand = []
    for d in pool:
        if _desc_close_to_any_selected_v12(d, selected):
            continue
        dist = float(np.linalg.norm(d.center - ex))
        if dist > tol:
            continue
        if not _complete_regular_shape_ok_v16(d, pattern, expected_area, slot_area_limits):
            continue
        # Center agreement dominates; then prefer cleaner/primary candidates.
        score = dist / max(tol, 1e-9)
        if pattern == 'via':
            score += 0.30 * (1.0 - min(1.0, d.circularity)) + 0.20 * (1.0 - min(1.0, d.solidity))
        else:
            reg = slot_regularity_metrics(d)
            score += 0.18 * (1.0 - min(1.0, reg['regularity_template_iou']))
            score += 0.12 * (1.0 - min(1.0, reg['regularity_solidity']))
        if _source_of(d) == 'primary':
            score -= 0.08
        cand.append((score, d))
    if not cand:
        return None
    cand.sort(key=lambda z: z[0])
    return cand[0][1]


def recover_local_object_v16_fixed(
    gray: np.ndarray,
    expected_xy: Tuple[float, float],
    expected_w: float,
    expected_h: float,
    pattern: str,
    fixed_pitch: float,
    expected_area: float,
    slot_area_limits: Optional[Tuple[float, float]],
) -> Optional[ShapeDesc]:
    """Local blur-tolerant search at ONE fixed predicted center; shape QC stays strict."""
    H, W = gray.shape
    cx, cy = expected_xy
    # Give blurred slot/via enough context, but do not let ROI reach the next lattice site.
    half_w = max(12.0, 0.85 * ARRAY_LOCAL_ROI_SIZE_FACTOR * max(5.0, expected_w) / 2.0)
    half_h = max(12.0, 0.85 * ARRAY_LOCAL_ROI_SIZE_FACTOR * max(5.0, expected_h) / 2.0)
    if np.isfinite(fixed_pitch) and fixed_pitch > 0:
        cap = 0.46 * fixed_pitch
        half_w = min(max(half_w, expected_w), cap)
        half_h = min(max(half_h, expected_h), cap)
    x0 = max(0, int(math.floor(cx - half_w)))
    x1 = min(W, int(math.ceil(cx + half_w + 1)))
    y0 = max(0, int(math.floor(cy - half_h)))
    y1 = min(H, int(math.ceil(cy + half_h + 1)))
    if x1 - x0 < 9 or y1 - y0 < 9:
        return None
    if _local_contrast_at_site_v14(gray, expected_xy, expected_w, expected_h) < V16_LOCAL_CONTRAST_GRAY:
        return None

    roi = gray[y0:y1, x0:x1]
    masks: Dict[str, np.ndarray] = {'v16_primary': segment_features(roi)}
    # A few enhanced masks help blur, but acceptance is governed by strict geometry below.
    try:
        em = enhanced_masks_v12(roi)
        for k in ('enh_clahe_otsu', 'enh_adaptive', 'enh_local_12'):
            if k in em:
                masks[k] = em[k]
    except Exception:
        pass
    percs = V16_LOCAL_PERCENTILES_DARK if FEATURE_POLARITY.lower() == 'dark' else V16_LOCAL_PERCENTILES_BRIGHT
    for q in percs:
        thr = float(np.percentile(roi, q))
        if FEATURE_POLARITY.lower() == 'dark':
            mm = (roi.astype(np.float32) <= thr).astype(np.uint8) * 255
        else:
            mm = (roi.astype(np.float32) >= thr).astype(np.uint8) * 255
        masks[f'v16_pct_{q:g}'] = _postprocess_binary_v12(mm, 1)

    center_tol = V16_LOCAL_CENTER_TOL_FRAC * max(half_w, half_h)
    candidates = []
    m = V16_LOCAL_ROI_BORDER_MARGIN_PX
    for mm in masks.values():
        contours, _ = cv2.findContours(mm, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        for c in contours:
            if contour_pixel_count(c) < MIN_COMPONENT_AREA_PX:
                continue
            bx, by, bw, bh = cv2.boundingRect(c)
            # Critical completeness test: an accepted recovered structure must have a
            # closed visible boundary fully inside the local crop, not be clipped.
            if bx <= m or by <= m or bx + bw >= roi.shape[1] - m or by + bh >= roi.shape[0] - m:
                continue
            cg = c.astype(np.int32).copy()
            cg[:, 0, 0] += x0
            cg[:, 0, 1] += y0
            d = describe_contour(cg)
            if d is None:
                continue
            if float(np.linalg.norm(d.center - np.asarray(expected_xy, dtype=float))) > center_tol:
                continue
            if not _complete_regular_shape_ok_v16(d, pattern, expected_area, slot_area_limits):
                continue
            _mark_source(d, 'fixed_center_recovered')
            candidates.append(d)
    if not candidates:
        return None
    return min(candidates, key=lambda d: _local_candidate_score_v12(d, expected_xy, expected_area, pattern))



def _largest_fixed_grid_component_v16(
    items: List[ShapeDesc],
    hp: float,
    vp: float,
    allow_diag: bool,
) -> List[ShapeDesc]:
    """Keep the largest connected set whose CENTER differences match the one frozen grid pitch."""
    if len(items) < 3:
        return list(items)
    orth = float(np.median([hp, vp]))
    diag = orth / V16_DIAGONAL_DIVISOR if allow_diag else np.nan
    targets = [
        np.array([ hp, 0.0]), np.array([-hp, 0.0]),
        np.array([0.0,  vp]), np.array([0.0, -vp]),
    ]
    if allow_diag and np.isfinite(diag):
        targets.extend([
            np.array([ diag,  diag]), np.array([ diag, -diag]),
            np.array([-diag,  diag]), np.array([-diag, -diag]),
        ])
    tol = max(4.0, V16_SITE_TOL_FRAC * min(hp, vp))
    P = np.asarray([d.center for d in items], dtype=float)
    adj = [[] for _ in items]
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            dv = P[j] - P[i]
            if min(float(np.linalg.norm(dv - t)) for t in targets) <= tol:
                adj[i].append(j); adj[j].append(i)
    seen = set(); comps = []
    for i in range(len(items)):
        if i in seen:
            continue
        stack = [i]; seen.add(i); comp = []
        while stack:
            k = stack.pop(); comp.append(k)
            for q in adj[k]:
                if q not in seen:
                    seen.add(q); stack.append(q)
        comps.append(comp)
    comps.sort(key=lambda c: (-len(c), min(c)))
    best = comps[0] if comps else []
    # Do not throw everything away if the pitch estimate was under-supported.
    if len(best) < 2:
        return list(items)
    return [items[i] for i in best]


def recursive_fixed_center_recovery_v16(
    gray: np.ndarray,
    pool: List[ShapeDesc],
    initial_selected: List[ShapeDesc],
    pattern: str,
    image_name: str,
    slot_area_limits: Optional[Tuple[float, float]] = None,
) -> Tuple[List[ShapeDesc], Dict[str, object]]:
    """
    Fixed recursive center search.

    The H/V pitch is estimated ONCE from initial accepted centers. New structures use
    exactly the same vectors; recovered centers never cause a smaller/re-fitted pitch.
    '-1' filenames additionally use four 45-degree vectors whose x/y component is
    fixed_orthogonal_pitch / 1.414.
    """
    allow_diag = _image_allows_diagonal_v16(image_name)
    meta = {
        'array_assist_used': False,
        'array_expected_positions': 0,
        'array_recovered_count': 0,
        'array_x_pitch_px': np.nan,
        'array_y_pitch_px': np.nan,
        'array_occupancy_seed': np.nan,
        'array_initial_recovered_count': 0,
        'array_final_recovered_count': 0,
        'array_final_predicted_positions': 0,
        'array_recovery_attempts': 0,
        'array_final_recovery_attempts': 0,
        'array_recovery_skipped_due_cap': 0,
        'array_lattice_type': 'fixed_HV_45' if allow_diag else 'fixed_HV',
        'array_basis_a_dx_px': np.nan, 'array_basis_a_dy_px': np.nan,
        'array_basis_b_dx_px': np.nan, 'array_basis_b_dy_px': np.nan,
        'array_basis_a_length_px': np.nan, 'array_basis_b_length_px': np.nan,
        'array_basis_angle_deg': 45.0 if allow_diag else 90.0,
        'center_recursive_used': False,
        'center_horizontal_pitch_px': np.nan,
        'center_vertical_pitch_px': np.nan,
        'center_diagonal_component_px': np.nan,
        'center_recursive_attempts': 0,
        'center_recursive_recovered': 0,
        'center_recursive_boundary_hits': 0,
        'fixed_diagonal_enabled': bool(allow_diag),
        'fixed_pitch_frozen': True,
    }
    selected = list(initial_selected)
    if len(selected) < 2:
        return selected, meta

    # IMPORTANT: estimate once from INITIAL accepted centers only.
    hp, vp = estimate_fixed_pitches_v16(initial_selected, gray.shape)
    if not np.isfinite(hp) and not np.isfinite(vp):
        return selected, meta
    if not np.isfinite(hp): hp = vp
    if not np.isfinite(vp): vp = hp
    if hp < V16_PITCH_MIN_PX or vp < V16_PITCH_MIN_PX:
        return selected, meta

    # Use the frozen pitch to remove isolated shape-like decoys. This is NOT a second
    # pitch fit: the same hp/vp estimated above are kept unchanged for all later search.
    grid_seed = _largest_fixed_grid_component_v16(initial_selected, hp, vp, allow_diag)
    selected = list(grid_seed)
    if len(selected) < 2:
        return selected, meta

    orth_pitch = float(np.median([hp, vp]))
    diag_comp = orth_pitch / V16_DIAGONAL_DIVISOR if allow_diag else np.nan
    meta.update({
        'array_assist_used': True,
        'array_x_pitch_px': float(hp),
        'array_y_pitch_px': float(vp),
        'center_recursive_used': True,
        'center_horizontal_pitch_px': float(hp),
        'center_vertical_pitch_px': float(vp),
        'center_diagonal_component_px': float(diag_comp) if allow_diag else np.nan,
        'array_occupancy_seed': float(len(selected) / max(1, len(initial_selected))),
    })

    H, W = gray.shape
    widths = [max(1.0, _contour_xy_bounds(d)[1] - _contour_xy_bounds(d)[0]) for d in selected]
    heights = [max(1.0, _contour_xy_bounds(d)[3] - _contour_xy_bounds(d)[2]) for d in selected]
    areas = [contour_pixel_count(d.contour) for d in selected]
    med_w = float(np.median(widths))
    med_h = float(np.median(heights))
    # Use p85 height as ROI support for four different slot lengths, but size QC remains strict.
    search_h = float(np.percentile(heights, 85)) if pattern == 'slot' and len(heights) >= 3 else med_h
    med_area = float(np.median(areas))
    fixed_pitch = float(min(hp, vp))
    tol = max(4.0, V16_SITE_TOL_FRAC * fixed_pitch)
    occ_tol = max(4.0, V16_OCCUPIED_TOL_FRAC * fixed_pitch)
    visit_cell = max(3.0, 0.55 * tol)

    directions = [
        np.array([ hp, 0.0]), np.array([-hp, 0.0]),
        np.array([0.0,  vp]), np.array([0.0, -vp]),
    ]
    if allow_diag and np.isfinite(diag_comp) and diag_comp >= 2.0:
        directions.extend([
            np.array([ diag_comp,  diag_comp]), np.array([ diag_comp, -diag_comp]),
            np.array([-diag_comp,  diag_comp]), np.array([-diag_comp, -diag_comp]),
        ])

    from collections import deque
    queue = deque(selected)
    visited = set()
    attempts = recovered = boundary_hits = 0

    def occupied(xy: np.ndarray) -> bool:
        return any(float(np.linalg.norm(d.center - xy)) <= occ_tol for d in selected)

    while queue and attempts < V16_RECURSIVE_MAX_ATTEMPTS and recovered < V16_RECURSIVE_MAX_RECOVERED:
        base = queue.popleft()
        for dv in directions:
            expected = np.asarray(base.center, dtype=float) + dv
            x, y = float(expected[0]), float(expected[1])
            margin_x = max(V16_BOUNDARY_MARGIN_PX, 0.55 * med_w)
            margin_y = max(V16_BOUNDARY_MARGIN_PX, 0.55 * search_h)
            if x <= margin_x or x >= W - 1 - margin_x or y <= margin_y or y >= H - 1 - margin_y:
                boundary_hits += 1
                continue
            key = (int(round(x / visit_cell)), int(round(y / visit_cell)))
            if key in visited:
                continue
            visited.add(key)
            if occupied(expected):
                continue
            attempts += 1
            if attempts == 1 or attempts % V16_RECURSIVE_PROGRESS_EVERY == 0:
                print(
                    f"      [{image_name}] V16 fixed-center search {attempts}/{V16_RECURSIVE_MAX_ATTEMPTS}; "
                    f"recovered={recovered}; diagonal={'ON' if allow_diag else 'OFF'} ...",
                    flush=True,
                )

            # Keep one frozen pitch, but allow dimensions to follow the actual source row.
            bx0, bx1, by0, by1 = _contour_xy_bounds(base)
            bw = max(med_w, float(bx1 - bx0))
            bh = max(search_h, float(by1 - by0)) if pattern == 'slot' else max(med_h, float(by1 - by0))
            expected_area = med_area

            d = _select_pool_candidate_fixed_site_v16(
                pool, selected, (x, y), 1.15 * tol, pattern, expected_area, slot_area_limits
            )
            if d is None:
                d = recover_local_object_v16_fixed(
                    gray, (x, y), bw, bh, pattern, fixed_pitch, expected_area, slot_area_limits
                )
                if d is not None:
                    pool.append(d)
            if d is None or _desc_close_to_any_selected_v12(d, selected):
                continue
            if _source_of(d) != 'fixed_center_recovered':
                _mark_source(d, 'fixed_center_selected')
            selected.append(d)
            queue.append(d)  # new actual center continues with the SAME frozen vectors
            recovered += 1

    meta.update({
        'array_recovered_count': int(recovered),
        'array_final_recovered_count': int(recovered),
        'array_final_predicted_positions': int(len(visited)),
        'array_recovery_attempts': int(attempts),
        'array_final_recovery_attempts': int(attempts),
        'center_recursive_attempts': int(attempts),
        'center_recursive_recovered': int(recovered),
        'center_recursive_boundary_hits': int(boundary_hits),
    })
    return selected, meta


def _dedup_desc_list_v16(items: List[ShapeDesc]) -> List[ShapeDesc]:
    out = []
    for d in sorted(items, key=lambda z: (z.center[1], z.center[0])):
        if not _desc_close_to_any_selected_v12(d, out):
            out.append(d)
    return out


def _seed_area_consistency_v16(items: List[ShapeDesc], lo_factor: float = 0.22, hi_factor: float = 4.5) -> List[ShapeDesc]:
    if len(items) < 3:
        return list(items)
    aa = np.asarray([contour_pixel_count(d.contour) for d in items], dtype=float)
    med = float(np.median(aa))
    keep = [d for d, a in zip(items, aa) if lo_factor * med <= a <= hi_factor * med]
    return keep if len(keep) >= 2 else list(items)


def _assign_slot_rows_v16(selected: List[ShapeDesc]) -> List[int]:
    """Four different slot lengths are a stable identity even when '-1' columns are diagonally offset."""
    if len(selected) < 4:
        return []
    lengths = np.asarray([max(1.0, _contour_xy_bounds(d)[3] - _contour_xy_bounds(d)[2]) for d in selected], dtype=np.float32)
    # If lengths are sufficiently distinct, cluster by Y-length; then name clusters 1..4
    # by their median physical Y center (top to bottom), preserving the user's numbering.
    data = lengths.reshape(-1, 1)
    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 100, 0.01)
    try:
        _, labs, _ = cv2.kmeans(data, 4, None, criteria, 30, cv2.KMEANS_PP_CENTERS)
        labs = labs.ravel().astype(int)
        cluster_median_y = []
        for k in range(4):
            ys = [selected[i].center[1] for i in range(len(selected)) if labs[i] == k]
            if not ys:
                raise ValueError('empty length cluster')
            cluster_median_y.append((float(np.median(ys)), k))
        order = {k: rid for rid, (_, k) in enumerate(sorted(cluster_median_y), 1)}
        return [int(order[int(k)]) for k in labs]
    except Exception:
        labs = cluster_1d_into_four_rows(np.asarray([d.center[1] for d in selected], dtype=np.float32))
        return [int(x) for x in labs]


def measure_via_v16(
    gray: np.ndarray, pool: List[ShapeDesc], px_nm: float, image_name: str, condition: str
) -> Tuple[List[dict], Optional[ShapeDesc], List[ShapeDesc], List[dict], Dict[str, object]]:
    strict = _dedup_desc_list_v16([d for d in pool if _via_strict_metrics_v12(d)[0]])
    strict = _seed_area_consistency_v16(strict)
    if len(strict) < 2:
        # Fallback seeds remain shape-complete; they are not the old ultra-loose V15 candidates.
        fallback = [d for d in pool if _complete_regular_shape_ok_v16(d, 'via', np.nan, None)]
        strict = _dedup_desc_list_v16(_seed_area_consistency_v16(fallback))
    selected, meta = recursive_fixed_center_recovery_v16(gray, pool, strict, 'via', image_name, None)

    rows = []
    for oid, d in enumerate(sorted(selected, key=lambda z: (z.center[1], z.center[0])), 1):
        pts = d.contour.reshape(-1, 2).astype(np.float64)
        w = float(pts[:, 0].max() - pts[:, 0].min())
        h = float(pts[:, 1].max() - pts[:, 1].min())
        rows.append({
            'condition': condition, 'image': image_name, 'pattern': 'via', 'object_id': oid,
            'width_nm': w * px_nm, 'height_nm': h * px_nm,
            'height_width_ratio': h / max(w, 1e-12), 'circularity': d.circularity,
            'solidity': d.solidity, 'axis_ratio': d.aspect, 'selection_source': _source_of(d),
            'pixel_size_nm': px_nm,
        })

    rejected = []
    for d in pool:
        if _desc_close_to_any_selected_v12(d, selected):
            continue
        ok, codes, reasons = _via_strict_metrics_v12(d)
        if _source_of(d) != 'primary' and not _complete_regular_shape_ok_v16(d, 'via', np.nan, None):
            continue
        if ok and not codes:
            codes = ['via_not_on_fixed_center_array']; reasons = ['regular via candidate was not selected by the frozen center-pitch array']
        elif not codes:
            codes = ['via_not_selected']; reasons = ['candidate was not selected']
        rejected.append(_make_reject_from_contour(
            d.contour, 'via_v16_filter', '+'.join(codes), '; '.join(reasons), desc=d,
            source=_source_of(d), circularity=d.circularity, solidity=d.solidity, axis_ratio=d.aspect,
        ))
    representative = max(selected, key=lambda d: (d.circularity, d.solidity)) if selected else None
    return rows, representative, selected, rejected, meta


def measure_slot_v16(
    gray: np.ndarray, pool: List[ShapeDesc], px_nm: float, image_name: str, condition: str
) -> Tuple[List[dict], List[ShapeDesc], Dict[int, ShapeDesc], List[dict], Dict[str, object]]:
    strict, rejects, reg_by_id, mean_area, med_area, _ = _slot_initial_filter_v12(pool)
    strict = _dedup_desc_list_v16(strict)
    if np.isfinite(mean_area) and mean_area > 0:
        area_lo = max(SLOT_MIN_AREA_PX, mean_area * SLOT_AREA_MIN_MEAN_FACTOR)
        area_hi = mean_area * SLOT_AREA_MAX_MEAN_FACTOR
    else:
        aa = np.asarray([contour_pixel_count(d.contour) for d in strict], dtype=float)
        mean_area = float(np.mean(aa)) if len(aa) else np.nan
        med_area = float(np.median(aa)) if len(aa) else np.nan
        area_lo = max(SLOT_MIN_AREA_PX, mean_area * SLOT_AREA_MIN_MEAN_FACTOR) if np.isfinite(mean_area) else SLOT_MIN_AREA_PX
        area_hi = mean_area * SLOT_AREA_MAX_MEAN_FACTOR if np.isfinite(mean_area) else np.inf

    if len(strict) < 2:
        fallback = []
        for d in pool:
            pa = contour_pixel_count(d.contour)
            if pa < area_lo or pa > area_hi:
                continue
            if _complete_regular_shape_ok_v16(d, 'slot', mean_area, (area_lo, area_hi)):
                fallback.append(d)
        strict = _dedup_desc_list_v16(fallback)

    selected, meta = recursive_fixed_center_recovery_v16(
        gray, pool, strict, 'slot', image_name, (area_lo, area_hi)
    )
    row_ids = _assign_slot_rows_v16(selected)
    if len(selected) < SLOT_REQUIRE_AT_LEAST_N or len(row_ids) != len(selected):
        for d in selected:
            rejects.append(_make_reject_from_contour(
                d.contour, 'slot_group_filter', 'slot_too_few_valid_objects',
                f'only {len(selected)} selected short-slot(s); need >= {SLOT_REQUIRE_AT_LEAST_N}',
                desc=d, source=_source_of(d),
            ))
        selected, row_ids = [], []

    rows, reps = [], {}
    if selected:
        order = sorted(range(len(selected)), key=lambda k: (row_ids[k], selected[k].center[0], selected[k].center[1]))
        selected = [selected[k] for k in order]
        row_ids = [row_ids[k] for k in order]
        for oid, (d, rid) in enumerate(zip(selected, row_ids), 1):
            pts = d.contour.reshape(-1, 2).astype(np.float64)
            width_px = float(pts[:, 0].max() - pts[:, 0].min())
            length_px = float(pts[:, 1].max() - pts[:, 1].min())
            reg = reg_by_id.get(id(d), slot_regularity_metrics(d))
            pa = contour_pixel_count(d.contour)
            rows.append({
                'condition': condition, 'image': image_name, 'pattern': 'slot', 'object_id': oid,
                'slot_row': int(rid), 'slot_width_nm': width_px * px_nm, 'slot_length_nm': length_px * px_nm,
                'slot_aspect_ratio': length_px / max(width_px, 1e-12),
                'slot_angle_from_y_deg': abs(90.0 - d.angle_from_x_deg),
                'slot_pixel_area_px': pa, 'slot_mean_candidate_area_px': mean_area,
                'slot_reference_median_area_px': med_area,
                'slot_area_ratio_to_mean': pa / mean_area if np.isfinite(mean_area) and mean_area > 0 else np.nan,
                'slot_regularity_solidity': reg['regularity_solidity'],
                'slot_regularity_convexity': reg['regularity_convexity'],
                'slot_rect_iou': reg['regularity_rect_iou'], 'slot_ellipse_iou': reg['regularity_ellipse_iou'],
                'slot_best_template_iou': reg['regularity_template_iou'], 'selection_source': _source_of(d),
                'pixel_size_nm': px_nm,
            })
            if rid not in reps:
                reps[rid] = d

    rejects = [r for r in rejects if r.get('_desc') is None or not _desc_close_to_any_selected_v12(r['_desc'], selected)]
    rejected_ids = {id(r.get('_desc')) for r in rejects if r.get('_desc') is not None}
    for d in pool:
        if _desc_close_to_any_selected_v12(d, selected) or id(d) in rejected_ids:
            continue
        pa = contour_pixel_count(d.contour)
        reg = slot_regularity_metrics(d)
        if _source_of(d) != 'primary' and not _complete_regular_shape_ok_v16(d, 'slot', mean_area, (area_lo, area_hi)):
            continue
        codes, reasons = [], []
        if pa < SLOT_MIN_AREA_PX:
            codes.append('slot_area_too_small'); reasons.append(f'filled pixel area {pa} px < {SLOT_MIN_AREA_PX} px')
        if np.isfinite(mean_area) and mean_area > 0:
            if pa < area_lo:
                codes.append('slot_area_below_0.1x_mean'); reasons.append(f'area {pa}px < {area_lo:.1f}px')
            if pa > area_hi:
                codes.append('slot_area_above_10x_mean'); reasons.append(f'area {pa}px > {area_hi:.1f}px')
        if reg['regularity_solidity'] < SLOT_V16_RECOVER_MIN_SOLIDITY:
            codes.append('slot_low_regularity_solidity'); reasons.append(f"solidity {reg['regularity_solidity']:.3f} too low")
        if reg['regularity_convexity'] < SLOT_V16_RECOVER_MIN_CONVEXITY:
            codes.append('slot_low_convexity'); reasons.append(f"convexity {reg['regularity_convexity']:.3f} too low")
        if reg['regularity_template_iou'] < SLOT_V16_RECOVER_MIN_TEMPLATE_IOU:
            codes.append('slot_irregular_template_match'); reasons.append(f"template IoU {reg['regularity_template_iou']:.3f} too low")
        if not codes:
            codes = ['slot_not_on_fixed_center_array']; reasons = ['regular short-slot candidate was not selected by the frozen center-pitch array']
        rejects.append(_make_reject_from_contour(
            d.contour, 'slot_v16_filter', '+'.join(codes), '; '.join(reasons), desc=d,
            source=_source_of(d), pixel_area_px=pa, **reg,
        ))
    return rows, selected, reps, rejects, meta



# ============================================================
# V17 SLOT: FOUR-ROW TOLERANT ARRAY MODEL
# ============================================================
def _slot_shape_ok_v17(
    d: ShapeDesc,
    area_limits: Tuple[float, float],
    seed: bool = False,
) -> bool:
    """Regular complete short-slot shape; aspect ratio is deliberately irrelevant."""
    pa = contour_pixel_count(d.contour)
    lo, hi = area_limits
    if pa < max(SLOT_MIN_AREA_PX, lo) or pa > hi:
        return False
    reg = slot_regularity_metrics(d)
    if seed:
        return (
            reg['regularity_solidity'] >= V17_SLOT_SEED_MIN_SOLIDITY
            and reg['regularity_convexity'] >= V17_SLOT_SEED_MIN_CONVEXITY
            and reg['regularity_template_iou'] >= V17_SLOT_SEED_MIN_TEMPLATE_IOU
        )
    return (
        reg['regularity_solidity'] >= V17_SLOT_RECOVER_MIN_SOLIDITY
        and reg['regularity_convexity'] >= V17_SLOT_RECOVER_MIN_CONVEXITY
        and reg['regularity_template_iou'] >= V17_SLOT_RECOVER_MIN_TEMPLATE_IOU
    )


def _fit_four_row_centers_v17(
    candidates: List[ShapeDesc], image_h: int, image_name: str = ""
) -> Optional[Dict[str, object]]:
    """
    Fit exactly four approximately equally spaced Y row centers.

    Unlike k-means(k=4), this can infer one completely missing row because every
    observed Y may be assigned to one of four arithmetic-progression levels.
    """
    if len(candidates) < 2:
        return None
    ys = np.asarray([float(d.center[1]) for d in candidates], dtype=float)
    min_pitch = V17_SLOT_MIN_ROW_PITCH_PX
    max_pitch = max(min_pitch + 1.0, V17_SLOT_MAX_ROW_PITCH_FRAC * image_h)

    hypotheses: List[float] = []
    for i in range(len(ys)):
        for j in range(i + 1, len(ys)):
            dy = abs(float(ys[j] - ys[i]))
            if dy < min_pitch * 0.65:
                continue
            for mult in (1, 2, 3):
                p = dy / mult
                if min_pitch <= p <= max_pitch:
                    hypotheses.append(float(p))
    # Fallback hypotheses help when only one pair of rows survives segmentation.
    hypotheses.extend([image_h / 5.0, image_h / 6.0, image_h / 7.0])
    hypotheses = [p for p in hypotheses if min_pitch <= p <= max_pitch]
    if not hypotheses:
        return None
    hypotheses = sorted(set(round(p, 3) for p in hypotheses))

    # V32: keep every original hypothesis and the original score/tie ordering.
    # Count support in four Y intervals by binary search first. Since the residual
    # penalty is nonnegative, 5*occupied_rows + support is an upper bound: models
    # below the best score cannot win and need no per-candidate distance matrix.
    # Batch both stages to avoid millions of tiny NumPy/Python calls. Memory is
    # bounded by batch size, rather than by all pitches * anchors * candidates.
    started = last_progress = time.perf_counter()
    if image_name:
        print(f"      [{image_name}] slot row fit: candidates={len(ys)}, pitches={len(hypotheses)}", flush=True)
    sorted_ys = np.sort(ys)
    levels = np.arange(4, dtype=float)
    anchors = np.repeat(ys, 4)
    anchor_levels = np.tile(levels, len(ys))
    best = None
    best_order = None
    for start in range(0, len(hypotheses), 32):
        pitches = np.asarray(hypotheses[start:start + 32], dtype=float)
        ps = np.repeat(pitches, len(anchors))
        origins = (anchors[None, :] - pitches[:, None] * anchor_levels[None, :]).ravel()
        centers = origins[:, None] + ps[:, None] * levels[None, :]
        tolerances = np.maximum(5.0, V17_SLOT_ROW_MODEL_RESIDUAL_FRAC * ps)
        inside = np.sum((centers >= -0.10 * image_h) & (centers <= 1.10 * image_h), axis=1)
        # Widen only this preliminary bound for roundoff and interval-boundary
        # ties. Exact acceptance still uses the original res <= tol below.
        pad = 8.0 * np.finfo(float).eps * (np.abs(centers) + tolerances[:, None] + 1.0)
        lo = np.searchsorted(sorted_ys, centers - tolerances[:, None] - pad, side='left')
        hi = np.searchsorted(sorted_ys, centers + tolerances[:, None] + pad, side='right')
        counts = hi - lo
        upper = 5 * np.count_nonzero(counts, axis=1) + np.sum(counts, axis=1)
        possible = np.flatnonzero((inside >= 3) & (np.sum(counts, axis=1) >= 2))
        # Strongest bounds first let later chunks skip losing models. Explicit
        # original order below preserves the old first-wins behavior on ties.
        possible = possible[np.argsort(-upper[possible], kind='stable')]
        for offset in range(0, len(possible), 128):
            ids = possible[offset:offset + 128]
            if best is not None:
                ids = ids[upper[ids] >= best[0]]
            if not len(ids):
                break
            if best is not None:
                # Dense, nearly perfect rows can give many models the same
                # support bound. To match the best score, at least half their
                # inliers must also fit within this maximum median residual.
                # Use a conservative lower bound on support and upper bound on
                # close points; skip this test when row intervals can overlap.
                max_median = (upper[ids] - best[0] + 1e-12) * ps[ids] / 1.6
                tight = (max_median < tolerances[ids]) & (2.0 * tolerances[ids] < ps[ids])
                positions = np.flatnonzero(tight)
                bounded_ids = ids[positions]
                if len(bounded_ids):
                    cc = centers[bounded_ids]
                    tt = tolerances[bounded_ids, None]
                    pp = pad[bounded_ids]
                    support_lo = np.sum(
                        np.searchsorted(sorted_ys, cc + tt - pp, side='right')
                        - np.searchsorted(sorted_ys, cc - tt + pp, side='left'), axis=1)
                    radius = max_median[positions, None]
                    close_hi = np.sum(
                        np.searchsorted(sorted_ys, cc + radius + pp, side='right')
                        - np.searchsorted(sorted_ys, cc - radius - pp, side='left'), axis=1)
                    keep = np.ones(len(ids), dtype=bool)
                    keep[positions] = close_hi >= (support_lo + 1) // 2
                    ids = ids[keep]
                if not len(ids):
                    continue
            dist = np.abs(ys[None, :, None] - centers[ids, None, :])
            ridx = np.argmin(dist, axis=2)
            res = np.take_along_axis(dist, ridx[:, :, None], axis=2)[:, :, 0]
            inl = res <= tolerances[ids, None]
            n_in = np.sum(inl, axis=1)
            valid = n_in >= 2
            ids, ridx, res, inl, n_in = ids[valid], ridx[valid], res[valid], inl[valid], n_in[valid]
            if not len(ids):
                continue
            occ = sum(np.any((ridx == r) & inl, axis=1).astype(int) for r in range(4))
            med_res = np.nanmedian(np.where(inl, res, np.nan), axis=1)
            scores = 5.0 * occ + 1.0 * n_in - 1.6 * med_res / np.maximum(ps[ids], 1e-9)
            winner = int(np.lexsort((-ids, -med_res, n_in, occ, scores))[-1])
            idx = int(ids[winner])
            order = start * len(anchors) + idx
            rec = (float(scores[winner]), int(occ[winner]), int(n_in[winner]),
                   -float(med_res[winner]), float(origins[idx]), float(ps[idx]),
                   ridx[winner].copy(), inl[winner].copy())
            if best is None or rec[:4] > best[:4] or (rec[:4] == best[:4] and order < best_order):
                best, best_order = rec, order
        now = time.perf_counter()
        if image_name and now - last_progress >= 5.0:
            print(f"      [{image_name}] slot row fit: pitches={min(start + 32, len(hypotheses))}/{len(hypotheses)}, elapsed={now - started:.1f}s", flush=True)
            last_progress = now
    if image_name:
        print(f"      [{image_name}] slot row fit done: elapsed={time.perf_counter() - started:.2f}s", flush=True)
    if best is None:
        return None

    _, occ, n_in, _, y0, p, ridx, inl = best
    # Refine with robust row medians. A missing row remains inferred from the fitted pitch.
    row_obs = []
    for r in range(4):
        vals = ys[(ridx == r) & inl]
        if len(vals):
            row_obs.append((r, float(np.median(vals))))
    if len(row_obs) >= 2:
        rr = np.asarray([a for a, _ in row_obs], dtype=float)
        yy = np.asarray([b for _, b in row_obs], dtype=float)
        A = np.column_stack([np.ones_like(rr), rr])
        coef, *_ = np.linalg.lstsq(A, yy, rcond=None)
        y0_ref = float(coef[0]); p_ref = float(coef[1])
        if min_pitch <= p_ref <= max_pitch:
            y0, p = y0_ref, p_ref
    centers = y0 + p * np.arange(4, dtype=float)
    tol = max(5.0, V17_SLOT_ROW_CENTER_TOL_FRAC * p)
    return {
        'used': True,
        'row_centers': centers,
        'row_pitch': float(p),
        'row_tol': float(tol),
        'occupied_seed_rows': int(occ),
        'seed_inliers': int(n_in),
    }


def _assign_to_row_centers_v17(
    candidates: List[ShapeDesc], row_centers: np.ndarray, row_tol: float
) -> Dict[int, List[ShapeDesc]]:
    out = {0: [], 1: [], 2: [], 3: []}
    for d in candidates:
        dy = np.abs(row_centers - float(d.center[1]))
        r = int(np.argmin(dy))
        if float(dy[r]) <= row_tol:
            out[r].append(d)
    return out


def _slot_horizontal_pitch_v17(row_map: Dict[int, List[ShapeDesc]], image_w: int) -> float:
    """One fixed horizontal center pitch for the whole image; never subdivided later."""
    nearest = []
    for r in range(4):
        xs = sorted(float(d.center[0]) for d in row_map.get(r, []))
        if len(xs) < 2:
            continue
        dif = np.diff(xs)
        dif = dif[(dif >= V16_PITCH_MIN_PX) & (dif <= V16_PITCH_MAX_IMAGE_FRAC * image_w)]
        nearest.extend(float(v) for v in dif)
    if not nearest:
        return np.nan

    # Cluster observed adjacent spacings. If missing objects created 2x gaps, prefer the
    # smallest well-supported cluster instead of inventing a smaller sub-pitch.
    arr = sorted(nearest)
    clusters: List[List[float]] = []
    for v in arr:
        placed = False
        for c in clusters:
            med = float(np.median(c))
            if abs(v - med) <= max(2.5, 0.18 * med):
                c.append(v); placed = True; break
        if not placed:
            clusters.append([v])
    max_support = max(len(c) for c in clusters)
    eligible = [c for c in clusters if len(c) >= max(2, int(math.ceil(0.35 * max_support)))]
    if not eligible:
        eligible = clusters
    eligible.sort(key=lambda c: (float(np.median(c)), -len(c)))
    return float(np.median(eligible[0]))


def _wrap_phase_residual_v17(x: float, phase: float, pitch: float) -> float:
    return float(((x - phase + 0.5 * pitch) % pitch) - 0.5 * pitch)


def _fit_phase_v17(values: List[float], pitch: float) -> float:
    if not values or not np.isfinite(pitch) or pitch <= 0:
        return np.nan
    vals = np.asarray(values, dtype=float)
    best_phase = float(vals[0] % pitch); best_score = float('inf')
    for anchor in vals:
        ph0 = float(anchor % pitch)
        res = np.asarray([_wrap_phase_residual_v17(v, ph0, pitch) for v in vals], dtype=float)
        ph = float((ph0 + np.median(res)) % pitch)
        rr = np.abs([_wrap_phase_residual_v17(v, ph, pitch) for v in vals])
        score = float(np.median(rr)) + 0.25 * float(np.mean(rr))
        if score < best_score:
            best_score, best_phase = score, ph
    return best_phase


def _fit_slot_phase_model_v17(
    row_map: Dict[int, List[ShapeDesc]], pitch: float, allow_diag: bool
) -> Tuple[float, float]:
    """
    Return (base_phase, per-row x shift).
    Normal files use one common phase. '-1' files may use 0 or +/- pitch/1.414
    per row, selecting whichever best matches already found centers.
    """
    obs = []
    for r in range(4):
        for d in row_map.get(r, []):
            obs.append((r, float(d.center[0])))
    if not obs:
        return np.nan, 0.0
    shifts = [0.0]
    if allow_diag:
        dc = pitch / V16_DIAGONAL_DIVISOR
        shifts.extend([dc, -dc])
    best = None
    for shift in shifts:
        adjusted = [x - r * shift for r, x in obs]
        phase = _fit_phase_v17(adjusted, pitch)
        if not np.isfinite(phase):
            continue
        residuals = [abs(_wrap_phase_residual_v17(x - r * shift, phase, pitch)) for r, x in obs]
        score = float(np.median(residuals)) + 0.25 * float(np.mean(residuals))
        rec = (score, abs(shift), phase, shift)
        if best is None or rec[:2] < best[:2]:
            best = rec
    return (float(best[2]), float(best[3])) if best is not None else (np.nan, 0.0)


def _expected_x_sites_v17(
    image_w: int, phase: float, row_shift: float, row_id0: int, pitch: float, margin: float
) -> List[float]:
    ph = float((phase + row_id0 * row_shift) % pitch)
    k0 = int(math.floor((margin - ph) / pitch)) - 1
    k1 = int(math.ceil((image_w - 1 - margin - ph) / pitch)) + 1
    xs = []
    for k in range(k0, k1 + 1):
        x = ph + k * pitch
        if margin <= x <= image_w - 1 - margin:
            xs.append(float(x))
    # Keep the visible row practical. The user expects roughly ten, not hundreds.
    if len(xs) > V17_SLOT_MAX_SITES_PER_ROW:
        mid = 0.5 * (image_w - 1)
        xs = sorted(xs, key=lambda x: abs(x - mid))[:V17_SLOT_MAX_SITES_PER_ROW]
        xs.sort()
    return xs


def _slot_candidate_site_score_v17(
    d: ShapeDesc, ex: float, ey: float, x_tol: float, y_tol: float
) -> float:
    dx = abs(float(d.center[0]) - ex) / max(x_tol, 1e-9)
    dy = abs(float(d.center[1]) - ey) / max(y_tol, 1e-9)
    reg = slot_regularity_metrics(d)
    return (
        0.55 * dx + 0.45 * dy
        + 0.15 * (1.0 - min(1.0, reg['regularity_template_iou']))
        + 0.10 * (1.0 - min(1.0, reg['regularity_solidity']))
    )


def _pick_pool_slot_at_site_v17(
    pool: List[ShapeDesc], used: List[ShapeDesc], ex: float, ey: float,
    x_tol: float, y_tol: float, area_limits: Tuple[float, float]
) -> Optional[ShapeDesc]:
    cand = []
    for d in pool:
        if _desc_close_to_any_selected_v12(d, used):
            continue
        if abs(float(d.center[0]) - ex) > x_tol or abs(float(d.center[1]) - ey) > y_tol:
            continue
        if not _slot_shape_ok_v17(d, area_limits, seed=False):
            continue
        cand.append((_slot_candidate_site_score_v17(d, ex, ey, x_tol, y_tol), d))
    if not cand:
        return None
    cand.sort(key=lambda z: z[0])
    return cand[0][1]


def _recover_local_slot_v17(
    gray: np.ndarray,
    expected_xy: Tuple[float, float],
    x_pitch: float,
    row_pitch: float,
    expected_w: float,
    expected_h: float,
    area_limits: Tuple[float, float],
    aggressive: bool = False,
) -> Optional[ShapeDesc]:
    """Low-contrast local search at a known 4-row lattice site; geometry QC stays strict."""
    H, W = gray.shape
    cx, cy = map(float, expected_xy)
    half_w = max(12.0, min(0.48 * x_pitch, max(1.8 * expected_w, 14.0)))
    # Give long short-slots enough Y context; center-row gating prevents neighboring rows.
    half_h = max(14.0, min(0.72 * row_pitch, max(0.95 * expected_h + 6.0, 18.0)))
    x0 = max(0, int(math.floor(cx - half_w))); x1 = min(W, int(math.ceil(cx + half_w + 1)))
    y0 = max(0, int(math.floor(cy - half_h))); y1 = min(H, int(math.ceil(cy + half_h + 1)))
    if x1 - x0 < 9 or y1 - y0 < 9:
        return None
    if not aggressive:
        if _local_contrast_at_site_v14(gray, expected_xy, expected_w, expected_h) < V17_SLOT_RECOVER_LOCAL_CONTRAST_GRAY:
            return None
    roi = gray[y0:y1, x0:x1]
    masks: Dict[str, np.ndarray] = {'v17_primary': segment_features(roi)}
    try:
        em = enhanced_masks_v12(roi)
        for k in ('enh_clahe_otsu', 'enh_adaptive', 'enh_local_12'):
            if k in em:
                masks[k] = em[k]
    except Exception:
        pass
    percs = V17_SLOT_LOCAL_PERCENTILES_DARK if FEATURE_POLARITY.lower() == 'dark' else V17_SLOT_LOCAL_PERCENTILES_BRIGHT
    for q in percs:
        thr = float(np.percentile(roi, q))
        if FEATURE_POLARITY.lower() == 'dark':
            mm = (roi.astype(np.float32) <= thr).astype(np.uint8) * 255
        else:
            mm = (roi.astype(np.float32) >= thr).astype(np.uint8) * 255
        masks[f'v17_pct_{q:g}'] = _postprocess_binary_v12(mm, 1)

    x_tol = max(6.0, V17_SLOT_RECOVER_CENTER_X_FRAC * x_pitch)
    y_tol = max(6.0, V17_SLOT_RECOVER_CENTER_Y_FRAC * row_pitch)
    m = V17_SLOT_RECOVER_ROI_BORDER_MARGIN_PX
    candidates = []
    for mm in masks.values():
        contours, _ = cv2.findContours(mm, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        for c in contours:
            if contour_pixel_count(c) < SLOT_MIN_AREA_PX:
                continue
            bx, by, bw, bh = cv2.boundingRect(c)
            # Do not accept clipped/incomplete boundaries, even in aggressive mode.
            if bx <= m or by <= m or bx + bw >= roi.shape[1] - m or by + bh >= roi.shape[0] - m:
                continue
            cg = c.astype(np.int32).copy()
            cg[:, 0, 0] += x0; cg[:, 0, 1] += y0
            d = describe_contour(cg)
            if d is None:
                continue
            if abs(float(d.center[0]) - cx) > x_tol or abs(float(d.center[1]) - cy) > y_tol:
                continue
            if not _slot_shape_ok_v17(d, area_limits, seed=False):
                continue
            _mark_source(d, 'slot_4row_recovered')
            score = _slot_candidate_site_score_v17(d, cx, cy, x_tol, y_tol)
            candidates.append((score, d))
    if not candidates:
        return None
    candidates.sort(key=lambda z: z[0])
    return candidates[0][1]


def _slot_four_row_array_v17(
    gray: np.ndarray,
    pool: List[ShapeDesc],
    strict_seed: List[ShapeDesc],
    image_name: str,
    area_limits: Tuple[float, float],
) -> Tuple[List[ShapeDesc], List[int], Dict[str, object]]:
    """Select/recover a 4-row, roughly ten-per-row short-slot array."""
    H, W = gray.shape
    allow_diag = _image_allows_diagonal_v16(image_name)
    # Add shape-complete relaxed candidates only for model fitting; weird contours stay out.
    relaxed = [d for d in pool if _slot_shape_ok_v17(d, area_limits, seed=True)]
    model_cands = _dedup_desc_list_v16(strict_seed + relaxed)
    row_model = _fit_four_row_centers_v17(model_cands, H)
    meta = {
        'array_assist_used': False,
        'array_lattice_type': 'slot_4row_H_45' if allow_diag else 'slot_4row_H',
        'array_x_pitch_px': np.nan, 'array_y_pitch_px': np.nan,
        'center_horizontal_pitch_px': np.nan, 'center_vertical_pitch_px': np.nan,
        'center_diagonal_component_px': np.nan,
        'fixed_diagonal_enabled': bool(allow_diag), 'fixed_pitch_frozen': True,
        'center_recursive_used': False, 'center_recursive_attempts': 0,
        'center_recursive_recovered': 0, 'center_recursive_boundary_hits': 0,
        'array_recovered_count': 0, 'array_final_recovered_count': 0,
        'array_expected_positions': 0, 'array_final_predicted_positions': 0,
        'array_recovery_attempts': 0, 'array_final_recovery_attempts': 0,
        'slot_expected_rows': 4,
        'slot_row_pitch_px': np.nan,
        'slot_horizontal_pitch_px': np.nan,
        'slot_phase_shift_per_row_px': np.nan,
        'slot_seed_rows_found': 0,
        'slot_row1_count': 0, 'slot_row2_count': 0, 'slot_row3_count': 0, 'slot_row4_count': 0,
        'slot_row1_expected_sites': 0, 'slot_row2_expected_sites': 0,
        'slot_row3_expected_sites': 0, 'slot_row4_expected_sites': 0,
    }
    if row_model is None:
        return [], [], meta
    row_centers = np.asarray(row_model['row_centers'], dtype=float)
    row_pitch = float(row_model['row_pitch'])
    row_tol = float(row_model['row_tol'])
    row_map = _assign_to_row_centers_v17(model_cands, row_centers, row_tol)
    hp = _slot_horizontal_pitch_v17(row_map, W)
    if not np.isfinite(hp) or hp < V16_PITCH_MIN_PX:
        # Fall back to the V16 H pitch estimator, but keep the four-row model.
        hp = _fixed_axis_pitch_from_centers_v16(model_cands, 'x', gray.shape)
    if not np.isfinite(hp) or hp < V16_PITCH_MIN_PX:
        return [], [], meta

    phase, row_shift = _fit_slot_phase_model_v17(row_map, hp, allow_diag)
    if not np.isfinite(phase):
        return [], [], meta

    # Typical dimensions. Row-specific height is used when available; missing rows get
    # a generous global reference so long/short members are not clipped by the search ROI.
    widths = [max(1.0, _contour_xy_bounds(d)[1] - _contour_xy_bounds(d)[0]) for d in model_cands]
    heights = [max(1.0, _contour_xy_bounds(d)[3] - _contour_xy_bounds(d)[2]) for d in model_cands]
    med_w = float(np.median(widths)) if widths else max(6.0, 0.20 * hp)
    global_h = float(np.percentile(heights, 85)) if heights else max(10.0, 0.45 * row_pitch)
    margin = max(4.0, 0.55 * med_w)

    expected_by_row: Dict[int, List[float]] = {}
    for r in range(4):
        xs = _expected_x_sites_v17(W, phase, row_shift, r, hp, margin)
        expected_by_row[r] = xs
        meta[f'slot_row{r+1}_expected_sites'] = len(xs)
    total_sites = sum(len(v) for v in expected_by_row.values())

    selected: List[ShapeDesc] = []
    row_ids: List[int] = []
    used_sites = set()
    x_tol = max(6.0, V17_SLOT_SITE_X_TOL_FRAC * hp)
    y_tol = max(6.0, V17_SLOT_SITE_Y_TOL_FRAC * row_pitch)

    # First use every already-segmented, shape-complete candidate whose CENTER is only
    # approximately near a predicted row/site. No strict neighbor connectivity required.
    for r in range(4):
        ey = float(row_centers[r])
        for si, ex in enumerate(expected_by_row[r]):
            d = _pick_pool_slot_at_site_v17(pool, selected, ex, ey, x_tol, y_tol, area_limits)
            if d is None:
                continue
            if _source_of(d) == 'primary':
                _mark_source(d, 'slot_4row_selected')
            selected.append(d); row_ids.append(r + 1); used_sites.add((r, si))

    recovered = attempts = 0
    # Search all missing theoretical sites once. This is what recovers rows that were
    # nearly absent from the original segmentation.
    for pass_no in (1, 2):
        counts_now = {r: sum(1 for rr in row_ids if rr == r + 1) for r in range(4)}
        if pass_no == 2 and not V17_SLOT_AGGRESSIVE_RECOVERY_FOR_MISSING_ROWS:
            break
        for r in range(4):
            # Second pass is reserved for genuinely weak rows.
            if pass_no == 2 and counts_now[r] >= V17_SLOT_ROW_WEAK_COUNT:
                continue
            row_heights = [max(1.0, _contour_xy_bounds(d)[3] - _contour_xy_bounds(d)[2])
                           for d, rr in zip(selected, row_ids) if rr == r + 1]
            expected_h = float(np.median(row_heights)) if row_heights else global_h
            ey = float(row_centers[r])
            for si, ex in enumerate(expected_by_row[r]):
                if (r, si) in used_sites:
                    continue
                attempts += 1
                if attempts == 1 or attempts % V17_SLOT_PROGRESS_EVERY == 0:
                    print(
                        f"      [{image_name}] V17 4-row short-slot recovery {attempts}; "
                        f"pass={pass_no}; row={r+1}; recovered={recovered} ...",
                        flush=True,
                    )
                d = _pick_pool_slot_at_site_v17(pool, selected, ex, ey, 1.10 * x_tol, 1.10 * y_tol, area_limits)
                if d is None:
                    d = _recover_local_slot_v17(
                        gray, (ex, ey), hp, row_pitch, med_w, expected_h,
                        area_limits, aggressive=(pass_no == 2)
                    )
                    if d is not None:
                        pool.append(d)
                if d is None or _desc_close_to_any_selected_v12(d, selected):
                    continue
                if _source_of(d) not in {'slot_4row_recovered'}:
                    _mark_source(d, 'slot_4row_selected')
                selected.append(d); row_ids.append(r + 1); used_sites.add((r, si)); recovered += 1

    # Stable order: row top->bottom, then x left->right.
    order = sorted(range(len(selected)), key=lambda k: (row_ids[k], selected[k].center[0], selected[k].center[1]))
    selected = [selected[k] for k in order]
    row_ids = [row_ids[k] for k in order]
    counts = [sum(1 for rr in row_ids if rr == r + 1) for r in range(4)]

    meta.update({
        'array_assist_used': True,
        'array_x_pitch_px': float(hp), 'array_y_pitch_px': float(row_pitch),
        'center_horizontal_pitch_px': float(hp), 'center_vertical_pitch_px': float(row_pitch),
        'center_diagonal_component_px': float(abs(row_shift)) if allow_diag else np.nan,
        'center_recursive_used': True,
        'center_recursive_attempts': int(attempts), 'center_recursive_recovered': int(recovered),
        'array_recovered_count': int(recovered), 'array_final_recovered_count': int(recovered),
        'array_expected_positions': int(total_sites), 'array_final_predicted_positions': int(total_sites),
        'array_recovery_attempts': int(attempts), 'array_final_recovery_attempts': int(attempts),
        'slot_row_pitch_px': float(row_pitch), 'slot_horizontal_pitch_px': float(hp),
        'slot_phase_shift_per_row_px': float(row_shift),
        'slot_seed_rows_found': int(row_model['occupied_seed_rows']),
        'slot_row1_count': int(counts[0]), 'slot_row2_count': int(counts[1]),
        'slot_row3_count': int(counts[2]), 'slot_row4_count': int(counts[3]),
        'array_occupancy_seed': float(len(model_cands) / max(1, total_sites)),
    })
    return selected, row_ids, meta


def measure_slot_v17(
    gray: np.ndarray, pool: List[ShapeDesc], px_nm: float, image_name: str, condition: str
) -> Tuple[List[dict], List[ShapeDesc], Dict[int, ShapeDesc], List[dict], Dict[str, object]]:
    """V17 short-slot measurement: four known rows + tolerant equal-X-pitch completion."""
    strict, rejects, reg_by_id, mean_area, med_area, _ = _slot_initial_filter_v12(pool)
    strict = _dedup_desc_list_v16(strict)
    if np.isfinite(mean_area) and mean_area > 0:
        area_lo = max(SLOT_MIN_AREA_PX, mean_area * SLOT_AREA_MIN_MEAN_FACTOR)
        area_hi = mean_area * SLOT_AREA_MAX_MEAN_FACTOR
    else:
        aa = np.asarray([contour_pixel_count(d.contour) for d in strict], dtype=float)
        mean_area = float(np.mean(aa)) if len(aa) else np.nan
        med_area = float(np.median(aa)) if len(aa) else np.nan
        area_lo = max(SLOT_MIN_AREA_PX, mean_area * SLOT_AREA_MIN_MEAN_FACTOR) if np.isfinite(mean_area) else SLOT_MIN_AREA_PX
        area_hi = mean_area * SLOT_AREA_MAX_MEAN_FACTOR if np.isfinite(mean_area) else np.inf
    area_limits = (float(area_lo), float(area_hi))

    # If strict QC found few seeds, bring in only complete regular relaxed shapes.
    if len(strict) < 4:
        extra = [d for d in pool if _slot_shape_ok_v17(d, area_limits, seed=True)]
        strict = _dedup_desc_list_v16(strict + extra)

    selected, row_ids, meta = _slot_four_row_array_v17(gray, pool, strict, image_name, area_limits)

    rows: List[dict] = []
    reps: Dict[int, ShapeDesc] = {}
    for oid, (d, rid) in enumerate(zip(selected, row_ids), 1):
        pts = d.contour.reshape(-1, 2).astype(np.float64)
        width_px = float(pts[:, 0].max() - pts[:, 0].min())
        length_px = float(pts[:, 1].max() - pts[:, 1].min())
        reg = reg_by_id.get(id(d), slot_regularity_metrics(d))
        pa = contour_pixel_count(d.contour)
        rows.append({
            'condition': condition, 'image': image_name, 'pattern': 'slot', 'object_id': oid,
            'slot_row': int(rid), 'slot_width_nm': width_px * px_nm, 'slot_length_nm': length_px * px_nm,
            'slot_aspect_ratio': length_px / max(width_px, 1e-12),
            'slot_angle_from_y_deg': abs(90.0 - d.angle_from_x_deg),
            'slot_pixel_area_px': pa, 'slot_mean_candidate_area_px': mean_area,
            'slot_reference_median_area_px': med_area,
            'slot_area_ratio_to_mean': pa / mean_area if np.isfinite(mean_area) and mean_area > 0 else np.nan,
            'slot_regularity_solidity': reg['regularity_solidity'],
            'slot_regularity_convexity': reg['regularity_convexity'],
            'slot_rect_iou': reg['regularity_rect_iou'], 'slot_ellipse_iou': reg['regularity_ellipse_iou'],
            'slot_best_template_iou': reg['regularity_template_iou'], 'selection_source': _source_of(d),
            'pixel_size_nm': px_nm,
        })
        if rid not in reps:
            reps[rid] = d

    # Remove old strict-filter rejects that actually got recovered/selected.
    rejects = [r for r in rejects if r.get('_desc') is None or not _desc_close_to_any_selected_v12(r['_desc'], selected)]
    rejected_ids = {id(r.get('_desc')) for r in rejects if r.get('_desc') is not None}

    # Report unused candidates without the old overly strict "not_on_fixed_center_array" reason.
    row_centers = None
    if meta.get('array_assist_used', False) and np.isfinite(meta.get('slot_row_pitch_px', np.nan)):
        # Reconstruct four centers from selected rows where possible for diagnostics only.
        rc = []
        for r in range(1, 5):
            yy = [float(d.center[1]) for d, rr in zip(selected, row_ids) if rr == r]
            rc.append(float(np.median(yy)) if yy else np.nan)
        row_centers = np.asarray(rc, dtype=float)

    for d in pool:
        if _desc_close_to_any_selected_v12(d, selected) or id(d) in rejected_ids:
            continue
        pa = contour_pixel_count(d.contour)
        reg = slot_regularity_metrics(d)
        # Ignore enhanced threshold fragments that are not even shape-complete.
        if _source_of(d) != 'primary' and not _slot_shape_ok_v17(d, area_limits, seed=True):
            continue
        codes, reasons = [], []
        if pa < SLOT_MIN_AREA_PX:
            codes.append('slot_area_too_small'); reasons.append(f'filled pixel area {pa} px < {SLOT_MIN_AREA_PX} px')
        if np.isfinite(mean_area) and mean_area > 0:
            if pa < area_lo:
                codes.append('slot_area_below_0.1x_mean'); reasons.append(f'area {pa}px < {area_lo:.1f}px')
            if pa > area_hi:
                codes.append('slot_area_above_10x_mean'); reasons.append(f'area {pa}px > {area_hi:.1f}px')
        if reg['regularity_solidity'] < V17_SLOT_SEED_MIN_SOLIDITY:
            codes.append('slot_low_regularity_solidity'); reasons.append(f"solidity {reg['regularity_solidity']:.3f} too low")
        if reg['regularity_convexity'] < V17_SLOT_SEED_MIN_CONVEXITY:
            codes.append('slot_low_convexity'); reasons.append(f"convexity {reg['regularity_convexity']:.3f} too low")
        if reg['regularity_template_iou'] < V17_SLOT_SEED_MIN_TEMPLATE_IOU:
            codes.append('slot_irregular_template_match'); reasons.append(f"template IoU {reg['regularity_template_iou']:.3f} too low")
        if not codes:
            codes = ['slot_center_not_used_by_4row_model']
            reasons = ['regular short-slot candidate center was not the best match to any expected site in the tolerant 4-row array']
        rejects.append(_make_reject_from_contour(
            d.contour, 'slot_v17_4row_filter', '+'.join(codes), '; '.join(reasons), desc=d,
            source=_source_of(d), pixel_area_px=pa, **reg,
        ))

    # Do NOT throw away all data merely because one row remains weak. Keep measured rows,
    # but expose row counts in image_status so the user can inspect the remaining miss.
    return rows, selected, reps, rejects, meta


def measure_via_v14(
    gray: np.ndarray,
    pool: List[ShapeDesc],
    px_nm: float,
    image_name: str,
    condition: str,
) -> Tuple[List[dict], Optional[ShapeDesc], List[ShapeDesc], List[dict], Dict[str, object]]:
    strict = [d for d in pool if _via_strict_metrics_v12(d)[0]]
    relaxed = [d for d in pool if _via_grid_relaxed_ok_v12(d)]
    ms = infer_lattice_v14(strict, gray.shape, "via") if LATTICE_ASSIST_VIA else {"used": False}
    mr = infer_lattice_v14(relaxed, gray.shape, "via") if LATTICE_ASSIST_VIA else {"used": False}
    model, seeds = (mr, relaxed) if _lattice_model_score_v14(mr) > _lattice_model_score_v14(ms) else (ms, strict)

    selected: List[ShapeDesc] = []
    used_ids: set = set()
    recovered = 0
    final_recovered = 0
    attempts = 0
    final_attempts = 0

    if model.get("used", False):
        med_w = float(model["median_w"])
        med_h = float(model["median_h"])
        area_ref = float(np.median([contour_pixel_count(d.contour) for d in seeds])) if seeds else max(30.0, med_w * med_h * 0.7)
        tol = _lattice_match_tol_v14(model, max(med_w, med_h))
        nearest_pitch = min(float(model["basis_a_len"]), float(model["basis_b_len"]))

        # Pass 1: all internal lattice sites, using normal relaxed array standards.
        for site in lattice_sites_v14(model, gray.shape, 0, "via"):
            xy = (site["x"], site["y"])
            d = _select_best_near_lattice_site_v14(pool, used_ids, xy, tol, "via", area_ref)
            if d is None and attempts < MAX_LOCAL_RECOVERY_SITES:
                attempts += 1
                if attempts == 1 or attempts % LOCAL_RECOVERY_PROGRESS_EVERY == 0:
                    print(f"      [{image_name}] lattice via recovery {attempts}/{MAX_LOCAL_RECOVERY_SITES} ...", flush=True)
                d = recover_local_object_v14(gray, xy, med_w, med_h, "via", nearest_pitch, area_ref, None, False)
                if d is not None:
                    recovered += 1
                    pool.append(d)
            if d is not None and not _desc_close_to_any_selected_v12(d, selected):
                selected.append(d); used_ids.add(id(d))

        # Pass 2: refit from accepted points, predict missing equal-pitch sites (including
        # one ring outside current bounds) and lower morphology standards locally.
        model = _refine_lattice_from_selected_v14(model, selected)
        final_sites = lattice_sites_v14(model, gray.shape, LATTICE_FINAL_EXTEND_RINGS, "via")
        for site in final_sites:
            xy = (site["x"], site["y"])
            if any(float(np.linalg.norm(d.center - np.asarray(xy))) <= tol for d in selected):
                continue
            if final_attempts >= LATTICE_FINAL_MAX_SITES:
                break
            final_attempts += 1
            if final_attempts == 1 or final_attempts % LATTICE_FINAL_PROGRESS_EVERY == 0:
                print(f"      [{image_name}] FINAL lattice via search {final_attempts}/{LATTICE_FINAL_MAX_SITES} ...", flush=True)
            d = _select_best_near_lattice_site_v14(pool, used_ids, xy, 1.15 * tol, "via", area_ref, None, True)
            if d is None:
                d = recover_local_object_v14(gray, xy, med_w, med_h, "via", nearest_pitch, area_ref, None, True)
                if d is not None:
                    pool.append(d)
            if d is not None and not _desc_close_to_any_selected_v12(d, selected):
                selected.append(d); used_ids.add(id(d)); final_recovered += 1
    else:
        selected = list(strict)

    # V15 final center-based recursive completion. Predictions start from actual
    # accepted contour CENTERS and every newly recovered center can expand further.
    center_meta = {}
    if CENTER_RECURSIVE_SEARCH_VIA and len(selected) >= 2:
        selected, center_meta = recursive_center_recovery_v15(
            gray, pool, selected, "via", model, image_name, None
        )
        if model.get("used", False):
            model = _refine_lattice_from_selected_v14(model, selected)

    rows: List[dict] = []
    for oid, d in enumerate(sorted(selected, key=lambda z: (z.center[1], z.center[0])), 1):
        pts = d.contour.reshape(-1, 2).astype(np.float64)
        w = float(pts[:, 0].max() - pts[:, 0].min())
        h = float(pts[:, 1].max() - pts[:, 1].min())
        rows.append({
            "condition": condition, "image": image_name, "pattern": "via", "object_id": oid,
            "width_nm": w * px_nm, "height_nm": h * px_nm,
            "height_width_ratio": h / max(w, 1e-12), "circularity": d.circularity,
            "solidity": d.solidity, "axis_ratio": d.aspect, "selection_source": _source_of(d),
            "pixel_size_nm": px_nm,
        })

    rejected: List[dict] = []
    for d in pool:
        if _desc_close_to_any_selected_v12(d, selected):
            continue
        strict_ok, codes, reasons = _via_strict_metrics_v12(d)
        if _source_of(d) != "primary" and not _via_grid_relaxed_ok_v12(d):
            continue
        if strict_ok and model.get("used", False):
            codes = ["via_off_array_grid"]; reasons = ["via-like candidate does not coincide with fitted 2D lattice"]
        elif not codes:
            codes = ["via_not_selected"]; reasons = ["candidate was not selected"]
        rejected.append(_make_reject_from_contour(
            d.contour, "via_v14_filter", "+".join(codes), "; ".join(reasons), desc=d,
            source=_source_of(d), circularity=d.circularity, solidity=d.solidity, axis_ratio=d.aspect,
        ))

    representative = max(selected, key=lambda d: (d.circularity, d.solidity)) if selected else None
    internal_n = len(lattice_sites_v14(model, gray.shape, 0, "via")) if model.get("used", False) else 0
    final_n = len(lattice_sites_v14(model, gray.shape, LATTICE_FINAL_EXTEND_RINGS, "via")) if model.get("used", False) else 0
    meta = {
        "array_assist_used": bool(model.get("used", False)),
        "array_expected_positions": internal_n,
        "array_recovered_count": recovered + final_recovered,
        "array_initial_recovered_count": recovered,
        "array_final_recovered_count": final_recovered,
        "array_final_predicted_positions": final_n,
        "array_x_pitch_px": np.nan,
        "array_y_pitch_px": np.nan,
        "array_occupancy_seed": float(model.get("occupancy", np.nan)),
        "array_recovery_attempts": attempts,
        "array_final_recovery_attempts": final_attempts,
        "array_recovery_skipped_due_cap": 0,
        "array_lattice_type": model.get("lattice_type", "none"),
        "array_basis_a_dx_px": float(np.asarray(model.get("a", [np.nan, np.nan]))[0]),
        "array_basis_a_dy_px": float(np.asarray(model.get("a", [np.nan, np.nan]))[1]),
        "array_basis_b_dx_px": float(np.asarray(model.get("b", [np.nan, np.nan]))[0]),
        "array_basis_b_dy_px": float(np.asarray(model.get("b", [np.nan, np.nan]))[1]),
        "array_basis_a_length_px": float(model.get("basis_a_len", np.nan)),
        "array_basis_b_length_px": float(model.get("basis_b_len", np.nan)),
        "array_basis_angle_deg": float(model.get("basis_angle_deg", np.nan)),
    }
    meta.update(center_meta)
    return rows, representative, selected, rejected, meta


def measure_slot_v14(
    gray: np.ndarray,
    pool: List[ShapeDesc],
    px_nm: float,
    image_name: str,
    condition: str,
) -> Tuple[List[dict], List[ShapeDesc], Dict[int, ShapeDesc], List[dict], Dict[str, object]]:
    strict, rejects, reg_by_id, mean_area, med_area, area_hi = _slot_initial_filter_v12(pool)
    if np.isfinite(mean_area) and mean_area > 0:
        area_lo = max(SLOT_MIN_AREA_PX, mean_area * SLOT_AREA_MIN_MEAN_FACTOR)
        area_hi = mean_area * SLOT_AREA_MAX_MEAN_FACTOR
    else:
        area_lo, area_hi = SLOT_MIN_AREA_PX, np.inf

    relaxed = []
    for d in pool:
        pa = contour_pixel_count(d.contour)
        ok, _, _, _ = _slot_regularity_ok_v12(d, relaxed=True)
        if pa >= SLOT_MIN_AREA_PX and ok:
            relaxed.append(d)
    ms = infer_lattice_v14(strict, gray.shape, "slot") if LATTICE_ASSIST_SLOT else {"used": False}
    mr = infer_lattice_v14(relaxed, gray.shape, "slot") if LATTICE_ASSIST_SLOT else {"used": False}
    model, seeds = (mr, relaxed) if _lattice_model_score_v14(mr) > _lattice_model_score_v14(ms) else (ms, strict)

    if (not np.isfinite(mean_area) or mean_area <= 0) and seeds:
        aa = np.asarray([contour_pixel_count(d.contour) for d in seeds], dtype=float)
        med_area = float(np.median(aa)); mean_area = float(np.mean(aa))
        area_lo = max(SLOT_MIN_AREA_PX, mean_area * SLOT_AREA_MIN_MEAN_FACTOR)
        area_hi = mean_area * SLOT_AREA_MAX_MEAN_FACTOR

    selected: List[ShapeDesc] = []
    row_ids: List[int] = []
    used_ids: set = set()
    recovered = 0
    final_recovered = 0
    attempts = 0
    final_attempts = 0

    if model.get("used", False):
        row_dims = _slot_row_dims_v14(model, seeds, mean_area)
        tol = _lattice_match_tol_v14(model, float(model.get("median_w", 10.0)))
        nearest_pitch = min(float(model["basis_a_len"]), float(model["basis_b_len"]))
        limits = (area_lo, area_hi)

        for site in lattice_sites_v14(model, gray.shape, 0, "slot"):
            rid = int(site["row_id"])
            rw, rh, rarea = row_dims.get(rid, (float(model["median_w"]), float(model["median_h"]), float(mean_area)))
            xy = (site["x"], site["y"])
            d = _select_best_near_lattice_site_v14(pool, used_ids, xy, tol, "slot", rarea, limits, False)
            if d is None and attempts < MAX_LOCAL_RECOVERY_SITES:
                attempts += 1
                if attempts == 1 or attempts % LOCAL_RECOVERY_PROGRESS_EVERY == 0:
                    print(f"      [{image_name}] lattice slot recovery {attempts}/{MAX_LOCAL_RECOVERY_SITES} ...", flush=True)
                d = recover_local_object_v14(gray, xy, rw, rh, "slot", nearest_pitch, rarea, limits, False)
                if d is not None:
                    recovered += 1; pool.append(d); reg_by_id[id(d)] = slot_regularity_metrics(d)
            if d is not None and not _desc_close_to_any_selected_v12(d, selected):
                selected.append(d); row_ids.append(rid); used_ids.add(id(d))

        model = _refine_lattice_from_selected_v14(model, selected)
        row_dims = _slot_row_dims_v14(model, selected if selected else seeds, mean_area)
        final_sites = lattice_sites_v14(model, gray.shape, LATTICE_FINAL_EXTEND_RINGS, "slot")
        for site in final_sites:
            rid = int(site["row_id"])
            # Ignore any predicted row outside S1..S4 (should not occur for slot).
            if rid < 1 or rid > 4:
                continue
            xy = (site["x"], site["y"])
            if any(float(np.linalg.norm(d.center - np.asarray(xy))) <= tol for d in selected):
                continue
            if final_attempts >= LATTICE_FINAL_MAX_SITES:
                break
            final_attempts += 1
            if final_attempts == 1 or final_attempts % LATTICE_FINAL_PROGRESS_EVERY == 0:
                print(f"      [{image_name}] FINAL lattice slot search {final_attempts}/{LATTICE_FINAL_MAX_SITES} ...", flush=True)
            rw, rh, rarea = row_dims.get(rid, (float(model["median_w"]), float(model["median_h"]), float(mean_area)))
            d = _select_best_near_lattice_site_v14(pool, used_ids, xy, 1.15 * tol, "slot", rarea, limits, True)
            if d is None:
                d = recover_local_object_v14(gray, xy, rw, rh, "slot", nearest_pitch, rarea, limits, True)
                if d is not None:
                    pool.append(d); reg_by_id[id(d)] = slot_regularity_metrics(d)
            if d is not None and not _desc_close_to_any_selected_v12(d, selected):
                selected.append(d); row_ids.append(rid); used_ids.add(id(d)); final_recovered += 1
    else:
        selected = list(strict)
        if len(selected) >= 4:
            labs = cluster_1d_into_four_rows(np.asarray([d.center[1] for d in selected], dtype=np.float32))
            row_ids = [int(x) for x in labs]

    # V15 recursive CENTER completion. Slot area hard limits are retained, but
    # morphology is deliberately looser at these strongly predicted positions.
    center_meta = {}
    if CENTER_RECURSIVE_SEARCH_SLOT and len(selected) >= 2:
        selected, center_meta = recursive_center_recovery_v15(
            gray, pool, selected, "slot", model, image_name, (area_lo, area_hi)
        )
        if model.get("used", False):
            model = _refine_lattice_from_selected_v14(model, selected)
        row_ids = assign_slot_rows_v15(model, selected)

    if len(selected) < SLOT_REQUIRE_AT_LEAST_N or len(row_ids) != len(selected):
        for d in selected:
            rejects.append(_make_reject_from_contour(
                d.contour, "slot_group_filter", "slot_too_few_valid_objects",
                f"only {len(selected)} selected slot(s); need >= {SLOT_REQUIRE_AT_LEAST_N}", desc=d,
                source=_source_of(d),
            ))
        selected, row_ids = [], []

    rows: List[dict] = []
    reps: Dict[int, ShapeDesc] = {}
    if selected:
        order = sorted(range(len(selected)), key=lambda k: (row_ids[k], selected[k].center[0], selected[k].center[1]))
        selected = [selected[k] for k in order]
        row_ids = [row_ids[k] for k in order]
        for oid, (d, rid) in enumerate(zip(selected, row_ids), 1):
            pts = d.contour.reshape(-1, 2).astype(np.float64)
            width_px = float(pts[:, 0].max() - pts[:, 0].min())
            length_px = float(pts[:, 1].max() - pts[:, 1].min())
            reg = reg_by_id.get(id(d), slot_regularity_metrics(d))
            pa = contour_pixel_count(d.contour)
            rows.append({
                "condition": condition, "image": image_name, "pattern": "slot", "object_id": oid,
                "slot_row": int(rid), "slot_width_nm": width_px * px_nm, "slot_length_nm": length_px * px_nm,
                "slot_aspect_ratio": length_px / max(width_px, 1e-12),
                "slot_angle_from_y_deg": abs(90.0 - d.angle_from_x_deg),
                "slot_pixel_area_px": pa, "slot_mean_candidate_area_px": mean_area,
                "slot_reference_median_area_px": med_area,
                "slot_area_ratio_to_mean": pa / mean_area if np.isfinite(mean_area) and mean_area > 0 else np.nan,
                "slot_regularity_solidity": reg["regularity_solidity"],
                "slot_regularity_convexity": reg["regularity_convexity"],
                "slot_rect_iou": reg["regularity_rect_iou"], "slot_ellipse_iou": reg["regularity_ellipse_iou"],
                "slot_best_template_iou": reg["regularity_template_iou"], "selection_source": _source_of(d),
                "pixel_size_nm": px_nm,
            })
            if rid not in reps:
                reps[rid] = d

    rejects = [r for r in rejects if r.get("_desc") is None or not _desc_close_to_any_selected_v12(r["_desc"], selected)]
    rejected_desc_ids = {id(r.get("_desc")) for r in rejects if r.get("_desc") is not None}
    for d in pool:
        if _desc_close_to_any_selected_v12(d, selected) or id(d) in rejected_desc_ids:
            continue
        pa = contour_pixel_count(d.contour)
        # Do not flood logs with weak enhanced fragments that never looked structure-like.
        ok_rel, reg_rel, _, _ = _slot_regularity_ok_v12(d, relaxed=True)
        if _source_of(d) != "primary" and not (pa >= SLOT_MIN_AREA_PX and ok_rel):
            continue
        codes, reasons = [], []
        if pa < SLOT_MIN_AREA_PX:
            codes.append("slot_area_too_small"); reasons.append(f"filled pixel area {pa} px < {SLOT_MIN_AREA_PX} px")
        if np.isfinite(mean_area) and mean_area > 0:
            if pa < area_lo:
                codes.append("slot_area_below_0.1x_mean"); reasons.append(f"area {pa}px < {area_lo:.1f}px")
            if pa > area_hi:
                codes.append("slot_area_above_10x_mean"); reasons.append(f"area {pa}px > {area_hi:.1f}px")
        if not codes and model.get("used", False):
            codes = ["slot_off_array_grid"]; reasons = ["slot-like candidate does not coincide with fitted 2D lattice"]
        elif not codes:
            codes = ["slot_not_selected"]; reasons = ["candidate was not selected"]
        rejects.append(_make_reject_from_contour(
            d.contour, "slot_v14_filter", "+".join(codes), "; ".join(reasons), desc=d,
            source=_source_of(d), pixel_area_px=pa, **reg_rel,
        ))

    internal_n = len(lattice_sites_v14(model, gray.shape, 0, "slot")) if model.get("used", False) else 0
    final_n = len(lattice_sites_v14(model, gray.shape, LATTICE_FINAL_EXTEND_RINGS, "slot")) if model.get("used", False) else 0
    meta = {
        "array_assist_used": bool(model.get("used", False)),
        "array_expected_positions": internal_n,
        "array_recovered_count": recovered + final_recovered,
        "array_initial_recovered_count": recovered,
        "array_final_recovered_count": final_recovered,
        "array_final_predicted_positions": final_n,
        "array_x_pitch_px": np.nan, "array_y_pitch_px": np.nan,
        "array_occupancy_seed": float(model.get("occupancy", np.nan)),
        "array_recovery_attempts": attempts, "array_final_recovery_attempts": final_attempts,
        "array_recovery_skipped_due_cap": 0,
        "array_lattice_type": model.get("lattice_type", "none"),
        "array_basis_a_dx_px": float(np.asarray(model.get("a", [np.nan, np.nan]))[0]),
        "array_basis_a_dy_px": float(np.asarray(model.get("a", [np.nan, np.nan]))[1]),
        "array_basis_b_dx_px": float(np.asarray(model.get("b", [np.nan, np.nan]))[0]),
        "array_basis_b_dy_px": float(np.asarray(model.get("b", [np.nan, np.nan]))[1]),
        "array_basis_a_length_px": float(model.get("basis_a_len", np.nan)),
        "array_basis_b_length_px": float(model.get("basis_b_len", np.nan)),
        "array_basis_angle_deg": float(model.get("basis_angle_deg", np.nan)),
    }
    meta.update(center_meta)
    return rows, selected, reps, rejects, meta


def measure_via_v12(
    gray: np.ndarray,
    pool: List[ShapeDesc],
    px_nm: float,
    image_name: str,
    condition: str,
) -> Tuple[List[dict], Optional[ShapeDesc], List[ShapeDesc], List[dict], Dict[str, object]]:
    strict = [d for d in pool if _via_strict_metrics_v12(d)[0]]
    relaxed_grid_seeds = [d for d in pool if _via_grid_relaxed_ok_v12(d)]
    grid_strict = infer_via_grid_v12(strict, gray.shape) if ARRAY_ASSIST_VIA else {"used": False}
    grid_relaxed = (
        infer_via_grid_v12(relaxed_grid_seeds, gray.shape)
        if ARRAY_ASSIST_VIA and len(relaxed_grid_seeds) >= ARRAY_MIN_GRID_SEEDS else {"used": False}
    )
    if _grid_model_score_v12(grid_relaxed, len(relaxed_grid_seeds)) > _grid_model_score_v12(grid_strict, len(strict)):
        grid, grid_seeds = grid_relaxed, relaxed_grid_seeds
    else:
        grid, grid_seeds = grid_strict, strict
    selected: List[ShapeDesc] = []
    recovered_count = 0
    used_ids: set = set()

    if grid.get("used", False):
        xp = np.asarray(grid["x_positions"], dtype=float)
        yp = np.asarray(grid["y_positions"], dtype=float)
        x_pitch = float(grid["x_pitch"])
        y_pitch = float(grid["y_pitch"])
        med_w = float(grid["median_w"])
        med_h = float(grid["median_h"])
        area_ref = float(np.median([contour_pixel_count(d.contour) for d in grid_seeds]))
        tol_x = _grid_match_tolerance_v12(x_pitch, med_w)
        tol_y = _grid_match_tolerance_v12(y_pitch, med_h)

        recovery_attempts = 0
        recovery_skipped = 0
        for y in yp:
            for x in xp:
                d = _select_best_near_site_v12(
                    pool, used_ids, (x, y), tol_x, tol_y,
                    "via", area_ref, None,
                )
                if d is None:
                    if recovery_attempts < MAX_LOCAL_RECOVERY_SITES:
                        recovery_attempts += 1
                        if recovery_attempts == 1 or recovery_attempts % LOCAL_RECOVERY_PROGRESS_EVERY == 0:
                            print(
                                f"      [{image_name}] local via recovery {recovery_attempts}/{MAX_LOCAL_RECOVERY_SITES} ...",
                                flush=True,
                            )
                        d = recover_local_object_v12(
                            gray, (x, y), med_w, med_h, "via",
                            x_pitch, y_pitch, area_ref,
                        )
                        if d is not None:
                            recovered_count += 1
                            pool.append(d)
                    else:
                        recovery_skipped += 1
                if d is not None and not _desc_close_to_any_selected_v12(d, selected):
                    selected.append(d)
                    used_ids.add(id(d))
    else:
        selected = list(strict)

    rows: List[dict] = []
    for i, d in enumerate(sorted(selected, key=lambda z: (z.center[1], z.center[0])), 1):
        pts = d.contour.reshape(-1, 2).astype(np.float64)
        width_px_xy = float(pts[:, 0].max() - pts[:, 0].min())
        height_px_xy = float(pts[:, 1].max() - pts[:, 1].min())
        rows.append({
            "condition": condition, "image": image_name, "pattern": "via", "object_id": i,
            "width_nm": width_px_xy * px_nm,
            "height_nm": height_px_xy * px_nm,
            "height_width_ratio": height_px_xy / max(width_px_xy, 1e-12),
            "circularity": d.circularity, "solidity": d.solidity, "axis_ratio": d.aspect,
            "selection_source": _source_of(d), "pixel_size_nm": px_nm,
        })

    rejected: List[dict] = []
    for d in pool:
        if _desc_close_to_any_selected_v12(d, selected):
            continue
        strict_ok, codes, reasons = _via_strict_metrics_v12(d)
        # Enhanced masks are exploratory; do not flood QC output with every weak
        # threshold fragment. Only log an enhanced object if it is itself via-like.
        if _source_of(d) != "primary" and not _via_grid_relaxed_ok_v12(d):
            continue
        if strict_ok and grid.get("used", False):
            codes = ["via_off_array_grid"]
            reasons = ["via-like candidate does not coincide with a fitted equal-pitch array site"]
        elif not codes:
            codes = ["via_not_selected"]
            reasons = ["candidate was not selected"]
        rejected.append(_make_reject_from_contour(
            d.contour, "via_v13_filter", "+".join(codes), "; ".join(reasons), desc=d,
            source=_source_of(d), circularity=d.circularity, solidity=d.solidity, axis_ratio=d.aspect,
        ))

    representative = max(selected, key=lambda d: (d.circularity, d.solidity)) if selected else None
    meta = {
        "array_assist_used": bool(grid.get("used", False)),
        "array_expected_positions": int(len(grid.get("x_positions", [])) * len(grid.get("y_positions", []))) if grid.get("used", False) else 0,
        "array_recovered_count": recovered_count,
        "array_x_pitch_px": float(grid.get("x_pitch", np.nan)),
        "array_y_pitch_px": float(grid.get("y_pitch", np.nan)),
        "array_occupancy_seed": float(grid.get("occupancy", np.nan)),
        "array_recovery_attempts": int(locals().get("recovery_attempts", 0)),
        "array_recovery_skipped_due_cap": int(locals().get("recovery_skipped", 0)),
    }
    return rows, representative, selected, rejected, meta


def _slot_initial_filter_v12(pool: List[ShapeDesc]) -> Tuple[List[ShapeDesc], List[dict], Dict[int, Dict[str, float]], float, float, float]:
    """Strict regularity + 30-pixel + relative-area QC; returns candidates before array fitting."""
    eligible: List[ShapeDesc] = []
    rejects: List[dict] = []
    reg_by_id: Dict[int, Dict[str, float]] = {}
    for d in pool:
        pixel_area = contour_pixel_count(d.contour)
        ok_reg, reg, codes, reasons = _slot_regularity_ok_v12(d, relaxed=False)
        reg_by_id[id(d)] = reg
        if pixel_area < SLOT_MIN_AREA_PX:
            codes = ["slot_area_too_small"] + codes
            reasons = [f"filled pixel area {pixel_area} px < {SLOT_MIN_AREA_PX} px"] + reasons
        if codes:
            # Primary-mask candidates are real QC observations. Secondary enhanced
            # masks are exploratory and can create many fragments, so failed enhanced
            # candidates are silently ignored unless they later look structure-like.
            if _source_of(d) == "primary":
                rejects.append(_make_reject_from_contour(
                    d.contour, "slot_shape_filter", "+".join(codes), "; ".join(reasons), desc=d,
                    source=_source_of(d), pixel_area_px=pixel_area, **reg,
                ))
            continue
        eligible.append(d)

    mean_area = np.nan
    med_area = np.nan
    lo = 0.0
    hi = np.inf
    if eligible:
        areas = np.asarray([contour_pixel_count(d.contour) for d in eligible], dtype=float)
        med_area = float(np.median(areas))
        seed = areas[(areas >= med_area * SLOT_AREA_MIN_MEAN_FACTOR) & (areas <= med_area * SLOT_AREA_MAX_MEAN_FACTOR)]
        if len(seed) == 0:
            seed = areas
        mean_area = float(np.mean(seed))
        lo = mean_area * SLOT_AREA_MIN_MEAN_FACTOR
        hi = mean_area * SLOT_AREA_MAX_MEAN_FACTOR

    kept: List[ShapeDesc] = []
    for d in eligible:
        a = float(contour_pixel_count(d.contour))
        if SLOT_USE_RELATIVE_AREA_FILTER and np.isfinite(mean_area):
            if a < lo:
                if _source_of(d) == "primary":
                    rejects.append(_make_reject_from_contour(
                        d.contour, "slot_relative_area_filter", "slot_area_below_0.1x_mean",
                        f"filled pixel area {a:.0f} px < 0.1x mean {lo:.2f} px", desc=d,
                        source=_source_of(d), pixel_area_px=a, slot_mean_area_px=mean_area,
                        area_ratio_to_mean=a / max(mean_area, 1e-9),
                    ))
                continue
            if a > hi:
                if _source_of(d) == "primary":
                    rejects.append(_make_reject_from_contour(
                        d.contour, "slot_relative_area_filter", "slot_area_above_10x_mean",
                        f"filled pixel area {a:.0f} px > 10x mean {hi:.2f} px", desc=d,
                        source=_source_of(d), pixel_area_px=a, slot_mean_area_px=mean_area,
                        area_ratio_to_mean=a / max(mean_area, 1e-9),
                    ))
                continue
        kept.append(d)
    return kept, rejects, reg_by_id, mean_area, med_area, hi


def measure_slot_v12(
    gray: np.ndarray,
    pool: List[ShapeDesc],
    px_nm: float,
    image_name: str,
    condition: str,
) -> Tuple[List[dict], List[ShapeDesc], Dict[int, ShapeDesc], List[dict], Dict[str, object]]:
    strict, rejects, reg_by_id, mean_area, med_area, area_hi = _slot_initial_filter_v12(pool)
    area_lo = mean_area * SLOT_AREA_MIN_MEAN_FACTOR if np.isfinite(mean_area) else SLOT_MIN_AREA_PX
    if not np.isfinite(area_hi):
        area_hi = np.inf

    # For a blurry image, strict regularity may reject nearly everything. Use a relaxed
    # regularity set only to infer the array; final acceptance still requires proximity
    # to a supported array site and reasonable geometry.
    relaxed_grid_seeds = []
    for d in pool:
        pa = contour_pixel_count(d.contour)
        ok, _, _, _ = _slot_regularity_ok_v12(d, relaxed=True)
        if pa >= SLOT_MIN_AREA_PX and ok:
            relaxed_grid_seeds.append(d)
    grid_strict = infer_slot_grid_v12(strict, gray.shape) if ARRAY_ASSIST_SLOT else {"used": False}
    grid_relaxed = (
        infer_slot_grid_v12(relaxed_grid_seeds, gray.shape)
        if ARRAY_ASSIST_SLOT and len(relaxed_grid_seeds) >= ARRAY_MIN_GRID_SEEDS else {"used": False}
    )
    if _grid_model_score_v12(grid_relaxed, len(relaxed_grid_seeds)) > _grid_model_score_v12(grid_strict, len(strict)):
        grid, grid_seeds = grid_relaxed, relaxed_grid_seeds
    else:
        grid, grid_seeds = grid_strict, strict

    selected: List[ShapeDesc] = []
    selected_row_ids: List[int] = []
    recovered_count = 0
    used_ids: set = set()

    if grid.get("used", False):
        xp = np.asarray(grid["x_positions"], dtype=float)
        yp = np.asarray(grid["y_positions"], dtype=float)
        x_pitch = float(grid["x_pitch"])
        y_pitch = float(grid["y_pitch"])
        # Use high-confidence candidates for typical physical size even if the
        # ARRAY POSITIONS came from relaxed/enhanced candidates. This prevents tiny
        # secondary-threshold fragments from shrinking the recovery ROI.
        dim_ref = strict if len(strict) >= 2 else grid_seeds
        if dim_ref:
            med_w = float(np.median([
                max(1.0, _contour_xy_bounds(d)[1] - _contour_xy_bounds(d)[0]) for d in dim_ref
            ]))
            med_h_global = float(np.median([
                max(1.0, _contour_xy_bounds(d)[3] - _contour_xy_bounds(d)[2]) for d in dim_ref
            ]))
        else:
            med_w = float(grid["median_w"])
            med_h_global = float(grid["median_h"])

        if np.isfinite(mean_area) and mean_area > 0:
            area_ref_global = float(mean_area)
        else:
            area_ref_global = float(np.median([contour_pixel_count(d.contour) for d in grid_seeds]))
            mean_area = area_ref_global
            med_area = float(np.median([contour_pixel_count(d.contour) for d in grid_seeds]))
            area_lo = max(SLOT_MIN_AREA_PX, area_ref_global * SLOT_AREA_MIN_MEAN_FACTOR)
            area_hi = area_ref_global * SLOT_AREA_MAX_MEAN_FACTOR

        # Row-specific dimensions are helpful because slot 1..4 have different Y lengths.
        row_dims: Dict[int, Tuple[float, float, float]] = {}
        for rid, y in enumerate(yp, 1):
            near0 = [d for d in grid_seeds if abs(d.center[1] - y) <= _grid_match_tolerance_v12(y_pitch, med_h_global)]
            # Row-specific size is useful only when supported by plausible objects.
            # A single tiny secondary-threshold fragment must not define the local-search ROI.
            near = []
            for d in near0:
                ww = max(1.0, _contour_xy_bounds(d)[1] - _contour_xy_bounds(d)[0])
                pa = contour_pixel_count(d.contour)
                min_row_area = max(SLOT_MIN_AREA_PX, 0.20 * area_ref_global)
                if 0.50 * med_w <= ww <= 2.5 * med_w and pa >= min_row_area:
                    near.append(d)
            if len(near) >= 2:
                ws = [max(1.0, _contour_xy_bounds(d)[1] - _contour_xy_bounds(d)[0]) for d in near]
                hs = [max(1.0, _contour_xy_bounds(d)[3] - _contour_xy_bounds(d)[2]) for d in near]
                ars = [contour_pixel_count(d.contour) for d in near]
                row_dims[rid] = (float(np.median(ws)), float(np.median(hs)), float(np.median(ars)))
            else:
                row_dims[rid] = (med_w, med_h_global, area_ref_global)

        tol_x = _grid_match_tolerance_v12(x_pitch, med_w)
        # Center alignment in Y should be much tighter than slot length; cap the tolerance.
        tol_y = min(_grid_match_tolerance_v12(y_pitch, max(4.0, 0.25 * med_h_global)), 0.38 * y_pitch if np.isfinite(y_pitch) else 9999)

        recovery_attempts = 0
        recovery_skipped = 0
        for rid, y in enumerate(yp, 1):
            rw, rh, rarea = row_dims[rid]
            for x in xp:
                d = _select_best_near_site_v12(
                    pool, used_ids, (x, y), tol_x, tol_y,
                    "slot", rarea, (area_lo, area_hi),
                )
                if d is None:
                    if recovery_attempts < MAX_LOCAL_RECOVERY_SITES:
                        recovery_attempts += 1
                        if recovery_attempts == 1 or recovery_attempts % LOCAL_RECOVERY_PROGRESS_EVERY == 0:
                            print(
                                f"      [{image_name}] local slot recovery {recovery_attempts}/{MAX_LOCAL_RECOVERY_SITES} ...",
                                flush=True,
                            )
                        d = recover_local_object_v12(
                            gray, (x, y), rw, rh, "slot",
                            x_pitch, y_pitch, rarea,
                        )
                        if d is not None:
                            pa = contour_pixel_count(d.contour)
                            if pa < max(SLOT_MIN_AREA_PX, area_lo) or pa > area_hi:
                                d = None
                        if d is not None:
                            recovered_count += 1
                            pool.append(d)
                            ok, reg, _, _ = _slot_regularity_ok_v12(d, relaxed=True)
                            reg_by_id[id(d)] = reg
                    else:
                        recovery_skipped += 1
                if d is not None and not _desc_close_to_any_selected_v12(d, selected):
                    selected.append(d)
                    selected_row_ids.append(rid)
                    used_ids.add(id(d))
    else:
        # Legacy fallback if the equal-pitch grid cannot be fitted.
        selected = list(strict)
        if len(selected) >= 4:
            labs = cluster_1d_into_four_rows(np.asarray([d.center[1] for d in selected], dtype=np.float32))
            selected_row_ids = [int(x) for x in labs]
        else:
            selected_row_ids = []

    if len(selected) < SLOT_REQUIRE_AT_LEAST_N or len(selected_row_ids) != len(selected):
        # Keep earlier reason records and add a clear image-level structure reason.
        for d in selected:
            rejects.append(_make_reject_from_contour(
                d.contour, "slot_group_filter", "slot_too_few_valid_objects",
                f"only {len(selected)} selected slot(s); need >= {SLOT_REQUIRE_AT_LEAST_N}", desc=d,
                source=_source_of(d),
            ))
        selected = []
        selected_row_ids = []

    rows: List[dict] = []
    if selected:
        order = sorted(range(len(selected)), key=lambda i: (selected_row_ids[i], selected[i].center[0]))
        selected = [selected[i] for i in order]
        selected_row_ids = [selected_row_ids[i] for i in order]
        for oid, (d, rid) in enumerate(zip(selected, selected_row_ids), 1):
            pts = d.contour.reshape(-1, 2).astype(np.float64)
            width_px = float(pts[:, 0].max() - pts[:, 0].min())
            length_px = float(pts[:, 1].max() - pts[:, 1].min())
            reg = reg_by_id.get(id(d))
            if reg is None:
                reg = slot_regularity_metrics(d)
                reg_by_id[id(d)] = reg
            pa = contour_pixel_count(d.contour)
            rows.append({
                "condition": condition, "image": image_name, "pattern": "slot", "object_id": oid,
                "slot_row": int(rid), "slot_width_nm": width_px * px_nm,
                "slot_length_nm": length_px * px_nm,
                "slot_aspect_ratio": length_px / max(width_px, 1e-12),
                "slot_angle_from_y_deg": abs(90.0 - d.angle_from_x_deg),
                "slot_pixel_area_px": pa,
                "slot_mean_candidate_area_px": mean_area,
                "slot_reference_median_area_px": med_area,
                "slot_area_ratio_to_mean": pa / mean_area if np.isfinite(mean_area) and mean_area > 0 else np.nan,
                "slot_regularity_solidity": reg["regularity_solidity"],
                "slot_regularity_convexity": reg["regularity_convexity"],
                "slot_rect_iou": reg["regularity_rect_iou"],
                "slot_ellipse_iou": reg["regularity_ellipse_iou"],
                "slot_best_template_iou": reg["regularity_template_iou"],
                "selection_source": _source_of(d), "pixel_size_nm": px_nm,
            })

    # A relaxed candidate can be rescued by array position even if strict regularity
    # rejected it earlier. Do not show the same physical object as both GREEN and RED.
    rejects = [
        r for r in rejects
        if r.get("_desc") is None or not _desc_close_to_any_selected_v12(r["_desc"], selected)
    ]

    # Candidates that passed strict shape QC but do not sit on the fitted lattice are
    # explicitly rejected as off-grid. Existing shape rejects are preserved.
    rejected_desc_ids = {id(r.get("_desc")) for r in rejects if r.get("_desc") is not None}
    for d in pool:
        if _desc_close_to_any_selected_v12(d, selected):
            continue
        if id(d) in rejected_desc_ids:
            continue
        ok_reg, reg, codes, reasons = _slot_regularity_ok_v12(d, relaxed=False)
        pa = contour_pixel_count(d.contour)
        if pa < SLOT_MIN_AREA_PX:
            codes = ["slot_area_too_small"] + codes
            reasons = [f"filled pixel area {pa} px < {SLOT_MIN_AREA_PX} px"] + reasons
        if _source_of(d) != "primary":
            relaxed_ok, _, _, _ = _slot_regularity_ok_v12(d, relaxed=True)
            if not relaxed_ok or pa < SLOT_MIN_AREA_PX:
                continue
        if not codes and grid.get("used", False):
            codes = ["slot_off_array_grid"]
            reasons = ["regular slot-like candidate does not coincide with fitted equal-pitch array sites"]
        elif not codes:
            codes = ["slot_not_selected"]
            reasons = ["candidate was not selected"]
        rejects.append(_make_reject_from_contour(
            d.contour, "slot_v13_filter", "+".join(codes), "; ".join(reasons), desc=d,
            source=_source_of(d), pixel_area_px=pa, **reg,
        ))

    reps: Dict[int, ShapeDesc] = {}
    if selected:
        xc = float(np.median([d.center[0] for d in selected]))
        for rid in [1, 2, 3, 4]:
            ds = [d for d, r in zip(selected, selected_row_ids) if r == rid]
            if ds:
                reps[rid] = min(ds, key=lambda d: abs(d.center[0] - xc))

    meta = {
        "array_assist_used": bool(grid.get("used", False)),
        "array_expected_positions": int(len(grid.get("x_positions", [])) * len(grid.get("y_positions", []))) if grid.get("used", False) else 0,
        "array_recovered_count": recovered_count,
        "array_x_pitch_px": float(grid.get("x_pitch", np.nan)),
        "array_y_pitch_px": float(grid.get("y_pitch", np.nan)),
        "array_occupancy_seed": float(grid.get("occupancy", np.nan)),
        "array_recovery_attempts": int(locals().get("recovery_attempts", 0)),
        "array_recovery_skipped_due_cap": int(locals().get("recovery_skipped", 0)),
    }
    return rows, selected, reps, rejects, meta


def annotate_selected_vias_v12(img: np.ndarray, selected: List[ShapeDesc]) -> np.ndarray:
    if not ANNOTATE_SELECTED_OBJECTS or not selected:
        return img
    out = img.copy()
    green = (0, 255, 0)
    for i, d in enumerate(sorted(selected, key=lambda z: (z.center[1], z.center[0])), 1):
        cv2.drawContours(out, [d.contour.astype(np.int32)], -1, green, SELECTED_CONTOUR_THICKNESS, cv2.LINE_AA)
        src = _source_of(d)
        suffix = (":C" if src in {"center_recursive_recovered", "center_recursive_selected", "fixed_center_recovered", "fixed_center_selected", "slot_4row_recovered", "slot_4row_selected"} else
                  (":F" if src == "lattice_final_recovered" else
                   (":G" if src in {"grid_recovered", "lattice_recovered"} else
                    (":E" if src != "primary" else ""))))
        _draw_green_label(out, d.contour, f"A{i}{suffix}")
    return out


def annotate_selected_slots_v12(img: np.ndarray, selected: List[ShapeDesc], rows: List[dict]) -> np.ndarray:
    if not ANNOTATE_SELECTED_OBJECTS or not selected:
        return img
    out = img.copy()
    green = (0, 255, 0)
    for i, d in enumerate(selected, 1):
        cv2.drawContours(out, [d.contour.astype(np.int32)], -1, green, SELECTED_CONTOUR_THICKNESS, cv2.LINE_AA)
        rid = rows[i - 1].get("slot_row", "?") if i - 1 < len(rows) else "?"
        src = _source_of(d)
        suffix = (":C" if src in {"center_recursive_recovered", "center_recursive_selected", "fixed_center_recovered", "fixed_center_selected", "slot_4row_recovered", "slot_4row_selected"} else
                  (":F" if src == "lattice_final_recovered" else
                   (":G" if src in {"grid_recovered", "lattice_recovered"} else
                    (":E" if src != "primary" else ""))))
        _draw_green_label(out, d.contour, f"A{i}:S{rid}{suffix}")
    return out


# ============================================================
# 8) ANNOTATION HELPERS
# ============================================================
def to_color(gray: np.ndarray) -> np.ndarray:
    return cv2.cvtColor(gray, cv2.COLOR_GRAY2BGR)


def p2i(p: np.ndarray) -> Tuple[int, int]:
    return int(round(float(p[0]))), int(round(float(p[1])))


def draw_double_arrow(img: np.ndarray, p1: np.ndarray, p2: np.ndarray, text: str, text_offset=(5, -5)) -> None:
    cv2.arrowedLine(img, p2i(p1), p2i(p2), (255, 255, 255), 1, cv2.LINE_AA, tipLength=0.08)
    cv2.arrowedLine(img, p2i(p2), p2i(p1), (255, 255, 255), 1, cv2.LINE_AA, tipLength=0.08)
    mid = 0.5 * (p1 + p2)
    cv2.putText(img, text,
                (int(mid[0] + text_offset[0]), int(mid[1] + text_offset[1])),
                cv2.FONT_HERSHEY_SIMPLEX, ANNOTATION_FONT_SCALE,
                (255, 255, 255), ANNOTATION_THICKNESS, cv2.LINE_AA)


def annotate_via(gray: np.ndarray, rep: ShapeDesc, px_nm: float) -> np.ndarray:
    out = to_color(gray)
    pts = rep.contour.reshape(-1, 2).astype(np.float64)
    x0, x1 = float(pts[:, 0].min()), float(pts[:, 0].max())
    y0, y1 = float(pts[:, 1].min()), float(pts[:, 1].max())
    cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
    cv2.drawContours(out, [rep.contour], -1, (255, 255, 255), 1, cv2.LINE_AA)
    draw_double_arrow(out, np.array([x0, cy]), np.array([x1, cy]), f"X={(x1-x0)*px_nm:.2f} nm")
    draw_double_arrow(out, np.array([cx, y0]), np.array([cx, y1]), f"Y={(y1-y0)*px_nm:.2f} nm", (5, 15))
    return out


def annotate_trench(gray: np.ndarray, pair: dict, px_nm: float) -> np.ndarray:
    out = to_color(gray)
    top = pair["top"]
    bottom = pair["bottom"]
    cv2.drawContours(out, [top.contour, bottom.contour], -1, (255, 255, 255), 1, cv2.LINE_AA)
    x = pair["common_x"]
    p1 = np.array([x, pair["top_tip_y"]], dtype=float)
    p2 = np.array([x, pair["bottom_tip_y"]], dtype=float)
    draw_double_arrow(out, p1, p2, f"Y Tip={pair['gap_px']*px_nm:.2f} nm")
    y = pair["top_width_y"]
    pw1 = np.array([pair["top_width_xl"], y], dtype=float)
    pw2 = np.array([pair["top_width_xr"], y], dtype=float)
    draw_double_arrow(out, pw1, pw2, f"X={pair['top_width_px']*px_nm:.2f} nm", (5, 15))
    for pt in [p1, p2]:
        cv2.circle(out, tuple(np.rint(pt).astype(int)), 3, (255, 255, 255), -1, cv2.LINE_AA)
    return out


def annotate_slot(gray: np.ndarray, reps: Dict[int, ShapeDesc], px_nm: float) -> np.ndarray:
    out = to_color(gray)
    for row_id in [1, 2, 3, 4]:
        if row_id not in reps:
            continue
        d = reps[row_id]
        pts = d.contour.reshape(-1, 2).astype(np.float64)
        x0, x1 = float(pts[:, 0].min()), float(pts[:, 0].max())
        y0, y1 = float(pts[:, 1].min()), float(pts[:, 1].max())
        cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        cv2.drawContours(out, [d.contour], -1, (255, 255, 255), 1, cv2.LINE_AA)
        draw_double_arrow(out, np.array([cx, y0]), np.array([cx, y1]),
                          f"Y{row_id}={(y1-y0)*px_nm:.2f} nm", (5, 0))
        if row_id == 1:
            draw_double_arrow(out, np.array([x0, cy]), np.array([x1, cy]),
                              f"X={(x1-x0)*px_nm:.2f} nm", (5, 15))
    return out


def _draw_green_label(out: np.ndarray, contour: np.ndarray, text: str) -> None:
    green = (0, 255, 0)  # BGR
    x, y, w, h = cv2.boundingRect(contour.astype(np.int32))
    org = (max(0, x), max(12, y - 3))
    cv2.putText(out, text, org, cv2.FONT_HERSHEY_SIMPLEX,
                SELECTED_LABEL_FONT_SCALE, (0, 0, 0), 3, cv2.LINE_AA)
    cv2.putText(out, text, org, cv2.FONT_HERSHEY_SIMPLEX,
                SELECTED_LABEL_FONT_SCALE, green, 1, cv2.LINE_AA)


def annotate_selected_vias(img: np.ndarray, selected: List[ShapeDesc]) -> np.ndarray:
    """Outline every accepted via in GREEN and label A1, A2, ..."""
    if not ANNOTATE_SELECTED_OBJECTS or not selected:
        return img
    out = img.copy()
    green = (0, 255, 0)
    for i, d in enumerate(selected, 1):
        cv2.drawContours(out, [d.contour.astype(np.int32)], -1,
                         green, SELECTED_CONTOUR_THICKNESS, cv2.LINE_AA)
        _draw_green_label(out, d.contour, f"A{i}")
    return out


def annotate_selected_trench_pairs(img: np.ndarray, pairs: List[dict]) -> np.ndarray:
    """Outline every trench component used in a selected upper/lower pair in GREEN."""
    if not ANNOTATE_SELECTED_OBJECTS or not pairs:
        return img
    out = img.copy()
    green = (0, 255, 0)
    for i, p in enumerate(pairs, 1):
        for suffix, d in (("T", p["top"]), ("B", p["bottom"])):
            cv2.drawContours(out, [d.contour.astype(np.int32)], -1,
                             green, SELECTED_CONTOUR_THICKNESS, cv2.LINE_AA)
            _draw_green_label(out, d.contour, f"P{i}{suffix}")
    return out


def annotate_selected_slots(img: np.ndarray, selected: List[ShapeDesc], rows: List[dict]) -> np.ndarray:
    """Outline every accepted slot in GREEN; label includes its Y-cluster row S1..S4."""
    if not ANNOTATE_SELECTED_OBJECTS or not selected:
        return img
    out = img.copy()
    green = (0, 255, 0)
    for i, d in enumerate(selected, 1):
        cv2.drawContours(out, [d.contour.astype(np.int32)], -1,
                         green, SELECTED_CONTOUR_THICKNESS, cv2.LINE_AA)
        row_id = rows[i - 1].get("slot_row", "?") if i - 1 < len(rows) else "?"
        _draw_green_label(out, d.contour, f"A{i}:S{row_id}")
    return out


def accepted_via_descs(descs: List[ShapeDesc]) -> List[ShapeDesc]:
    """Mirror measure_via hard filters for visualization/counting."""
    out: List[ShapeDesc] = []
    for d in descs:
        if d.circularity < VIA_MIN_CIRCULARITY:
            continue
        if d.aspect > VIA_MAX_AXIS_RATIO:
            continue
        if d.solidity < VIA_MIN_SOLIDITY:
            continue
        if min(d.length_px, d.width_px) < VIA_MIN_DIAMETER_PX:
            continue
        out.append(d)
    return out


def _short_reject_code(code: str) -> str:
    mapping = {
        "area_too_small": "small",
        "area_too_large": "large",
        "touches_image_border": "border",
        "descriptor_failed": "desc",
        "via_low_circularity": "circ",
        "via_axis_ratio_too_high": "AR",
        "via_low_solidity": "solid",
        "via_too_small": "size",
        "via_off_array_grid": "off-grid",
        "via_not_selected": "not-sel",
        "trench_aspect_too_low": "AR",
        "trench_too_short": "short",
        "trench_not_vertical": "angle",
        "trench_no_compatible_partner": "no-pair",
        "trench_pair_conflict": "pair-conflict",
        "slot_area_too_small": "small",
        "slot_area_below_0.1x_mean": "area<0.1mean",
        "slot_area_above_10x_mean": "area>10mean",
        "slot_low_regularity_solidity": "concave",
        "slot_low_convexity": "jagged",
        "slot_irregular_template_match": "irregular",
        "slot_aspect_too_low": "AR",
        "slot_too_short": "short",
        "slot_not_vertical": "angle",
        "slot_too_few_valid_objects": "group<4",
        "slot_off_array_grid": "off-grid",
        "slot_not_selected": "not-sel",
    }
    # compound code: show the first reason compactly
    first = code.split("+")[0]
    return mapping.get(first, first[:12])


def annotate_rejections(img: np.ndarray, rejects: List[dict]) -> np.ndarray:
    """Overlay discarded candidates in RED. R-number maps directly to rejected_objects.csv."""
    if not ANNOTATE_REJECTED_OBJECTS or not rejects:
        return img
    out = img.copy()
    red = (0, 0, 255)  # OpenCV BGR
    shown = 0
    for r in rejects:
        if shown >= MAX_REJECTED_ANNOTATIONS_PER_IMAGE:
            break
        if r.get("stage") == "segmentation_qc" and not ANNOTATE_SEGMENTATION_REJECTS:
            continue
        c = r.get("_contour")
        if c is None or len(c) == 0:
            continue
        shown += 1
        reject_id = r.get("reject_id", shown)
        cv2.drawContours(out, [c.astype(np.int32)], -1, red, REJECTED_CONTOUR_THICKNESS, cv2.LINE_AA)
        x, y, w, h = cv2.boundingRect(c.astype(np.int32))
        cv2.line(out, (x, y), (x + w - 1, y + h - 1), red, REJECTED_X_THICKNESS, cv2.LINE_AA)
        cv2.line(out, (x + w - 1, y), (x, y + h - 1), red, REJECTED_X_THICKNESS, cv2.LINE_AA)
        label = f"R{reject_id}:{_short_reject_code(str(r.get('reason_code', 'reject')))}"
        org = (max(0, x), max(12, y - 3))
        # black shadow then red label for readability
        cv2.putText(out, label, org, cv2.FONT_HERSHEY_SIMPLEX,
                    REJECTED_LABEL_FONT_SCALE, (0, 0, 0), 3, cv2.LINE_AA)
        cv2.putText(out, label, org, cv2.FONT_HERSHEY_SIMPLEX,
                    REJECTED_LABEL_FONT_SCALE, red, 1, cv2.LINE_AA)
    return out


def format_rejection_summary(rejects: List[dict]) -> str:
    if not rejects:
        return ""
    counts = Counter(str(r.get("reason_code", "unknown")) for r in rejects)
    return "; ".join(f"{k}: {v}" for k, v in counts.most_common())


def exportable_rejection_record(r: dict) -> dict:
    return {k: v for k, v in r.items() if not k.startswith("_")}


# ============================================================
# 9) IMAGE-LEVEL AND CONDITION-LEVEL SUMMARIES
# ============================================================
def build_image_summary(objects_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if objects_df.empty:
        return pd.DataFrame()

    for (condition, image, pattern), g in objects_df.groupby(["condition", "image", "pattern"], dropna=False):
        base = {
            "condition": str(condition),
            "image": image,
            "pattern": pattern,
            "n_objects": len(g),
            "pixel_size_nm": float(g["pixel_size_nm"].iloc[0]),
        }
        if pattern == "via":
            base.update({
                "via_width_nm": g["width_nm"].mean(),
                "via_height_nm": g["height_nm"].mean(),
                "via_height_width_ratio": g["height_width_ratio"].mean(),
                "via_width_sd_nm": g["width_nm"].std(ddof=1) if len(g) > 1 else np.nan,
                "via_height_sd_nm": g["height_nm"].std(ddof=1) if len(g) > 1 else np.nan,
            })
        elif pattern == "trench":
            base.update({
                "trench_width_nm": g["trench_width_nm"].mean(),
                "tip_to_tip_nm": g["tip_to_tip_nm"].mean(),
                "trench_width_sd_nm": g["trench_width_nm"].std(ddof=1) if len(g) > 1 else np.nan,
                "tip_to_tip_sd_nm": g["tip_to_tip_nm"].std(ddof=1) if len(g) > 1 else np.nan,
            })
        elif pattern == "slot":
            base["slot_width_nm"] = g["slot_width_nm"].mean()
            base["slot_width_sd_nm"] = g["slot_width_nm"].std(ddof=1) if len(g) > 1 else np.nan
            for r in [1, 2, 3, 4]:
                gr = g[g["slot_row"] == r]
                base[f"slot_L{r}_nm"] = gr["slot_length_nm"].mean() if len(gr) else np.nan
                base[f"slot_L{r}_n"] = len(gr)
        rows.append(base)

    return pd.DataFrame(rows)


def build_condition_summary(image_summary: pd.DataFrame) -> pd.DataFrame:
    """Equal-weight images: average image means within each condition/pattern."""
    if image_summary.empty:
        return pd.DataFrame()

    rows = []
    numeric_cols = [c for c in image_summary.columns if c not in {"condition", "image", "pattern"}]
    for (condition, pattern), g in image_summary.groupby(["condition", "pattern"], dropna=False):
        row = {
            "condition": str(condition),
            "pattern": pattern,
            "n_images": len(g),
            "n_objects_total": int(g["n_objects"].sum()),
        }
        for c in numeric_cols:
            if c in {"n_objects", "pixel_size_nm"} or c.endswith("_n") or c.endswith("_sd_nm"):
                continue
            vals = pd.to_numeric(g[c], errors="coerce").dropna() if c in g.columns else pd.Series(dtype=float)
            if len(vals):
                row[c] = vals.mean()
                row[c + "__sd_between_images"] = vals.std(ddof=1) if len(vals) > 1 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


# ============================================================
# 10) PLOTTING
# ============================================================
def _configure_chinese_matplotlib() -> None:
    """Use an installed Chinese font, preferring Windows Microsoft YaHei."""
    try:
        from matplotlib import font_manager
        font_paths = [
            Path(r"C:\Windows\Fonts\msyh.ttc"),
            Path(r"C:\Windows\Fonts\msyhbd.ttc"),
            Path(r"C:\Windows\Fonts\simhei.ttf"),
            Path(r"C:\Windows\Fonts\simsun.ttc"),
        ]
        chosen = None
        for fp in font_paths:
            if fp.exists():
                try:
                    font_manager.fontManager.addfont(str(fp))
                    chosen = font_manager.FontProperties(fname=str(fp)).get_name()
                    break
                except Exception:
                    pass
        if chosen:
            plt.rcParams["font.sans-serif"] = [chosen, "Microsoft YaHei", "SimHei", "DejaVu Sans"]
        else:
            plt.rcParams["font.sans-serif"] = [
                "Microsoft YaHei", "SimHei", "Microsoft JhengHei", "Noto Sans CJK SC",
                "Arial Unicode MS", "DejaVu Sans",
            ]
    except Exception:
        plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False


def _condition_tick_labels_zh(sub: pd.DataFrame) -> List[str]:
    return [CONDITION_LABELS_ZH.get(str(x), str(x)) for x in sub["condition"].astype(str)]


def save_bar_with_error(df: pd.DataFrame, pattern: str, metrics: List[Tuple[str, str]], out_path: Path, title: str) -> None:
    _configure_chinese_matplotlib()
    sub = df[df["pattern"] == pattern].copy()
    if sub.empty:
        return
    sub["condition"] = pd.Categorical(sub["condition"], categories=CONDITIONS, ordered=True)
    sub = sub.sort_values("condition")

    x = np.arange(len(sub))
    width = 0.8 / max(len(metrics), 1)
    fig, ax = plt.subplots(figsize=(10.5, 5.5))
    for i, (col, label) in enumerate(metrics):
        if col not in sub.columns:
            continue
        y = pd.to_numeric(sub[col], errors="coerce").to_numpy(dtype=float)
        sd_col = col + "__sd_between_images"
        yerr = pd.to_numeric(sub[sd_col], errors="coerce").to_numpy(dtype=float) if sd_col in sub.columns else None
        xpos = x - 0.4 + width / 2 + i * width
        ax.bar(xpos, y, width=width, label=label, yerr=yerr, capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels(_condition_tick_labels_zh(sub), rotation=10, ha="right", fontsize=9)
    ax.set_xlabel("刻蚀条件")
    ax.set_ylabel("长度（nm）" if any("_nm" in c for c, _ in metrics) else "数值")
    ax.set_title(title)
    ax.legend()
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def save_via_ratio_plot(df: pd.DataFrame, out_path: Path) -> None:
    _configure_chinese_matplotlib()
    sub = df[df["pattern"] == "via"].copy()
    if sub.empty or "via_height_width_ratio" not in sub.columns:
        return
    sub["condition"] = pd.Categorical(sub["condition"], categories=CONDITIONS, ordered=True)
    sub = sub.sort_values("condition")
    x = np.arange(len(sub))
    y = sub["via_height_width_ratio"].to_numpy(float)
    sd_col = "via_height_width_ratio__sd_between_images"
    yerr = sub[sd_col].to_numpy(float) if sd_col in sub.columns else None
    fig, ax = plt.subplots(figsize=(10.5, 5.5))
    ax.bar(x, y, yerr=yerr, capsize=3)
    ax.axhline(1.0, linestyle="--", linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(_condition_tick_labels_zh(sub), rotation=10, ha="right", fontsize=9)
    ax.set_xlabel("刻蚀条件")
    ax.set_ylabel("y方向长度 / x方向长度")
    ax.set_title("圆孔 y/x 方向长度比对比")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def save_slot_length_plot(df: pd.DataFrame, out_path: Path) -> None:
    _configure_chinese_matplotlib()
    sub = df[df["pattern"] == "slot"].copy()
    if sub.empty:
        return
    sub["condition"] = pd.Categorical(sub["condition"], categories=CONDITIONS, ordered=True)
    sub = sub.sort_values("condition")
    x = np.arange(len(sub))
    fig, ax = plt.subplots(figsize=(10.5, 5.5))
    for r in [1, 2, 3, 4]:
        col = f"slot_L{r}_nm"
        if col not in sub.columns:
            continue
        y = sub[col].to_numpy(float)
        ax.plot(x, y, marker="o", label=f"短槽{r}")
    ax.set_xticks(x)
    ax.set_xticklabels(_condition_tick_labels_zh(sub), rotation=10, ha="right", fontsize=9)
    ax.set_xlabel("刻蚀条件")
    ax.set_ylabel("y方向长度（nm）")
    ax.set_title("不同短槽编号的y方向长度对比（从上到下1–4）")
    ax.legend(title="短槽编号")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def make_all_plots(condition_summary: pd.DataFrame, plot_dir: Path) -> None:
    plot_dir.mkdir(parents=True, exist_ok=True)
    save_bar_with_error(
        condition_summary,
        "via",
        [("via_width_nm", "x方向长度"), ("via_height_nm", "y方向长度")],
        plot_dir / "via_width_height.png",
        "圆孔 x/y 方向长度对比",
    )
    save_via_ratio_plot(condition_summary, plot_dir / "via_height_width_ratio.png")
    save_bar_with_error(
        condition_summary,
        "trench",
        [("trench_width_nm", "沟槽x方向宽度"), ("tip_to_tip_nm", "尖端间y方向距离")],
        plot_dir / "trench_width_tipgap.png",
        "沟槽x方向宽度与尖端间y方向距离对比",
    )
    save_bar_with_error(
        condition_summary,
        "slot",
        [("slot_width_nm", "短槽平均x方向长度")],
        plot_dir / "slot_width.png",
        "短槽平均x方向长度对比",
    )
    save_slot_length_plot(condition_summary, plot_dir / "slot_lengths_1_to_4.png")


# ============================================================
# 11) MAIN BATCH PIPELINE
# ============================================================
def find_tifs(condition_dir: Path) -> List[Path]:
    files = []
    for ext in ("*.tif", "*.tiff", "*.TIF", "*.TIFF"):
        files.extend(condition_dir.rglob(ext))
    return sorted(set(files))



# ============================================================
# V18: DIAMOND VIA + DEEP 4-ROW SHORT-SLOT RECOVERY
# ============================================================
# VIA rule supplied by user:
#   * filename stem contains "-1" -> diamond / rhombic array whose nearest-neighbor
#     center-to-center vectors are approximately +/-45 degrees to horizontal.
#     For these files DO NOT reject objects for missing H/V neighbors.
#   * other VIA files keep the stable V16 horizontal/vertical fixed-pitch logic.
V18_VIA_DIAG_ANGLE_TOL_DEG = 14.0
V18_VIA_DIAG_COMPONENT_REL_TOL = 0.22
V18_VIA_DIAG_SITE_TOL_FRAC = 0.25
V18_VIA_DIAG_OCCUPIED_TOL_FRAC = 0.20
V18_VIA_DIAG_MAX_ATTEMPTS = 600
V18_VIA_DIAG_MAX_RECOVERED = 260

# SLOT rule supplied by user:
#   * exactly four valid rows are a hard layout prior;
#   * area and geometry QC are NEVER loosened during recovery;
#   * only contrast/blur segmentation is progressively made more sensitive;
#   * broad low-gray regions may be used as POSITION proposals, but they cannot be
#     accepted until a complete regular contour passes the unchanged shape/area QC.
V18_SLOT_DEEP_RECOVERY_PASSES = 4
V18_SLOT_TARGET_APPROX_PER_ROW = 10
V18_SLOT_LOWGRAY_X_SEARCH_FRAC = 0.46
V18_SLOT_LOWGRAY_Y_SEARCH_FRAC = 0.44
V18_SLOT_LOWGRAY_BOX_W_FACTOR = 0.80
V18_SLOT_LOWGRAY_BOX_H_FACTOR = 0.72
V18_SLOT_LOWGRAY_GAUSS_SIGMAS = (0.0, 1.2, 2.2, 3.5, 5.0)
V18_SLOT_DIFF_SIGMAS = (5.0, 9.0, 14.0, 22.0)
V18_SLOT_DARK_PERCENTILES = (4.0, 7.0, 10.0, 13.0, 16.0, 20.0, 24.0, 28.0, 32.0, 36.0, 40.0, 45.0, 50.0, 55.0)
V18_SLOT_BRIGHT_PERCENTILES = tuple(100.0 - q for q in V18_SLOT_DARK_PERCENTILES)
V18_SLOT_INTERPOLATION_ALPHAS = (0.18, 0.28, 0.38, 0.50, 0.62, 0.74, 0.84)
V18_SLOT_MAX_LOCAL_MASKS = 70
V18_SLOT_MODEL_LOWGRAY_MAX_PROPOSALS = 180
V18_SLOT_LOWGRAY_MODEL_AREA_MIN = 20
V18_SLOT_LOWGRAY_MODEL_AREA_MAX_FRAC = 0.10
V18_SLOT_ROW_SEARCH_UNTIL_NO_PROGRESS = True
V18_SLOT_MAX_NO_PROGRESS_PASSES = 2


def _image_allows_diagonal_v16(image_name: str) -> bool:
    """V18 fix: '-1' means a standalone dash-delimited token, not the prefix of '-160'.

    Examples:
      ZEP210-160-1.tif -> True
      ZEP210-160-2.tif -> False
      ZEP210-120-2.tif -> False
    """
    parts = [x.strip() for x in Path(image_name).stem.split('-')]
    return any(x == '1' for x in parts)


def _via_diag_component_pitch_v18(items: List[ShapeDesc], image_shape: Tuple[int, int]) -> float:
    """One frozen diagonal x/y component from center pairs close to +/-45 degrees."""
    if len(items) < 2:
        return np.nan
    P = np.asarray([d.center for d in items], dtype=float)
    vals = []
    max_allowed = max(image_shape) * V16_PITCH_MAX_IMAGE_FRAC
    tan_lo = math.tan(math.radians(45.0 - V18_VIA_DIAG_ANGLE_TOL_DEG))
    tan_hi = math.tan(math.radians(45.0 + V18_VIA_DIAG_ANGLE_TOL_DEG))
    for i in range(len(P)):
        near = []
        for j in range(len(P)):
            if i == j:
                continue
            dx = abs(float(P[j,0] - P[i,0])); dy = abs(float(P[j,1] - P[i,1]))
            if dx < V16_PITCH_MIN_PX / 1.414 or dy < V16_PITCH_MIN_PX / 1.414:
                continue
            dist = math.hypot(dx, dy)
            if dist > max_allowed:
                continue
            ratio = dy / max(dx, 1e-9)
            if not (tan_lo <= ratio <= tan_hi):
                continue
            comp = 0.5 * (dx + dy)
            near.append(comp)
        if near:
            vals.append(min(near))
    if not vals:
        return np.nan
    return _robust_one_pitch_v16(vals)


def _largest_diamond_component_v18(items: List[ShapeDesc], comp: float) -> List[ShapeDesc]:
    """Largest center-connected component using ONLY the four diagonal neighbor vectors."""
    if len(items) < 3 or not np.isfinite(comp):
        return list(items)
    targets = [
        np.array([ comp, comp]), np.array([ comp,-comp]),
        np.array([-comp, comp]), np.array([-comp,-comp]),
    ]
    tol = max(4.0, V18_VIA_DIAG_SITE_TOL_FRAC * comp)
    P = np.asarray([d.center for d in items], dtype=float)
    adj=[[] for _ in items]
    for i in range(len(items)):
        for j in range(i+1,len(items)):
            dv=P[j]-P[i]
            if min(float(np.linalg.norm(dv-t)) for t in targets) <= tol:
                adj[i].append(j); adj[j].append(i)
    seen=set(); comps=[]
    for i in range(len(items)):
        if i in seen: continue
        st=[i]; seen.add(i); cc=[]
        while st:
            k=st.pop(); cc.append(k)
            for q in adj[k]:
                if q not in seen:
                    seen.add(q); st.append(q)
        comps.append(cc)
    comps.sort(key=lambda c:(-len(c),min(c)))
    if not comps or len(comps[0]) < 2:
        return list(items)
    return [items[i] for i in comps[0]]


def recursive_diamond_via_recovery_v18(
    gray: np.ndarray,
    pool: List[ShapeDesc],
    initial_selected: List[ShapeDesc],
    image_name: str,
) -> Tuple[List[ShapeDesc], Dict[str, object]]:
    """
    '-1' VIA images: four diagonal directions ONLY.

    The center pitch is estimated from the diagonal nearest-neighbor pairs once and is
    frozen. Horizontal/vertical connectivity is intentionally not part of acceptance.
    """
    meta = {
        'array_assist_used': False,
        'array_lattice_type': 'diamond_45_center_array',
        'array_expected_positions': 0,
        'array_recovered_count': 0,
        'array_x_pitch_px': np.nan, 'array_y_pitch_px': np.nan,
        'array_occupancy_seed': np.nan,
        'array_initial_recovered_count': 0, 'array_final_recovered_count': 0,
        'array_final_predicted_positions': 0,
        'array_recovery_attempts': 0, 'array_final_recovery_attempts': 0,
        'array_recovery_skipped_due_cap': 0,
        'array_basis_a_dx_px': np.nan, 'array_basis_a_dy_px': np.nan,
        'array_basis_b_dx_px': np.nan, 'array_basis_b_dy_px': np.nan,
        'array_basis_a_length_px': np.nan, 'array_basis_b_length_px': np.nan,
        'array_basis_angle_deg': 90.0,
        'center_recursive_used': False,
        'center_horizontal_pitch_px': np.nan, 'center_vertical_pitch_px': np.nan,
        'center_diagonal_component_px': np.nan,
        'center_recursive_attempts': 0, 'center_recursive_recovered': 0,
        'center_recursive_boundary_hits': 0,
        'fixed_diagonal_enabled': True, 'fixed_pitch_frozen': True,
    }
    selected=list(initial_selected)
    if len(selected)<2:
        return selected,meta
    comp=_via_diag_component_pitch_v18(initial_selected,gray.shape)
    if not np.isfinite(comp) or comp < V16_PITCH_MIN_PX/1.414:
        # If strict seeds are sparse, derive the diagonal component from every complete
        # regular via candidate, but final acceptance remains strict.
        relaxed=[d for d in pool if _complete_regular_shape_ok_v16(d,'via',np.nan,None)]
        relaxed=_dedup_desc_list_v16(_seed_area_consistency_v16(relaxed))
        comp=_via_diag_component_pitch_v18(relaxed,gray.shape)
        if np.isfinite(comp):
            selected=relaxed
    if not np.isfinite(comp):
        return selected,meta
    selected=_largest_diamond_component_v18(selected,comp)
    if len(selected)<2:
        return selected,meta

    nn=comp*V16_DIAGONAL_DIVISOR
    meta.update({
        'array_assist_used':True,
        # For a diamond lattice, these are the full horizontal/vertical repeat distances
        # between every second diagonal neighbor. The true nearest-neighbor component is
        # reported separately below.
        'array_x_pitch_px':float(2.0*comp), 'array_y_pitch_px':float(2.0*comp),
        'center_horizontal_pitch_px':float(2.0*comp), 'center_vertical_pitch_px':float(2.0*comp),
        'center_diagonal_component_px':float(comp),
        'center_recursive_used':True,
        'array_basis_a_dx_px':float(comp), 'array_basis_a_dy_px':float(comp),
        'array_basis_b_dx_px':float(comp), 'array_basis_b_dy_px':float(-comp),
        'array_basis_a_length_px':float(nn), 'array_basis_b_length_px':float(nn),
        'array_occupancy_seed':float(len(selected)/max(1,len(initial_selected))),
    })
    H,W=gray.shape
    widths=[max(1.0,_contour_xy_bounds(d)[1]-_contour_xy_bounds(d)[0]) for d in selected]
    heights=[max(1.0,_contour_xy_bounds(d)[3]-_contour_xy_bounds(d)[2]) for d in selected]
    areas=[contour_pixel_count(d.contour) for d in selected]
    med_w=float(np.median(widths)); med_h=float(np.median(heights)); med_area=float(np.median(areas))
    tol=max(4.0,V18_VIA_DIAG_SITE_TOL_FRAC*comp)
    occ_tol=max(4.0,V18_VIA_DIAG_OCCUPIED_TOL_FRAC*comp)
    directions=[np.array([ comp, comp]),np.array([ comp,-comp]),np.array([-comp, comp]),np.array([-comp,-comp])]
    from collections import deque
    q=deque(selected); visited=set(); attempts=recovered=boundary_hits=0
    visit_cell=max(3.0,0.55*tol)
    def occupied(xy):
        return any(float(np.linalg.norm(d.center-xy))<=occ_tol for d in selected)
    while q and attempts<V18_VIA_DIAG_MAX_ATTEMPTS and recovered<V18_VIA_DIAG_MAX_RECOVERED:
        base=q.popleft()
        for dv in directions:
            ex=np.asarray(base.center,dtype=float)+dv
            x,y=map(float,ex)
            mx=max(V16_BOUNDARY_MARGIN_PX,0.55*med_w); my=max(V16_BOUNDARY_MARGIN_PX,0.55*med_h)
            if x<=mx or x>=W-1-mx or y<=my or y>=H-1-my:
                boundary_hits+=1; continue
            key=(int(round(x/visit_cell)),int(round(y/visit_cell)))
            if key in visited: continue
            visited.add(key)
            if occupied(ex): continue
            attempts+=1
            if attempts==1 or attempts%V16_RECURSIVE_PROGRESS_EVERY==0:
                print(f"      [{image_name}] V18 45度菱形via补点 {attempts}/{V18_VIA_DIAG_MAX_ATTEMPTS}; recovered={recovered} ...",flush=True)
            d=_select_pool_candidate_fixed_site_v16(pool,selected,(x,y),1.15*tol,'via',med_area,None)
            if d is None:
                # fixed_pitch only limits ROI size here; use true diagonal neighbor distance.
                d=recover_local_object_v16_fixed(gray,(x,y),med_w,med_h,'via',nn,med_area,None)
                if d is not None: pool.append(d)
            if d is None or _desc_close_to_any_selected_v12(d,selected):
                continue
            if _source_of(d)!='fixed_center_recovered':
                _mark_source(d,'diamond_45_selected')
            else:
                _mark_source(d,'diamond_45_recovered')
            selected.append(d); q.append(d); recovered+=1
    meta.update({
        'array_recovered_count':int(recovered),'array_final_recovered_count':int(recovered),
        'array_final_predicted_positions':int(len(visited)),
        'array_recovery_attempts':int(attempts),'array_final_recovery_attempts':int(attempts),
        'center_recursive_attempts':int(attempts),'center_recursive_recovered':int(recovered),
        'center_recursive_boundary_hits':int(boundary_hits),
    })
    return selected,meta


def measure_via_v18(
    gray: np.ndarray, pool: List[ShapeDesc], px_nm: float, image_name: str, condition: str
) -> Tuple[List[dict], Optional[ShapeDesc], List[ShapeDesc], List[dict], Dict[str, object]]:
    """Normal VIA -> V16 H/V. '-1' VIA -> true 45-degree diamond neighbor array."""
    if not _image_allows_diagonal_v16(image_name):
        return measure_via_v16(gray,pool,px_nm,image_name,condition)
    strict=_dedup_desc_list_v16([d for d in pool if _via_strict_metrics_v12(d)[0]])
    strict=_seed_area_consistency_v16(strict)
    if len(strict)<2:
        fallback=[d for d in pool if _complete_regular_shape_ok_v16(d,'via',np.nan,None)]
        strict=_dedup_desc_list_v16(_seed_area_consistency_v16(fallback))
    selected,meta=recursive_diamond_via_recovery_v18(gray,pool,strict,image_name)
    rows=[]
    for oid,d in enumerate(sorted(selected,key=lambda z:(z.center[1],z.center[0])),1):
        pts=d.contour.reshape(-1,2).astype(np.float64)
        w=float(pts[:,0].max()-pts[:,0].min()); h=float(pts[:,1].max()-pts[:,1].min())
        rows.append({
            'condition':condition,'image':image_name,'pattern':'via','object_id':oid,
            'width_nm':w*px_nm,'height_nm':h*px_nm,'height_width_ratio':h/max(w,1e-12),
            'circularity':d.circularity,'solidity':d.solidity,'axis_ratio':d.aspect,
            'selection_source':_source_of(d),'pixel_size_nm':px_nm,
        })
    rejected=[]
    for d in pool:
        if _desc_close_to_any_selected_v12(d,selected): continue
        ok,codes,reasons=_via_strict_metrics_v12(d)
        if _source_of(d)!='primary' and not _complete_regular_shape_ok_v16(d,'via',np.nan,None):
            continue
        if ok and not codes:
            codes=['via_not_on_45deg_diamond_array']; reasons=['regular via candidate center was not consistent with the frozen +/-45 degree diamond nearest-neighbor array']
        elif not codes:
            codes=['via_not_selected']; reasons=['candidate was not selected']
        rejected.append(_make_reject_from_contour(
            d.contour,'via_v18_diamond_filter','+'.join(codes),'; '.join(reasons),desc=d,
            source=_source_of(d),circularity=d.circularity,solidity=d.solidity,axis_ratio=d.aspect,
        ))
    representative=max(selected,key=lambda d:(d.circularity,d.solidity)) if selected else None
    return rows,representative,selected,rejected,meta


def _slot_lowgray_model_proposals_v18(gray: np.ndarray, image_name: str = "") -> List[ShapeDesc]:
    """Very permissive LOW-GRAY localization proposals for ARRAY FITTING ONLY."""
    H,W=gray.shape
    g=gray.astype(np.uint8)
    out=[]
    max_area=V18_SLOT_LOWGRAY_MODEL_AREA_MAX_FRAC*H*W
    # Blur suppresses SEM high-frequency noise while preserving broad low-gray objects.
    for sigma in (1.2,2.5,4.0):
        if image_name:
            print(f"      [{image_name}] slot low-gray localization: sigma={sigma}, proposals={len(out)}", flush=True)
        b=cv2.GaussianBlur(g,(0,0),sigmaX=sigma,sigmaY=sigma) if sigma>0 else g
        qs=(12.0,18.0,24.0,30.0,36.0,42.0) if FEATURE_POLARITY.lower()=='dark' else (58.0,64.0,70.0,76.0,82.0,88.0)
        for qv in qs:
            th=float(np.percentile(b,qv))
            if FEATURE_POLARITY.lower()=='dark': mm=(b<=th).astype(np.uint8)*255
            else: mm=(b>=th).astype(np.uint8)*255
            mm=cv2.morphologyEx(mm,cv2.MORPH_CLOSE,np.ones((3,3),np.uint8),iterations=1)
            # Consume each full-image mask immediately; retaining all 18 masks
            # wastes memory and prevents the proposal cap from stopping generation.
            contours,_=cv2.findContours(mm,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
            for c in contours:
                pa=contour_pixel_count(c)
                if pa<V18_SLOT_LOWGRAY_MODEL_AREA_MIN or pa>max_area: continue
                x,y,w,h=cv2.boundingRect(c)
                if x<=1 or y<=1 or x+w>=W-2 or y+h>=H-2: continue
                d=describe_contour(c)
                if d is None: continue
                # Localization-only dedup; NO shape acceptance here.
                if not _desc_close_to_any_selected_v12(d,out):
                    _mark_source(d,'slot_lowgray_locator_only')
                    out.append(d)
                    if len(out)>=V18_SLOT_MODEL_LOWGRAY_MAX_PROPOSALS:
                        return out
    return out


def _slot_darkness_center_v18(
    gray: np.ndarray, expected_xy: Tuple[float,float], x_pitch: float, row_pitch: float,
    expected_w: float, expected_h: float
) -> Tuple[Tuple[float,float], float]:
    """Locate a broad dark/bright basin near a predicted site using local box-average intensity."""
    H,W=gray.shape; cx,cy=map(float,expected_xy)
    hw=max(10,int(round(V18_SLOT_LOWGRAY_X_SEARCH_FRAC*x_pitch)))
    hh=max(10,int(round(V18_SLOT_LOWGRAY_Y_SEARCH_FRAC*row_pitch)))
    x0=max(0,int(round(cx-hw))); x1=min(W,int(round(cx+hw+1)))
    y0=max(0,int(round(cy-hh))); y1=min(H,int(round(cy+hh+1)))
    if x1-x0<7 or y1-y0<7: return (cx,cy),np.nan
    roi=gray[y0:y1,x0:x1].astype(np.float32)
    kw=max(3,int(round(max(3.0,expected_w*V18_SLOT_LOWGRAY_BOX_W_FACTOR))))
    kh=max(3,int(round(max(3.0,expected_h*V18_SLOT_LOWGRAY_BOX_H_FACTOR))))
    kw=min(kw,max(3,roi.shape[1]//2*2-1)); kh=min(kh,max(3,roi.shape[0]//2*2-1))
    if kw%2==0: kw=max(3,kw-1)
    if kh%2==0: kh=max(3,kh-1)
    sm=cv2.GaussianBlur(roi,(0,0),sigmaX=1.4,sigmaY=1.4)
    mean_map=cv2.boxFilter(sm,-1,(kw,kh),normalize=True,borderType=cv2.BORDER_REFLECT)
    # Prevent selecting window centers that would make the slot cross ROI boundary.
    mx=max(1,kw//2); my=max(1,kh//2)
    work=mean_map.copy()
    if FEATURE_POLARITY.lower()=='dark':
        work[:my,:]=np.inf; work[-my:,:]=np.inf; work[:,:mx]=np.inf; work[:,-mx:]=np.inf
        iy,ix=np.unravel_index(int(np.argmin(work)),work.shape); score=float(work[iy,ix])
    else:
        work[:my,:]=-np.inf; work[-my:,:]=-np.inf; work[:,:mx]=-np.inf; work[:,-mx:]=-np.inf
        iy,ix=np.unravel_index(int(np.argmax(work)),work.shape); score=float(work[iy,ix])
    return (float(x0+ix),float(y0+iy)),score


def _v18_slot_masks_at_roi(roi: np.ndarray, aggressive_level: int) -> Dict[str,np.ndarray]:
    """Progressively blur/contrast-tolerant masks. Final geometry QC remains unchanged."""
    masks: Dict[str,np.ndarray]={}
    u=roi.astype(np.uint8)
    # Standard methods first.
    try: masks['primary']=segment_features(u)
    except Exception: pass
    try:
        clahe=cv2.createCLAHE(clipLimit=2.0+0.8*aggressive_level,tileGridSize=(8,8)).apply(u)
        masks['clahe_otsu']=segment_features(clahe)
    except Exception: pass
    # Smoothed percentile masks. Stronger levels add more smoothing and more percentiles.
    sigmas=V18_SLOT_LOWGRAY_GAUSS_SIGMAS[:min(len(V18_SLOT_LOWGRAY_GAUSS_SIGMAS),2+aggressive_level)]
    percs=V18_SLOT_DARK_PERCENTILES if FEATURE_POLARITY.lower()=='dark' else V18_SLOT_BRIGHT_PERCENTILES
    step=max(1,4-aggressive_level)
    percs=percs[::step]
    for sigma in sigmas:
        b=cv2.GaussianBlur(u,(0,0),sigmaX=sigma,sigmaY=sigma) if sigma>0 else u
        for qv in percs:
            th=float(np.percentile(b,qv))
            if FEATURE_POLARITY.lower()=='dark': mm=(b<=th).astype(np.uint8)*255
            else: mm=(b>=th).astype(np.uint8)*255
            # Close blur-broken edges, then a very light open to kill isolated noise.
            mm=cv2.morphologyEx(mm,cv2.MORPH_CLOSE,np.ones((3,3),np.uint8),iterations=1+(aggressive_level>=3))
            if aggressive_level>=2:
                mm=cv2.morphologyEx(mm,cv2.MORPH_OPEN,np.ones((2,2),np.uint8),iterations=1)
            masks[f'pct_s{sigma:g}_{qv:g}']=mm
            if len(masks)>=V18_SLOT_MAX_LOCAL_MASKS: return masks
    # Background-subtraction emphasizes a large low-gray basin even if global contrast is weak.
    for sigma in V18_SLOT_DIFF_SIGMAS[:min(len(V18_SLOT_DIFF_SIGMAS),1+aggressive_level)]:
        bg=cv2.GaussianBlur(u.astype(np.float32),(0,0),sigmaX=sigma,sigmaY=sigma)
        diff=(bg-u.astype(np.float32)) if FEATURE_POLARITY.lower()=='dark' else (u.astype(np.float32)-bg)
        for qv in (55,65,72,78,84,89,93):
            th=float(np.percentile(diff,qv))
            if th<=0: continue
            mm=(diff>=th).astype(np.uint8)*255
            mm=cv2.morphologyEx(mm,cv2.MORPH_CLOSE,np.ones((3,3),np.uint8),iterations=1)
            masks[f'diff_{sigma:g}_{qv}']=mm
            if len(masks)>=V18_SLOT_MAX_LOCAL_MASKS: return masks
    return masks


_V32_SLOT_RECOVERY_DIAGNOSTICS = Counter()


def _v32_slot_flatfield_response(roi: np.ndarray) -> Tuple[np.ndarray, float]:
    """Estimate background from the ROI rim, excluding foreground outliers.

    Subtracting a fitted plane preserves feature contrast under a brightness
    gradient; stretching the raw ROI alone cannot remove that gradient.
    """
    smooth = cv2.GaussianBlur(roi.astype(np.float32), (0, 0), 0.9)
    yy, xx = np.mgrid[-1:1:complex(roi.shape[0]), -1:1:complex(roi.shape[1])]
    rim = (np.abs(xx) >= 0.70) | (np.abs(yy) >= 0.85)
    design = np.column_stack([np.ones(np.count_nonzero(rim)), xx[rim], yy[rim]])
    values = smooth[rim].astype(float)
    stride = max(1, len(values) // 6000)
    design, values = design[::stride], values[::stride]
    keep = np.ones(len(values), dtype=bool)
    sign = 1.0 if FEATURE_POLARITY.lower() == 'dark' else -1.0
    for _ in range(3):
        if np.count_nonzero(keep) < 8:
            break
        coef, *_ = np.linalg.lstsq(design[keep], values[keep], rcond=None)
        residual = sign * (values - design @ coef)
        median = float(np.median(residual))
        spread = max(0.5, 1.4826 * float(np.median(np.abs(residual - median))))
        keep = (residual >= median - 2.5 * spread) & (residual <= median + 3.5 * spread)
    background = coef[0] + coef[1] * xx + coef[2] * yy
    noise_residual = roi.astype(float) - smooth
    noise = 1.4826 * float(np.median(np.abs(noise_residual - np.median(noise_residual))))
    return np.maximum(sign * (background - smooth), 0).astype(np.float32), noise


def _v32_slot_flatfield_candidates(gray: np.ndarray, pool: List[ShapeDesc]) -> List[ShapeDesc]:
    """Find real contours before lattice fitting, even with uneven illumination."""
    h, w = gray.shape
    smooth = cv2.GaussianBlur(gray, (0, 0), 0.9)
    residual = gray.astype(float) - smooth.astype(float)
    noise = 1.4826 * float(np.median(np.abs(residual - np.median(residual))))
    widths = [cv2.boundingRect(d.contour)[2] for d in pool
              if _slot_shape_ok_v17(d, (SLOT_MIN_AREA_PX, np.inf), seed=True)]
    width = float(np.median(widths)) if widths else max(5.0, w / 32.0)
    kernels = sorted(set(max(5, min(129, int(round(width * factor)) | 1)) for factor in (2., 3., 4.)))
    result = []
    for size in kernels:
        kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (size, 1))
        operation = cv2.MORPH_BLACKHAT if FEATURE_POLARITY.lower() == 'dark' else cv2.MORPH_TOPHAT
        response = cv2.morphologyEx(smooth, operation, kernel).astype(np.float32)
        peak = float(np.percentile(response, 99.5))
        if peak < max(2.0, 3.0 * noise):
            continue
        for fraction in (0.50, 0.35):
            threshold = max(1.5, 2.5 * noise, fraction * peak)
            mask = (response >= threshold).astype(np.uint8) * 255
            mask = _postprocess_binary_v12(mask, 1)
            candidates = extract_loose_descriptors_v12(mask, 'slot_v32_flatfield', 'slot', 160)
            for d in candidates:
                if not _slot_shape_ok_v17(d, (SLOT_MIN_AREA_PX, V18_SLOT_LOWGRAY_MODEL_AREA_MAX_FRAC * h * w), seed=True):
                    continue
                if not _desc_close_to_any_selected_v12(d, result):
                    result.append(d)
                    if len(result) >= 160:
                        return result
    return result


def _v32_recover_flatfield_slot(
    gray, expected_xy, x_pitch, row_pitch, expected_w, expected_h, area_limits,
) -> Optional[ShapeDesc]:
    h, w = gray.shape
    cx, cy = map(float, expected_xy)
    # Use the site cell, not a narrow ROI derived from a possibly fragmented seed.
    hw = max(14., 0.48 * x_pitch)
    hh = max(16., 0.72 * row_pitch)
    x0, x1 = max(0, int(cx - hw)), min(w, int(math.ceil(cx + hw + 1)))
    y0, y1 = max(0, int(cy - hh)), min(h, int(math.ceil(cy + hh + 1)))
    roi = gray[y0:y1, x0:x1]
    if min(roi.shape, default=0) < 9:
        _V32_SLOT_RECOVERY_DIAGNOSTICS['roi_too_small'] += 1
        return None
    response, noise = _v32_slot_flatfield_response(roi)
    peak = float(np.percentile(response, 99.0))
    if peak < max(2.0, 3.0 * noise):
        _V32_SLOT_RECOVERY_DIAGNOSTICS['weak_local_contrast'] += 1
        return None
    xt = max(6., V17_SLOT_RECOVER_CENTER_X_FRAC * x_pitch)
    yt = max(6., V17_SLOT_RECOVER_CENTER_Y_FRAC * row_pitch)
    best = None
    for fraction in (0.50, 0.35, 0.65):
        mask = (response >= max(1.0, 2.0 * noise, fraction * peak)).astype(np.uint8) * 255
        mask = _postprocess_binary_v12(mask, 1)
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if not contours:
            _V32_SLOT_RECOVERY_DIAGNOSTICS['no_contour'] += 1
        for contour in contours:
            if contour_pixel_count(contour) < SLOT_MIN_AREA_PX:
                continue
            bx, by, bw, bh = cv2.boundingRect(contour)
            if bx <= 1 or by <= 1 or bx + bw >= roi.shape[1] - 1 or by + bh >= roi.shape[0] - 1:
                _V32_SLOT_RECOVERY_DIAGNOSTICS['touches_roi_border'] += 1
                continue
            contour = contour.copy(); contour[:, 0, 0] += x0; contour[:, 0, 1] += y0
            d = describe_contour(contour)
            if d is None:
                continue
            if abs(float(d.center[0]) - cx) > xt or abs(float(d.center[1]) - cy) > yt:
                _V32_SLOT_RECOVERY_DIAGNOSTICS['outside_predicted_site'] += 1
                continue
            if not _slot_shape_ok_v17(d, area_limits, seed=False):
                _V32_SLOT_RECOVERY_DIAGNOSTICS['area_or_shape_qc'] += 1
                continue
            score = _slot_candidate_site_score_v17(d, cx, cy, xt, yt) + 0.08 * abs(fraction - 0.5)
            if best is None or score < best[0]:
                best = (score, d)
    if best is None:
        return None
    _V32_SLOT_RECOVERY_DIAGNOSTICS['flatfield_recovered'] += 1
    return _mark_source(best[1], 'slot_v32_flatfield_recovered')


def _recover_local_slot_v18_deep(
    gray, expected_xy, x_pitch, row_pitch, expected_w, expected_h, area_limits, aggressive_level,
) -> Optional[ShapeDesc]:
    d = _v32_recover_flatfield_slot(gray, expected_xy, x_pitch, row_pitch, expected_w, expected_h, area_limits)
    if d is not None:
        return d
    d = _recover_local_slot_v18_legacy(gray, expected_xy, x_pitch, row_pitch, expected_w, expected_h, area_limits, aggressive_level)
    _V32_SLOT_RECOVERY_DIAGNOSTICS['legacy_recovered' if d is not None else 'site_not_recovered'] += 1
    return d


def _recover_local_slot_v18_legacy(
    gray: np.ndarray,
    expected_xy: Tuple[float,float],
    x_pitch: float,
    row_pitch: float,
    expected_w: float,
    expected_h: float,
    area_limits: Tuple[float,float],
    aggressive_level: int,
) -> Optional[ShapeDesc]:
    """
    Deep recovery at one predicted 4-row site.

    Contrast/blur handling becomes increasingly permissive. Area, completeness and
    shape thresholds are exactly the V17 recovery thresholds at every level.
    """
    H,W=gray.shape
    # First relocate the center toward the broad low-gray basin. This is position-only.
    basin_xy,_=_slot_darkness_center_v18(gray,expected_xy,x_pitch,row_pitch,expected_w,expected_h)
    # Blend prediction and basin more strongly at deeper levels, but never leave the site cell.
    blend=(0.30,0.50,0.70,0.85)[min(aggressive_level,3)]
    cx=(1-blend)*float(expected_xy[0])+blend*float(basin_xy[0])
    cy=(1-blend)*float(expected_xy[1])+blend*float(basin_xy[1])
    x_tol=max(6.0,V17_SLOT_RECOVER_CENTER_X_FRAC*x_pitch)
    y_tol=max(6.0,V17_SLOT_RECOVER_CENTER_Y_FRAC*row_pitch)
    cx=float(np.clip(cx,float(expected_xy[0])-x_tol,float(expected_xy[0])+x_tol))
    cy=float(np.clip(cy,float(expected_xy[1])-y_tol,float(expected_xy[1])+y_tol))

    half_w=max(14.0,min(0.49*x_pitch,max(2.2*expected_w,18.0)))
    half_h=max(16.0,min(0.74*row_pitch,max(1.15*expected_h+8.0,22.0)))
    x0=max(0,int(math.floor(cx-half_w))); x1=min(W,int(math.ceil(cx+half_w+1)))
    y0=max(0,int(math.floor(cy-half_h))); y1=min(H,int(math.ceil(cy+half_h+1)))
    if x1-x0<9 or y1-y0<9: return None
    roi=gray[y0:y1,x0:x1]
    masks=_v18_slot_masks_at_roi(roi,aggressive_level)

    # Add thresholds interpolated between local foreground basin and surrounding background.
    if len(masks)<V18_SLOT_MAX_LOCAL_MASKS:
        locx=int(round(cx-x0)); locy=int(round(cy-y0))
        rw=max(2,int(round(expected_w/2))); rh=max(2,int(round(expected_h/2)))
        xa=max(0,locx-rw); xb=min(roi.shape[1],locx+rw+1)
        ya=max(0,locy-rh); yb=min(roi.shape[0],locy+rh+1)
        fg=float(np.median(roi[ya:yb,xa:xb])) if xb>xa and yb>ya else float(np.median(roi))
        bg=float(np.median(roi))
        for a in V18_SLOT_INTERPOLATION_ALPHAS:
            th=fg+a*(bg-fg)
            if FEATURE_POLARITY.lower()=='dark': mm=(roi.astype(np.float32)<=th).astype(np.uint8)*255
            else: mm=(roi.astype(np.float32)>=th).astype(np.uint8)*255
            mm=cv2.morphologyEx(mm,cv2.MORPH_CLOSE,np.ones((3,3),np.uint8),iterations=1+(aggressive_level>=3))
            masks[f'interp_{a:g}']=mm
            if len(masks)>=V18_SLOT_MAX_LOCAL_MASKS: break

    candidates=[]; m=V17_SLOT_RECOVER_ROI_BORDER_MARGIN_PX
    for mm in masks.values():
        contours,_=cv2.findContours(mm,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
        for c in contours:
            if contour_pixel_count(c)<SLOT_MIN_AREA_PX: continue
            bx,by,bw,bh=cv2.boundingRect(c)
            # Shape/edge completeness is a hard constraint and is NEVER loosened.
            if bx<=m or by<=m or bx+bw>=roi.shape[1]-m or by+bh>=roi.shape[0]-m: continue
            cg=c.astype(np.int32).copy(); cg[:,0,0]+=x0; cg[:,0,1]+=y0
            d=describe_contour(cg)
            if d is None: continue
            # Center may follow the low-gray basin, but still must belong to this theoretical site.
            if abs(float(d.center[0])-float(expected_xy[0]))>1.10*x_tol: continue
            if abs(float(d.center[1])-float(expected_xy[1]))>1.10*y_tol: continue
            # CRITICAL: same area + same shape requirements at every recovery level.
            if not _slot_shape_ok_v17(d,area_limits,seed=False): continue
            _mark_source(d,'slot_4row_deep_recovered')
            score=_slot_candidate_site_score_v17(d,float(expected_xy[0]),float(expected_xy[1]),1.10*x_tol,1.10*y_tol)
            # Prefer contours whose center also agrees with the broad low-gray basin.
            score += 0.12*float(np.linalg.norm(d.center-np.asarray(basin_xy,dtype=float)))/max(x_pitch,row_pitch,1.0)
            candidates.append((score,d))
    if not candidates: return None
    candidates.sort(key=lambda z:z[0])
    return candidates[0][1]


def _slot_four_row_array_v18(
    gray: np.ndarray,
    pool: List[ShapeDesc],
    strict_seed: List[ShapeDesc],
    image_name: str,
    area_limits: Tuple[float,float],
) -> Tuple[List[ShapeDesc],List[int],Dict[str,object]]:
    """Hard 4-row prior + repeated blur/contrast-only recovery until no progress."""
    H,W=gray.shape
    # Slot's own 4-row model remains horizontal-row based; '-1' is a VIA-specific diamond
    # requirement in V18. We do not force diagonal slot logic from the filename token.
    allow_diag=False
    relaxed=[d for d in pool if _slot_shape_ok_v17(d,area_limits,seed=True)]
    model_cands=_dedup_desc_list_v16(strict_seed+relaxed)
    row_model=_fit_four_row_centers_v17(model_cands,H,image_name)

    lowgray_loc=[]
    if row_model is None or len(model_cands)<8:
        lowgray_loc=_slot_lowgray_model_proposals_v18(gray)
        # LOW-GRAY proposals are used ONLY to infer where the 4 rows/pitch are.
        row_model=_fit_four_row_centers_v17(_dedup_desc_list_v16(model_cands+lowgray_loc),H,image_name)
    meta={
        'array_assist_used':False,'array_lattice_type':'slot_4row_H_deep_lowgray',
        'array_x_pitch_px':np.nan,'array_y_pitch_px':np.nan,
        'center_horizontal_pitch_px':np.nan,'center_vertical_pitch_px':np.nan,
        'center_diagonal_component_px':np.nan,'fixed_diagonal_enabled':False,'fixed_pitch_frozen':True,
        'center_recursive_used':False,'center_recursive_attempts':0,'center_recursive_recovered':0,
        'center_recursive_boundary_hits':0,'array_recovered_count':0,'array_final_recovered_count':0,
        'array_expected_positions':0,'array_final_predicted_positions':0,
        'array_recovery_attempts':0,'array_final_recovery_attempts':0,
        'slot_expected_rows':4,'slot_row_pitch_px':np.nan,'slot_horizontal_pitch_px':np.nan,
        'slot_phase_shift_per_row_px':0.0,'slot_seed_rows_found':0,
        'slot_row1_count':0,'slot_row2_count':0,'slot_row3_count':0,'slot_row4_count':0,
        'slot_row1_expected_sites':0,'slot_row2_expected_sites':0,'slot_row3_expected_sites':0,'slot_row4_expected_sites':0,
        'slot_lowgray_locator_proposals':int(len(lowgray_loc)),'slot_deep_recovery_passes':0,
    }
    if row_model is None:
        return [],[],meta
    row_centers=np.asarray(row_model['row_centers'],dtype=float); row_pitch=float(row_model['row_pitch']); row_tol=float(row_model['row_tol'])
    fit_cands=_dedup_desc_list_v16(model_cands+lowgray_loc)
    row_map=_assign_to_row_centers_v17(fit_cands,row_centers,row_tol)
    # Prefer shape-valid candidates for H pitch; only use low-gray locators if necessary.
    shape_row_map=_assign_to_row_centers_v17(model_cands,row_centers,row_tol)
    hp=_slot_horizontal_pitch_v17(shape_row_map,W)
    if not np.isfinite(hp) or hp<V16_PITCH_MIN_PX:
        hp=_slot_horizontal_pitch_v17(row_map,W)
    if not np.isfinite(hp) or hp<V16_PITCH_MIN_PX:
        hp=_fixed_axis_pitch_from_centers_v16(fit_cands,'x',gray.shape)
    if not np.isfinite(hp) or hp<V16_PITCH_MIN_PX:
        return [],[],meta

    # Slots are four horizontal rows; use one common horizontal phase.
    phase,_=_fit_slot_phase_model_v17(shape_row_map if any(shape_row_map[r] for r in range(4)) else row_map,hp,False)
    if not np.isfinite(phase):
        phase,_=_fit_slot_phase_model_v17(row_map,hp,False)
    if not np.isfinite(phase): return [],[],meta

    widths=[max(1.0,_contour_xy_bounds(d)[1]-_contour_xy_bounds(d)[0]) for d in model_cands]
    heights=[max(1.0,_contour_xy_bounds(d)[3]-_contour_xy_bounds(d)[2]) for d in model_cands]
    if not widths:
        # Low-gray boxes are localization only, but they provide a conservative ROI size fallback.
        widths=[max(5.0,0.18*hp)]
    if not heights: heights=[max(10.0,0.45*row_pitch)]
    med_w=float(np.median(widths)); global_h=float(np.percentile(heights,85))
    margin=max(4.0,0.55*med_w)
    expected_by_row={}
    for r in range(4):
        xs=_expected_x_sites_v17(W,phase,0.0,r,hp,margin)
        expected_by_row[r]=xs; meta[f'slot_row{r+1}_expected_sites']=len(xs)
    total_sites=sum(len(v) for v in expected_by_row.values())

    selected=[]; row_ids=[]; used_sites=set()
    x_tol=max(6.0,V17_SLOT_SITE_X_TOL_FRAC*hp); y_tol=max(6.0,V17_SLOT_SITE_Y_TOL_FRAC*row_pitch)
    # First consume all already segmented SHAPE-VALID candidates.
    for r in range(4):
        ey=float(row_centers[r])
        for si,ex in enumerate(expected_by_row[r]):
            d=_pick_pool_slot_at_site_v17(pool,selected,ex,ey,x_tol,y_tol,area_limits)
            if d is None: continue
            if _source_of(d)=='primary': _mark_source(d,'slot_4row_selected')
            selected.append(d); row_ids.append(r+1); used_sites.add((r,si))

    attempts=recovered=0; no_progress=0; pass_count=0
    for pass_no in range(1,V18_SLOT_DEEP_RECOVERY_PASSES+1):
        pass_count=pass_no; before=len(selected)
        counts_now=[sum(1 for rr in row_ids if rr==r+1) for r in range(4)]
        for r in range(4):
            # Every row must be actively searched. Later passes focus naturally on missing sites.
            row_heights=[max(1.0,_contour_xy_bounds(d)[3]-_contour_xy_bounds(d)[2]) for d,rr in zip(selected,row_ids) if rr==r+1]
            expected_h=float(np.median(row_heights)) if row_heights else global_h
            ey=float(row_centers[r])
            for si,ex in enumerate(expected_by_row[r]):
                if (r,si) in used_sites: continue
                attempts+=1
                if attempts==1 or attempts%V17_SLOT_PROGRESS_EVERY==0:
                    print(f"      [{image_name}] V18短槽4行深度补点 {attempts}; pass={pass_no}; row={r+1}; recovered={recovered} ...",flush=True)
                # Existing pool first: center tolerance can be broad, but shape/area cannot.
                d=_pick_pool_slot_at_site_v17(pool,selected,ex,ey,(1.10+0.08*pass_no)*x_tol,(1.10+0.06*pass_no)*y_tol,area_limits)
                if d is None:
                    d=_recover_local_slot_v18_deep(gray,(ex,ey),hp,row_pitch,med_w,expected_h,area_limits,pass_no-1)
                    if d is not None: pool.append(d)
                if d is None or _desc_close_to_any_selected_v12(d,selected): continue
                selected.append(d); row_ids.append(r+1); used_sites.add((r,si)); recovered+=1
        gained=len(selected)-before
        if gained==0: no_progress+=1
        else: no_progress=0
        counts_after=[sum(1 for rr in row_ids if rr==r+1) for r in range(4)]
        # Stop early only when every row is already close to the known ~10-per-row layout.
        if min(counts_after)>=max(7,V18_SLOT_TARGET_APPROX_PER_ROW-2):
            break
        if V18_SLOT_ROW_SEARCH_UNTIL_NO_PROGRESS and no_progress>=V18_SLOT_MAX_NO_PROGRESS_PASSES:
            break

    order=sorted(range(len(selected)),key=lambda k:(row_ids[k],selected[k].center[0],selected[k].center[1]))
    selected=[selected[k] for k in order]; row_ids=[row_ids[k] for k in order]
    counts=[sum(1 for rr in row_ids if rr==r+1) for r in range(4)]
    meta.update({
        'array_assist_used':True,'array_x_pitch_px':float(hp),'array_y_pitch_px':float(row_pitch),
        'center_horizontal_pitch_px':float(hp),'center_vertical_pitch_px':float(row_pitch),
        'center_recursive_used':True,'center_recursive_attempts':int(attempts),'center_recursive_recovered':int(recovered),
        'array_recovered_count':int(recovered),'array_final_recovered_count':int(recovered),
        'array_expected_positions':int(total_sites),'array_final_predicted_positions':int(total_sites),
        'array_recovery_attempts':int(attempts),'array_final_recovery_attempts':int(attempts),
        'slot_row_pitch_px':float(row_pitch),'slot_horizontal_pitch_px':float(hp),'slot_phase_shift_per_row_px':0.0,
        'slot_seed_rows_found':int(row_model['occupied_seed_rows']),
        'slot_row1_count':int(counts[0]),'slot_row2_count':int(counts[1]),'slot_row3_count':int(counts[2]),'slot_row4_count':int(counts[3]),
        'array_occupancy_seed':float(len(model_cands)/max(1,total_sites)),
        'slot_deep_recovery_passes':int(pass_count),
    })
    return selected,row_ids,meta


def measure_slot_v18(
    gray: np.ndarray, pool: List[ShapeDesc], px_nm: float, image_name: str, condition: str
) -> Tuple[List[dict],List[ShapeDesc],Dict[int,ShapeDesc],List[dict],Dict[str,object]]:
    """V18 short-slot: hard 4-row layout + blur/contrast-only deep recovery."""
    strict,rejects,reg_by_id,mean_area,med_area,_=_slot_initial_filter_v12(pool)
    strict=_dedup_desc_list_v16(strict)
    if np.isfinite(mean_area) and mean_area>0:
        area_lo=max(SLOT_MIN_AREA_PX,mean_area*SLOT_AREA_MIN_MEAN_FACTOR); area_hi=mean_area*SLOT_AREA_MAX_MEAN_FACTOR
    else:
        aa=np.asarray([contour_pixel_count(d.contour) for d in strict],dtype=float)
        # If strict segmentation is sparse, infer reference area only from shape-valid pool candidates.
        if len(aa)==0:
            shp=[]
            for d in pool:
                reg=slot_regularity_metrics(d); pa=contour_pixel_count(d.contour)
                if pa>=SLOT_MIN_AREA_PX and reg['regularity_solidity']>=V17_SLOT_SEED_MIN_SOLIDITY and reg['regularity_convexity']>=V17_SLOT_SEED_MIN_CONVEXITY and reg['regularity_template_iou']>=V17_SLOT_SEED_MIN_TEMPLATE_IOU:
                    shp.append(pa)
            aa=np.asarray(shp,dtype=float)
        mean_area=float(np.mean(aa)) if len(aa) else np.nan; med_area=float(np.median(aa)) if len(aa) else np.nan
        area_lo=max(SLOT_MIN_AREA_PX,mean_area*SLOT_AREA_MIN_MEAN_FACTOR) if np.isfinite(mean_area) else SLOT_MIN_AREA_PX
        area_hi=mean_area*SLOT_AREA_MAX_MEAN_FACTOR if np.isfinite(mean_area) else np.inf
    area_limits=(float(area_lo),float(area_hi))
    # Add only geometry-valid seeds; low-gray-only anchors are handled inside V18 model fitting.
    if len(strict)<4:
        extra=[d for d in pool if _slot_shape_ok_v17(d,area_limits,seed=True)]
        strict=_dedup_desc_list_v16(strict+extra)
    selected,row_ids,meta=_slot_four_row_array_v18(gray,pool,strict,image_name,area_limits)

    rows=[]; reps={}
    for oid,(d,rid) in enumerate(zip(selected,row_ids),1):
        pts=d.contour.reshape(-1,2).astype(np.float64)
        width_px=float(pts[:,0].max()-pts[:,0].min()); length_px=float(pts[:,1].max()-pts[:,1].min())
        reg=reg_by_id.get(id(d),slot_regularity_metrics(d)); pa=contour_pixel_count(d.contour)
        rows.append({
            'condition':condition,'image':image_name,'pattern':'slot','object_id':oid,'slot_row':int(rid),
            'slot_width_nm':width_px*px_nm,'slot_length_nm':length_px*px_nm,
            'slot_aspect_ratio':length_px/max(width_px,1e-12),'slot_angle_from_y_deg':abs(90.0-d.angle_from_x_deg),
            'slot_pixel_area_px':pa,'slot_mean_candidate_area_px':mean_area,'slot_reference_median_area_px':med_area,
            'slot_area_ratio_to_mean':pa/mean_area if np.isfinite(mean_area) and mean_area>0 else np.nan,
            'slot_regularity_solidity':reg['regularity_solidity'],'slot_regularity_convexity':reg['regularity_convexity'],
            'slot_rect_iou':reg['regularity_rect_iou'],'slot_ellipse_iou':reg['regularity_ellipse_iou'],
            'slot_best_template_iou':reg['regularity_template_iou'],'selection_source':_source_of(d),'pixel_size_nm':px_nm,
        })
        if rid not in reps: reps[rid]=d

    rejects=[r for r in rejects if r.get('_desc') is None or not _desc_close_to_any_selected_v12(r['_desc'],selected)]
    rejected_ids={id(r.get('_desc')) for r in rejects if r.get('_desc') is not None}
    for d in pool:
        if _desc_close_to_any_selected_v12(d,selected) or id(d) in rejected_ids: continue
        pa=contour_pixel_count(d.contour); reg=slot_regularity_metrics(d)
        if _source_of(d)!='primary' and not _slot_shape_ok_v17(d,area_limits,seed=False): continue
        codes=[]; reasons=[]
        if pa<SLOT_MIN_AREA_PX:
            codes.append('slot_area_too_small'); reasons.append(f'filled pixel area {pa} px < {SLOT_MIN_AREA_PX} px')
        if np.isfinite(mean_area) and mean_area>0:
            if pa<area_lo: codes.append('slot_area_below_0.1x_mean'); reasons.append(f'area {pa}px < {area_lo:.1f}px')
            if pa>area_hi: codes.append('slot_area_above_10x_mean'); reasons.append(f'area {pa}px > {area_hi:.1f}px')
        if reg['regularity_solidity']<V17_SLOT_RECOVER_MIN_SOLIDITY:
            codes.append('slot_low_regularity_solidity'); reasons.append(f"solidity {reg['regularity_solidity']:.3f} too low")
        if reg['regularity_convexity']<V17_SLOT_RECOVER_MIN_CONVEXITY:
            codes.append('slot_low_convexity'); reasons.append(f"convexity {reg['regularity_convexity']:.3f} too low")
        if reg['regularity_template_iou']<V17_SLOT_RECOVER_MIN_TEMPLATE_IOU:
            codes.append('slot_irregular_template_match'); reasons.append(f"template IoU {reg['regularity_template_iou']:.3f} too low")
        if not codes:
            codes=['slot_center_not_used_by_4row_model']; reasons=['shape/area-valid short-slot candidate was not the best center match for one of the four-row theoretical sites']
        rejects.append(_make_reject_from_contour(
            d.contour,'slot_v18_4row_deep_filter','+'.join(codes),'; '.join(reasons),desc=d,
            source=_source_of(d),pixel_area_px=pa,**reg,
        ))
    return rows,selected,reps,rejects,meta



# ============================================================
# V19: TWO-ROOT / EIGHT-CONDITION BEFORE-AFTER PIPELINE
# ============================================================

# Main paths are defined at the top of the file so they are easy to edit.

CONDITION_IDS_V19 = [str(i) for i in range(1, 9)]
CONDITION_LABELS_V19 = {str(i): f"Condition {i}" for i in range(1, 9)}
STAGE_ORDER_V19 = ["before", "after"]
STAGE_LABELS_ZH_V19 = {"before": "处理前", "after": "处理后"}
PATTERN_KEYS_V19 = ["via40", "via60", "trench160", "slot210"]
PATTERN_LABELS_ZH_V19 = {
    "via40": "Via40（标称图案）",
    "via60": "Via60（标称图案）",
    "trench160": "160开头相对trench",
    "slot210": "210开头slot最下排",
}
PREFIX_TO_PATTERN_V19 = {
    "40": "via40",
    "60": "via60",
    "160": "trench160",
    "210": "slot210",
}

# Image discovery. JPG/JPEG are the expected formats; TIFF/PNG support is retained.
IMAGE_SUFFIXES_V19 = {".jpg", ".jpeg", ".tif", ".tiff", ".png"}
OUTPUT_DIR_NAME_V19 = "sem_before_after_outputV25"
SAVE_BINARY_MASK_V19 = False
DRAW_REJECTED_OBJECTS_V19 = False

# Pixel/geometry output.
MEASUREMENT_EXTENT_TRIM_PERCENT_V19 = 0.0  # 0 = contour min/max; e.g. 1.0 trims 1% at both ends
PLOT_ERROR_BAR_V19 = "std"               # "std", "sem", or "none"
PLOT_DPI_V19 = 220

# More conservative post-QC for final objects. These work in addition to V18 QC.
V19_DYNAMIC_BORDER_MIN_PX = 4.0
V19_DYNAMIC_BORDER_SIZE_FACTOR = 0.12
V19_AREA_LO_FACTOR = 0.35
V19_AREA_HI_FACTOR = 2.85
V19_ARRAY_LINK_MAX_PITCH_FACTOR = 1.70
V19_ARRAY_LINK_MIN_PITCH_FACTOR = 0.30

# Via bottom-row / axis fitting.
V19_VIA_BOTTOM_ROW_TOL_DIAM_FACTOR = 0.58
V19_VIA_BOTTOM_ROW_MAX_BELOW_TOL_FACTOR = 1.8
V19_VIA_MIN_ROW_X_COMPONENT = 0.35
V19_VIA_MIN_BOTTOM_ROW_OBJECTS = 2

# Trench-specific enhanced recognition. This replaces the old one-mask-only pairing.
V19_TRENCH_MIN_VERTICAL_ASPECT = 1.25
V19_TRENCH_MAX_ANGLE_FROM_Y_DEG = 35.0
V19_TRENCH_MIN_LENGTH_PX = 10.0
V19_TRENCH_MIN_SOLIDITY = 0.42
V19_TRENCH_MIN_WIDTH_PX = 2.0
V19_TRENCH_MAX_X_OFFSET_FACTOR = 2.25
V19_TRENCH_MIN_X_OVERLAP_FRAC = 0.12
V19_TRENCH_MAX_WIDTH_MISMATCH_FRAC = 1.15
V19_TRENCH_MAX_GAP_WIDTH_FACTOR = 40.0
V19_TRENCH_GAP_OUTLIER_FRAC = 0.40
V19_TRENCH_VERTICAL_CLOSE_KERNELS = (3, 5)
V19_TRENCH_MAX_POOL = 240

# The generic V14 lattice is used for 40*/60* vias because their array may be rotated.
# Disable the older post-lattice H/V/45-degree recursive expansion; V19 does its own
# final strict QC and does not need that aggressive extra pass.
CENTER_RECURSIVE_SEARCH_VIA = False


def _read_gray_image_unicode_v19(path: Path) -> np.ndarray:
    """Unicode-safe OpenCV reader for Windows paths."""
    try:
        raw = np.fromfile(str(path), dtype=np.uint8)
        img = cv2.imdecode(raw, cv2.IMREAD_UNCHANGED)
    except Exception:
        img = None
    if img is None:
        return read_gray_image(path)
    if img.ndim == 3:
        if img.shape[2] == 4:
            img = cv2.cvtColor(img, cv2.COLOR_BGRA2GRAY)
        else:
            img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    img = img.astype(np.float32)
    lo, hi = np.percentile(img, [1.0, 99.0])
    if hi <= lo:
        lo, hi = float(img.min()), float(img.max())
    if hi <= lo:
        return np.zeros_like(img, dtype=np.uint8)
    return np.clip((img - lo) * 255.0 / (hi - lo), 0, 255).astype(np.uint8)


def _imwrite_unicode_v19(path: Path, image: np.ndarray) -> None:
    """Unicode-safe image writer for Windows paths."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ext = path.suffix.lower() or ".png"
    ok, buf = cv2.imencode(ext, image)
    if not ok:
        raise RuntimeError(f"OpenCV could not encode output image: {path}")
    buf.tofile(str(path))


def _find_sem_images_v19(condition_dir: Path) -> List[Path]:
    if not condition_dir.exists():
        return []
    return sorted(
        [p for p in condition_dir.iterdir() if p.is_file() and p.suffix.lower() in IMAGE_SUFFIXES_V19],
        key=lambda p: p.name.lower(),
    )


def _classify_pattern_key_v19(image_path: Path) -> Optional[str]:
    stem = image_path.stem.strip()
    # Longest first is explicit and avoids any future prefix ambiguity.
    for prefix in ("210", "160", "60", "40"):
        if stem.startswith(prefix):
            return PREFIX_TO_PATTERN_V19[prefix]
    return None


def _base_pattern_v19(pattern_key: str) -> str:
    if pattern_key in {"via40", "via60"}:
        return "via"
    if pattern_key == "trench160":
        return "trench"
    if pattern_key == "slot210":
        return "slot"
    return "unknown"


def _safe_stem_v19(path: Path) -> str:
    return re.sub(r"[^0-9A-Za-z._-]+", "_", path.stem).strip("_") or "image"


def _bbox_complete_v19(
    d: ShapeDesc,
    image_shape: Tuple[int, int],
    size_factor: float = V19_DYNAMIC_BORDER_SIZE_FACTOR,
) -> bool:
    """Reject clipped objects using both a fixed and object-size-scaled border margin."""
    h, w = image_shape[:2]
    x, y, bw, bh = cv2.boundingRect(d.contour.astype(np.int32))
    margin = max(V19_DYNAMIC_BORDER_MIN_PX, size_factor * max(1.0, min(bw, bh)))
    return bool(
        x > margin
        and y > margin
        and x + bw < w - margin
        and y + bh < h - margin
    )


def _projection_span_v19(contour: np.ndarray, axis: np.ndarray) -> float:
    pts = contour.reshape(-1, 2).astype(np.float64)
    axis = np.asarray(axis, dtype=np.float64)
    axis /= np.linalg.norm(axis) + 1e-12
    q = pts @ axis
    trim = float(np.clip(MEASUREMENT_EXTENT_TRIM_PERCENT_V19, 0.0, 20.0))
    if trim > 0 and len(q) >= 20:
        lo, hi = np.percentile(q, [trim, 100.0 - trim])
        return float(max(0.0, hi - lo))
    return float(max(0.0, q.max() - q.min()))


def _dedup_v19(items: List[ShapeDesc]) -> List[ShapeDesc]:
    return _dedup_desc_list_v16(items)


def _area_consistency_filter_v19(
    items: List[ShapeDesc],
    stage_name: str,
    lo_factor: float = V19_AREA_LO_FACTOR,
    hi_factor: float = V19_AREA_HI_FACTOR,
) -> Tuple[List[ShapeDesc], List[dict]]:
    if len(items) < 3:
        return list(items), []
    areas = np.asarray([contour_pixel_count(d.contour) for d in items], dtype=float)
    med = float(np.median(areas))
    if not np.isfinite(med) or med <= 0:
        return list(items), []
    kept: List[ShapeDesc] = []
    rejected: List[dict] = []
    for d, a in zip(items, areas):
        if lo_factor * med <= a <= hi_factor * med:
            kept.append(d)
        else:
            rejected.append(_make_reject_from_contour(
                d.contour,
                stage_name,
                "area_inconsistent_with_image_population",
                f"filled area {a:.1f}px is outside [{lo_factor:.2f}, {hi_factor:.2f}] x median {med:.1f}px",
                desc=d,
                pixel_area_px=float(a),
                median_area_px=med,
            ))
    return (kept if len(kept) >= 2 else list(items)), rejected


def _nearest_neighbor_pitch_v19(items: List[ShapeDesc]) -> float:
    if len(items) < 2:
        return np.nan
    p = np.asarray([d.center for d in items], dtype=float)
    nearest = []
    for i in range(len(p)):
        dist = np.linalg.norm(p - p[i], axis=1)
        dist[i] = np.inf
        v = float(np.min(dist))
        if np.isfinite(v) and v > 0:
            nearest.append(v)
    return float(np.median(nearest)) if nearest else np.nan


def _largest_center_component_v19(
    items: List[ShapeDesc],
    reject_stage: str,
) -> Tuple[List[ShapeDesc], List[dict], float]:
    """Keep the largest nearest-neighbor-connected array component."""
    if len(items) < 4:
        return list(items), [], _nearest_neighbor_pitch_v19(items)
    pitch = _nearest_neighbor_pitch_v19(items)
    if not np.isfinite(pitch) or pitch <= 0:
        return list(items), [], pitch
    p = np.asarray([d.center for d in items], dtype=float)
    lo = V19_ARRAY_LINK_MIN_PITCH_FACTOR * pitch
    hi = V19_ARRAY_LINK_MAX_PITCH_FACTOR * pitch
    adj = [[] for _ in items]
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            dd = float(np.linalg.norm(p[i] - p[j]))
            if lo <= dd <= hi:
                adj[i].append(j)
                adj[j].append(i)
    comps: List[List[int]] = []
    seen = set()
    for i in range(len(items)):
        if i in seen:
            continue
        stack = [i]
        seen.add(i)
        comp = []
        while stack:
            k = stack.pop()
            comp.append(k)
            for j in adj[k]:
                if j not in seen:
                    seen.add(j)
                    stack.append(j)
        comps.append(comp)
    comps.sort(key=lambda c: (len(c), sum(len(adj[i]) for i in c)), reverse=True)
    best = set(comps[0]) if comps else set(range(len(items)))
    # Do not collapse a small but valid data set to a singleton.
    if len(best) < max(3, int(math.ceil(0.45 * len(items)))):
        return list(items), [], pitch
    kept = [d for i, d in enumerate(items) if i in best]
    rejected = [
        _make_reject_from_contour(
            d.contour,
            reject_stage,
            "isolated_from_main_array",
            "regular candidate is isolated from the largest equal-spacing center component",
            desc=d,
            estimated_pitch_px=pitch,
        )
        for i, d in enumerate(items) if i not in best
    ]
    return kept, rejected, pitch


def _fit_tls_line_v19(points: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    pts = np.asarray(points, dtype=float)
    origin = np.mean(pts, axis=0)
    if len(pts) < 2:
        return origin, np.array([1.0, 0.0], dtype=float)
    _, _, vt = np.linalg.svd(pts - origin, full_matrices=False)
    u = vt[0].astype(float)
    u /= np.linalg.norm(u) + 1e-12
    if u[0] < 0:
        u = -u
    return origin, u


def _fit_via_bottom_row_axes_v19(
    items: List[ShapeDesc], image_shape: Tuple[int, int]
) -> Dict[str, object]:
    """Fit the image-lowest collinear via row, then create X and upward Y axes."""
    centers = np.asarray([d.center for d in items], dtype=float)
    if len(items) == 0:
        return {
            "origin": np.array([0.0, 0.0]),
            "x_axis": np.array([1.0, 0.0]),
            "y_axis_up": np.array([0.0, -1.0]),
            "bottom_indices": [],
            "row_angle_deg": 0.0,
            "row_fit_rms_px": np.nan,
            "fallback": True,
        }
    diam = float(np.median([
        max(1.0, 0.5 * ((_contour_xy_bounds(d)[1] - _contour_xy_bounds(d)[0]) +
                        (_contour_xy_bounds(d)[3] - _contour_xy_bounds(d)[2])))
        for d in items
    ]))
    tol = max(2.5, V19_VIA_BOTTOM_ROW_TOL_DIAM_FACTOR * diam)
    min_sep = max(4.0, 1.10 * diam)
    best = None
    n = len(items)
    for i in range(n):
        for j in range(i + 1, n):
            v = centers[j] - centers[i]
            dist = float(np.linalg.norm(v))
            if dist < min_sep:
                continue
            u = v / dist
            if u[0] < 0:
                u = -u
            # A bottom row must have a meaningful left-right image component.
            if abs(float(u[0])) < V19_VIA_MIN_ROW_X_COMPONENT:
                continue
            n_down = np.array([-u[1], u[0]], dtype=float)
            if n_down[1] < 0:
                n_down = -n_down
            t = centers @ n_down
            t0 = 0.5 * (t[i] + t[j])
            inliers = np.where(np.abs(t - t0) <= tol)[0]
            if len(inliers) < V19_VIA_MIN_BOTTOM_ROW_OBJECTS:
                continue
            t_in = t[inliers]
            below = float(np.max(t) - np.median(t_in))
            if below > V19_VIA_BOTTOM_ROW_MAX_BELOW_TOL_FACTOR * tol:
                continue
            residual = float(np.sqrt(np.mean((t_in - np.median(t_in)) ** 2)))
            along = np.sort(centers[inliers] @ u)
            span = float(along[-1] - along[0]) if len(along) >= 2 else 0.0
            # Lexicographic objective: row population, bottom-envelope quality,
            # left-right character, span, and residual.
            rec = (
                len(inliers),
                -below / max(tol, 1e-9),
                abs(float(u[0])),
                span / max(diam, 1e-9),
                -residual / max(tol, 1e-9),
                inliers,
            )
            if best is None or rec[:5] > best[:5]:
                best = rec

    fallback = False
    if best is None:
        fallback = True
        k = min(len(items), max(2, int(math.ceil(0.30 * len(items)))))
        idx = np.argsort(centers[:, 1])[-k:]
    else:
        idx = np.asarray(best[5], dtype=int)

    origin, x_axis = _fit_tls_line_v19(centers[idx])
    # Recompute the row after TLS refinement.
    n_down = np.array([-x_axis[1], x_axis[0]], dtype=float)
    if n_down[1] < 0:
        n_down = -n_down
    t = centers @ n_down
    t0 = float(origin @ n_down)
    idx2 = np.where(np.abs(t - t0) <= tol)[0]
    if len(idx2) >= V19_VIA_MIN_BOTTOM_ROW_OBJECTS:
        # Keep only a line near the true lower envelope.
        if float(np.max(t) - np.median(t[idx2])) <= V19_VIA_BOTTOM_ROW_MAX_BELOW_TOL_FACTOR * tol:
            idx = idx2
            origin, x_axis = _fit_tls_line_v19(centers[idx])

    y_up = np.array([-x_axis[1], x_axis[0]], dtype=float)
    if y_up[1] > 0:  # image Y grows downward; upward must have negative image-Y component
        y_up = -y_up
    if x_axis[0] < 0:
        x_axis = -x_axis
    residuals = (centers[idx] - origin) @ y_up if len(idx) else np.asarray([], dtype=float)
    rms = float(np.sqrt(np.mean(residuals ** 2))) if len(residuals) else np.nan
    angle = float(math.degrees(math.atan2(x_axis[1], x_axis[0])))
    return {
        "origin": origin,
        "x_axis": x_axis,
        "y_axis_up": y_up,
        "bottom_indices": [int(i) for i in idx],
        "row_angle_deg": angle,
        "row_fit_rms_px": rms,
        "row_tolerance_px": tol,
        "fallback": fallback,
    }


# -----------------------------
# V21 VIA: PHYSICAL-SIZE MATCHED DARK-DISK + BRIGHT-RING + RECTANGULAR LATTICE
# -----------------------------
# Why V21 exists:
#   V20 could lock onto the darkest INNER CORE of a via and therefore report objects
#   far smaller than the real hole. V21 uses the known physical scale as a strong prior:
#       40* image -> nominal via diameter ~40 nm
#       60* image -> nominal via diameter ~60 nm
#   PixelSize converts those sizes to pixels before any via localization is attempted.
#
# Recognition flow:
#   1) scale-matched filter: dark disk around one via diameter + brighter outer annulus;
#   2) local-maxima detection at that physical scale (tiny dark speckles are suppressed);
#   3) rectangular/affine center-lattice fit and missing-site recovery;
#   4) true edge measurement near the NOMINAL radius, using dark->bright gradient plus
#      half-gray crossing. Search is forbidden close to the center, so an inner dark core
#      cannot become the measured via diameter;
#   5) local half-level contour fallback if a 1-D edge is too fuzzy.
#
# Trench and slot logic below this section are intentionally unchanged.

V21_VIA_NOMINAL_DIAM_NM = {"via40": 40.0, "via60": 60.0}
# V32 adapts the reference diameter to the image. Keep a local template range:
# oversized annuli can overlap adjacent holes and pull the center off its site.
V21_VIA_SCALE_RATIOS = (0.80, 0.90, 1.00, 1.10, 1.22)
V21_VIA_INNER_RADIUS_DIAM_FRAC = 0.30
V21_VIA_RING_INNER_DIAM_FRAC = 0.56
V21_VIA_RING_OUTER_DIAM_FRAC = 0.82
V21_VIA_SCORE_LOCAL_DARK_WEIGHT = 0.12
V21_VIA_MAX_PEAKS = 150
V21_VIA_PEAK_SCORE_MIN = 0.27
V21_VIA_MIN_RAW_CONTRAST_GRAY = 0.35
V21_VIA_PEAK_NMS_DIAM_FRAC = 0.52
V21_VIA_BORDER_DIAM_FRAC = 0.60
V21_VIA_LEGACY_MIN_DIAM_FACTOR = 0.50
V21_VIA_LEGACY_MAX_DIAM_FACTOR = 2.00
V21_VIA_LATTICE_MIN_PITCH_DIAM_FACTOR = 1.12
V21_VIA_LATTICE_MAX_PITCH_DIAM_FACTOR = 8.0
V21_VIA_LATTICE_MIN_BASIS_ANGLE = 52.0
V21_VIA_LATTICE_MAX_BASIS_ANGLE = 128.0
V21_VIA_LATTICE_MAX_EXPECTED = 180
V21_VIA_SITE_REFINE_DIAM_FRAC = 0.34
V21_VIA_INTERNAL_SITE_SCORE_MIN = 0.28
V21_VIA_EXTERNAL_SITE_SCORE_MIN = 0.36
V21_VIA_EXTERNAL_RAW_FACTOR = 1.05
V21_VIA_EDGE_MIN_RADIUS_FACTOR = 0.62
V21_VIA_EDGE_MAX_RADIUS_FACTOR = 1.48
V21_VIA_EDGE_PRIOR_SIGMA_RADIUS_FACTOR = 0.38
V21_VIA_EDGE_GRADIENT_WEIGHT = 0.72
V21_VIA_EDGE_HALFMID_WEIGHT = 0.28
V21_VIA_EDGE_PERP_BAND_DIAM_FRAC = 0.13
V21_VIA_MEASURE_MIN_DIAM_FACTOR = 0.70
V21_VIA_MEASURE_MAX_DIAM_FACTOR = 1.55
V21_VIA_LOCAL_CONTOUR_ROI_DIAM_FACTOR = 1.05
V21_VIA_LOCAL_CONTOUR_CLOSE_FRAC = 0.07
V21_VIA_MIN_FINAL_OBJECTS = 3


def _robust01_v20(arr: np.ndarray, lo_q: float = 5.0, hi_q: float = 97.0) -> np.ndarray:
    """Shared robust normalization; the V20 trench code below also uses this name."""
    a = np.asarray(arr, dtype=np.float32)
    if a.size == 0:
        return np.zeros_like(a, dtype=np.float32)
    lo, hi = np.percentile(a, [lo_q, hi_q])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo + 1e-6:
        lo, hi = float(np.min(a)), float(np.max(a))
    if hi <= lo + 1e-6:
        return np.zeros_like(a, dtype=np.float32)
    return np.clip((a - lo) / (hi - lo), 0.0, 1.0).astype(np.float32)


def _via_nominal_diam_px_v21(pattern_key: str, px_nm: float) -> Tuple[float, float]:
    nominal_nm = float(V21_VIA_NOMINAL_DIAM_NM.get(pattern_key, 40.0))
    nominal_px = nominal_nm / max(float(px_nm), 1e-9)
    return nominal_nm, max(6.0, nominal_px)


def _disk_annulus_kernel_v21(diam_px: float) -> np.ndarray:
    """Kernel response = mean(outer annulus) - mean(inner disk).

    A real dark via of approximately diam_px gives a strong positive response. A tiny
    dark speck does not fill the inner disk and is therefore strongly suppressed.
    """
    d = max(6.0, float(diam_px))
    r0 = max(1.5, V21_VIA_INNER_RADIUS_DIAM_FRAC * d)
    r1 = max(r0 + 1.0, V21_VIA_RING_INNER_DIAM_FRAC * d)
    r2 = max(r1 + 1.0, V21_VIA_RING_OUTER_DIAM_FRAC * d)
    half = int(math.ceil(r2)) + 1
    yy, xx = np.mgrid[-half:half + 1, -half:half + 1]
    rr = np.sqrt(xx.astype(np.float64) ** 2 + yy.astype(np.float64) ** 2)
    inner = rr <= r0
    ring = (rr >= r1) & (rr <= r2)
    k = np.zeros(rr.shape, dtype=np.float32)
    ni = int(np.count_nonzero(inner)); nr = int(np.count_nonzero(ring))
    if ni > 0:
        k[inner] = -1.0 / ni
    if nr > 0:
        k[ring] = 1.0 / nr
    return k


def _via_scale_response_v21(
    gray: np.ndarray, nominal_diam_px: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict[str, object]]:
    """Return scale-matched score, raw annulus-inner contrast and winning diameter map."""
    g = gray.astype(np.float32)
    min_dim = min(gray.shape[:2])
    # Very mild blur suppresses pixel noise without moving the physical edge much.
    sigma = max(0.65, min(1.6, 0.025 * nominal_diam_px))
    gs = cv2.GaussianBlur(g, (0, 0), sigmaX=sigma, sigmaY=sigma)

    norm_responses: List[np.ndarray] = []
    raw_responses: List[np.ndarray] = []
    scale_diams: List[float] = []
    for ratio in V21_VIA_SCALE_RATIOS:
        d = max(6.0, nominal_diam_px * float(ratio))
        kernel = _disk_annulus_kernel_v21(d)
        raw = cv2.filter2D(gs, cv2.CV_32F, kernel, borderType=cv2.BORDER_REFLECT_101)
        # Negative response is not dark-hole evidence.
        raw = np.maximum(raw, 0.0)
        norm = _robust01_v20(raw, 25.0, 99.65)
        raw_responses.append(raw)
        norm_responses.append(norm)
        scale_diams.append(d)

    stack = np.stack(norm_responses, axis=0)
    winner = np.argmax(stack, axis=0)
    best = np.max(stack, axis=0)
    raw_stack = np.stack(raw_responses, axis=0)
    best_raw = np.take_along_axis(raw_stack, winner[None, ...], axis=0)[0]
    dvals = np.asarray(scale_diams, dtype=np.float32)
    best_diam = dvals[winner]

    # Weak auxiliary evidence: local darkness at approximately one via scale. This is
    # intentionally low weight, unlike V20, so a tiny darkest core cannot dominate.
    bg_sigma = max(3.0, 0.70 * nominal_diam_px)
    local_bg = cv2.GaussianBlur(gs, (0, 0), sigmaX=bg_sigma, sigmaY=bg_sigma)
    local_dark = _robust01_v20(np.maximum(local_bg - gs, 0.0), 20.0, 99.3)
    score = (1.0 - V21_VIA_SCORE_LOCAL_DARK_WEIGHT) * best + V21_VIA_SCORE_LOCAL_DARK_WEIGHT * local_dark
    score = cv2.GaussianBlur(score.astype(np.float32), (0, 0), sigmaX=0.65, sigmaY=0.65)
    score = np.clip(score, 0.0, 1.0)

    lo, hi = np.percentile(g, [5.0, 95.0])
    robust_range = max(1.0, float(hi - lo))
    raw_floor = max(V21_VIA_MIN_RAW_CONTRAST_GRAY, 0.008 * robust_range)
    return score, best_raw, best_diam, {
        "search_diameter_px": float(nominal_diam_px),
        "robust_gray_range": robust_range,
        "raw_contrast_floor": float(raw_floor),
        "scale_ratios": ",".join(f"{x:.2f}" for x in V21_VIA_SCALE_RATIOS),
    }


def _make_circle_desc_v21(center: np.ndarray, diam_px: float, source: str) -> Optional[ShapeDesc]:
    cxy = np.asarray(center, dtype=float)
    r = max(2.5, 0.50 * float(diam_px))
    pts = []
    for a in np.linspace(0.0, 2.0 * math.pi, 48, endpoint=False):
        pts.append([cxy[0] + r * math.cos(a), cxy[1] + r * math.sin(a)])
    contour = np.asarray(np.rint(pts), dtype=np.int32).reshape(-1, 1, 2)
    d = describe_contour(contour)
    if d is None:
        return None
    d.center = cxy.astype(np.float64)
    _mark_source(d, source)
    return d


def _peak_components_v21(mask: np.ndarray, score: np.ndarray) -> List[Tuple[float, int, int]]:
    n, labels = cv2.connectedComponents(mask.astype(np.uint8), connectivity=8)
    out: List[Tuple[float, int, int]] = []
    for lab in range(1, n):
        ys, xs = np.where(labels == lab)
        if len(xs) == 0:
            continue
        vals = score[ys, xs]
        k = int(np.argmax(vals))
        out.append((float(vals[k]), int(xs[k]), int(ys[k])))
    out.sort(key=lambda z: z[0], reverse=True)
    return out


def _via_scale_peaks_v21(
    gray: np.ndarray,
    score: np.ndarray,
    raw: np.ndarray,
    best_diam: np.ndarray,
    nominal_diam_px: float,
    raw_floor: float,
) -> List[ShapeDesc]:
    H, W = gray.shape
    k = int(max(3, round(V21_VIA_PEAK_NMS_DIAM_FRAC * nominal_diam_px)))
    if k % 2 == 0:
        k += 1
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k, k))
    dil = cv2.dilate(score, kernel)
    maxima = (score >= dil - 1e-6) & (score >= V21_VIA_PEAK_SCORE_MIN) & (raw >= raw_floor)
    margin = max(3, int(math.ceil(V21_VIA_BORDER_DIAM_FRAC * nominal_diam_px)))
    maxima[:margin, :] = False; maxima[-margin:, :] = False
    maxima[:, :margin] = False; maxima[:, -margin:] = False
    peaks = _peak_components_v21(maxima, score)

    descs: List[ShapeDesc] = []
    for sc, x, y in peaks[:V21_VIA_MAX_PEAKS]:
        dloc = float(best_diam[y, x]) if np.isfinite(best_diam[y, x]) else nominal_diam_px
        # Detection scale is allowed to vary, but the synthetic candidate itself uses
        # a near-physical via size; it is only a center marker, not the final measurement.
        d = _make_circle_desc_v21(np.array([x, y], dtype=float), dloc, "v21_scale_matched_peak")
        if d is None:
            continue
        setattr(d, "_v21_score", float(sc))
        setattr(d, "_v21_raw", float(raw[y, x]))
        setattr(d, "_v21_detect_diam_px", float(dloc))
        descs.append(d)
    return descs


def _via_legacy_size_candidates_v21(
    pool: List[ShapeDesc],
    score: np.ndarray,
    raw: np.ndarray,
    nominal_diam_px: float,
    raw_floor: float,
    image_shape: Tuple[int, int],
) -> List[ShapeDesc]:
    """Keep useful V18/V19 candidates, but only when their scale is physically plausible."""
    H, W = image_shape[:2]
    out: List[ShapeDesc] = []
    margin = max(3.0, V21_VIA_BORDER_DIAM_FRAC * nominal_diam_px)
    for d in pool:
        x, y = map(float, d.center)
        if x <= margin or y <= margin or x >= W - 1 - margin or y >= H - 1 - margin:
            continue
        # Equivalent diameter is more stable than jagged min/max contour extents.
        pa = max(1.0, float(contour_pixel_count(d.contour)))
        eqd = math.sqrt(4.0 * pa / math.pi)
        if not (V21_VIA_LEGACY_MIN_DIAM_FACTOR * nominal_diam_px <= eqd <= V21_VIA_LEGACY_MAX_DIAM_FACTOR * nominal_diam_px):
            continue
        xi = int(np.clip(round(x), 0, W - 1)); yi = int(np.clip(round(y), 0, H - 1))
        sc = float(score[yi, xi]); rw = float(raw[yi, xi])
        # Old candidates are permitted slightly below the scale-peak score threshold;
        # the physical-size and lattice checks will decide their fate.
        if sc < 0.18 or rw < 0.75 * raw_floor:
            continue
        setattr(d, "_v21_score", sc)
        setattr(d, "_v21_raw", rw)
        setattr(d, "_v21_detect_diam_px", eqd)
        _mark_source(d, "v21_legacy_size_valid")
        out.append(d)
    return out


def _merge_via_center_candidates_v21(items: List[ShapeDesc], nominal_diam_px: float) -> List[ShapeDesc]:
    items = sorted(items, key=lambda d: float(getattr(d, "_v21_score", 0.0)), reverse=True)
    out: List[ShapeDesc] = []
    tol = max(2.5, 0.32 * nominal_diam_px)
    for d in items:
        hit = None
        for i, e in enumerate(out):
            if float(np.linalg.norm(d.center - e.center)) <= tol:
                hit = i
                break
        if hit is None:
            out.append(d)
        elif float(getattr(d, "_v21_score", 0.0)) > float(getattr(out[hit], "_v21_score", 0.0)):
            out[hit] = d
    return out[:V21_VIA_MAX_PEAKS]


def _choose_via_lattice_v21(
    candidates: List[ShapeDesc], nominal_diam_px: float
) -> Tuple[Dict[str, object], List[ShapeDesc]]:
    if len(candidates) < 4:
        return {"used": False, "lattice_type": "none"}, []
    ranked = sorted(candidates, key=lambda d: float(getattr(d, "_v21_score", 0.0)), reverse=True)
    counts = [10, 14, 20, 28, 38, 52, 70, 95, 125, len(ranked)]
    counts = sorted(set(max(4, min(len(ranked), c)) for c in counts))
    best_model: Dict[str, object] = {"used": False, "lattice_type": "none"}
    best_seeds: List[ShapeDesc] = []
    best_score = -np.inf
    for n in counts:
        seeds = ranked[:n]
        model = infer_lattice_v14(seeds, (10_000, 10_000), "via")
        # infer_lattice_v14 only uses image_shape for final site bookkeeping; re-fit with
        # the real shape later is unnecessary because origin/basis come from centers.
        if not model.get("used", False):
            continue
        la = float(model.get("basis_a_len", np.nan)); lb = float(model.get("basis_b_len", np.nan))
        ang = float(model.get("basis_angle_deg", np.nan))
        if not (np.isfinite(la) and np.isfinite(lb) and np.isfinite(ang)):
            continue
        if min(la, lb) < V21_VIA_LATTICE_MIN_PITCH_DIAM_FACTOR * nominal_diam_px:
            continue
        if max(la, lb) > V21_VIA_LATTICE_MAX_PITCH_DIAM_FACTOR * nominal_diam_px:
            continue
        if not (V21_VIA_LATTICE_MIN_BASIS_ANGLE <= ang <= V21_VIA_LATTICE_MAX_BASIS_ANGLE):
            continue
        expected = int(model.get("expected_sites", 0)); supported = int(model.get("supported_sites", 0))
        occ = float(model.get("occupancy", 0.0))
        if expected <= 0 or expected > V21_VIA_LATTICE_MAX_EXPECTED or supported < 5:
            continue
        mean_peak = float(np.mean([getattr(d, "_v21_score", 0.0) for d in seeds])) if seeds else 0.0
        # Prefer near-rectangular, well-populated arrays with many actual scale-matched peaks.
        q = (
            4.4 * occ
            + 0.065 * min(supported, 90)
            - 0.008 * expected
            - 0.018 * abs(ang - 90.0)
            + 0.55 * mean_peak
        )
        if q > best_score:
            best_score = q
            best_model = model
            best_seeds = seeds
    return best_model, best_seeds


def _via_refine_site_center_v21(
    score: np.ndarray,
    raw: np.ndarray,
    site_xy: np.ndarray,
    nominal_diam_px: float,
    nearest_pitch: float,
) -> Tuple[np.ndarray, float, float]:
    H, W = score.shape
    r = V21_VIA_SITE_REFINE_DIAM_FRAC * nominal_diam_px
    if np.isfinite(nearest_pitch) and nearest_pitch > 0:
        r = min(r, 0.22 * nearest_pitch)
    r = max(2.0, r)
    cx, cy = map(float, site_xy)
    x0 = max(0, int(math.floor(cx - r))); x1 = min(W, int(math.ceil(cx + r + 1)))
    y0 = max(0, int(math.floor(cy - r))); y1 = min(H, int(math.ceil(cy + r + 1)))
    if x1 <= x0 or y1 <= y0:
        return np.asarray(site_xy, dtype=float), 0.0, 0.0
    sub = score[y0:y1, x0:x1]
    yy, xx = np.mgrid[y0:y1, x0:x1]
    inside = (xx - cx) ** 2 + (yy - cy) ** 2 <= r ** 2
    if not np.any(inside):
        return np.asarray(site_xy, dtype=float), 0.0, 0.0
    masked = np.where(inside, sub, -np.inf)
    iy, ix = np.unravel_index(int(np.argmax(masked)), masked.shape)
    x = x0 + ix; y = y0 + iy
    return np.array([float(x), float(y)]), float(score[y, x]), float(raw[y, x])


def _sample_axis_profile_v20(
    gray: np.ndarray,
    center: np.ndarray,
    axis: np.ndarray,
    half_span: float,
    perp_half_band: float,
) -> Tuple[np.ndarray, np.ndarray]:
    axis = np.asarray(axis, dtype=np.float64)
    axis /= np.linalg.norm(axis) + 1e-12
    perp = np.array([-axis[1], axis[0]], dtype=np.float64)
    n_t = max(15, int(math.ceil(2.0 * half_span)) + 1)
    t = np.linspace(-half_span, half_span, n_t, dtype=np.float32)
    n_s = max(3, int(2 * math.ceil(perp_half_band) + 1))
    s = np.linspace(-perp_half_band, perp_half_band, n_s, dtype=np.float32)
    tt, ss = np.meshgrid(t, s)
    cx, cy = float(center[0]), float(center[1])
    map_x = (cx + tt * axis[0] + ss * perp[0]).astype(np.float32)
    map_y = (cy + tt * axis[1] + ss * perp[1]).astype(np.float32)
    vals = cv2.remap(
        gray.astype(np.float32), map_x, map_y,
        interpolation=cv2.INTER_LINEAR,
        borderMode=cv2.BORDER_REFLECT_101,
    )
    prof = np.median(vals, axis=0).astype(np.float64)
    if len(prof) >= 5:
        prof = cv2.GaussianBlur(prof.reshape(1, -1).astype(np.float32), (0, 0), sigmaX=0.85).ravel().astype(np.float64)
    return t.astype(np.float64), prof


def _one_side_edge_v21(
    t: np.ndarray,
    prof: np.ndarray,
    side: str,
    nominal_radius: float,
) -> Tuple[Optional[float], Dict[str, float]]:
    """Find one physical via edge, explicitly excluding the inner-core radius."""
    if len(t) < 9:
        return None, {"contrast": 0.0}
    R = max(3.0, float(nominal_radius))
    rmin = V21_VIA_EDGE_MIN_RADIUS_FACTOR * R
    rmax = V21_VIA_EDGE_MAX_RADIUS_FACTOR * R
    core = np.abs(t) <= 0.30 * R
    if np.count_nonzero(core) < 2:
        core = np.abs(t) <= max(1.0, 0.20 * R)
    dark = float(np.percentile(prof[core], 60.0)) if np.any(core) else float(prof[np.argmin(np.abs(t))])

    if side == "left":
        r = -t
        side_mask = t < 0
        signed_grad = -np.gradient(prof, t)
    else:
        r = t
        side_mask = t > 0
        signed_grad = np.gradient(prof, t)

    search = side_mask & (r >= rmin) & (r <= rmax)
    if np.count_nonzero(search) < 3:
        return None, {"contrast": 0.0}

    bgmask = side_mask & (r >= 1.12 * R) & (r <= min(1.70 * R, float(np.max(np.abs(t)))))
    if np.count_nonzero(bgmask) >= 3:
        bg = float(np.median(prof[bgmask]))
    else:
        idx_side = np.where(side_mask)[0]
        take = max(2, int(0.18 * len(idx_side)))
        pick = idx_side[:take] if side == "left" else idx_side[-take:]
        bg = float(np.median(prof[pick])) if len(pick) else float(np.median(prof))
    contrast = bg - dark
    if contrast <= 0.15:
        return None, {"contrast": contrast}

    idx = np.where(search)[0]
    rr = r[idx]
    gg = np.maximum(signed_grad[idx], 0.0)
    gscale = float(np.percentile(np.abs(np.gradient(prof, t)), 90.0)) + 1e-6
    gnorm = np.clip(gg / gscale, 0.0, 3.0)
    prior = np.exp(-0.5 * ((rr - R) / max(1.0, V21_VIA_EDGE_PRIOR_SIGMA_RADIUS_FACTOR * R)) ** 2)
    edge_score = 0.82 * gnorm + 0.18 * prior
    gi = int(idx[int(np.argmax(edge_score))])

    half = dark + 0.50 * contrast
    # Half-gray location only within the physical-radius search interval.
    half_i = int(idx[int(np.argmin(np.abs(prof[idx] - half)))])
    prior_t = -R if side == "left" else R
    edge_t = 0.65 * float(t[gi]) + 0.25 * float(t[half_i]) + 0.10 * float(prior_t)
    return edge_t, {
        "contrast": contrast,
        "gradient": float(signed_grad[gi]),
        "grad_edge": float(t[gi]),
        "halfmid": float(t[half_i]),
        "nominal_radius": R,
    }


def _measure_via_axis_width_v21(
    gray: np.ndarray,
    center: np.ndarray,
    axis: np.ndarray,
    nominal_diam_px: float,
    nearest_pitch: float,
) -> Tuple[float, Dict[str, float]]:
    D = max(6.0, float(nominal_diam_px)); R = 0.5 * D
    half_span = max(1.78 * R, 0.90 * D)
    if np.isfinite(nearest_pitch) and nearest_pitch > 0:
        half_span = min(half_span, 0.47 * nearest_pitch)
        half_span = max(half_span, 1.52 * R)
    perp = max(1.0, V21_VIA_EDGE_PERP_BAND_DIAM_FRAC * D)
    t, prof = _sample_axis_profile_v20(gray, center, axis, half_span, perp)
    left, ml = _one_side_edge_v21(t, prof, "left", R)
    right, mr = _one_side_edge_v21(t, prof, "right", R)
    if left is None or right is None or right <= left:
        return np.nan, {"contrast": min(ml.get("contrast", 0.0), mr.get("contrast", 0.0))}
    width = float(right - left)
    if not (V21_VIA_MEASURE_MIN_DIAM_FACTOR * D <= width <= V21_VIA_MEASURE_MAX_DIAM_FACTOR * D):
        return np.nan, {
            "contrast": min(ml.get("contrast", 0.0), mr.get("contrast", 0.0)),
            "raw_width": width,
        }
    return width, {
        "contrast": min(ml.get("contrast", 0.0), mr.get("contrast", 0.0)),
        "left_edge_t": float(left), "right_edge_t": float(right),
        "left_gradient": float(ml.get("gradient", np.nan)),
        "right_gradient": float(mr.get("gradient", np.nan)),
    }


def _local_halflevel_contour_v21(
    gray: np.ndarray,
    center: np.ndarray,
    nominal_diam_px: float,
) -> Optional[np.ndarray]:
    """Fallback segmentation at the local core/background half-gray level.

    Because the background estimate comes from outside the nominal radius, this fallback
    also avoids the V20 failure mode of contouring only the darkest inner core.
    """
    H, W = gray.shape
    D = max(6.0, float(nominal_diam_px)); R = 0.5 * D
    half = max(5, int(math.ceil(V21_VIA_LOCAL_CONTOUR_ROI_DIAM_FACTOR * D)))
    cx, cy = map(float, center)
    x0 = max(0, int(math.floor(cx - half))); x1 = min(W, int(math.ceil(cx + half + 1)))
    y0 = max(0, int(math.floor(cy - half))); y1 = min(H, int(math.ceil(cy + half + 1)))
    if x1 - x0 < 7 or y1 - y0 < 7:
        return None
    roi = gray[y0:y1, x0:x1].astype(np.float32)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    rr = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
    core = rr <= 0.30 * R
    ring = (rr >= 1.10 * R) & (rr <= 1.62 * R)
    if np.count_nonzero(core) < 4 or np.count_nonzero(ring) < 8:
        return None
    dark = float(np.percentile(gray[yy[core], xx[core]], 60.0))
    bg = float(np.median(gray[yy[ring], xx[ring]]))
    if bg - dark <= 0.15:
        return None
    thr = dark + 0.52 * (bg - dark)
    mask = (roi <= thr).astype(np.uint8) * 255
    k = int(max(1, round(V21_VIA_LOCAL_CONTOUR_CLOSE_FRAC * D)))
    if k % 2 == 0:
        k += 1
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k, k))
    mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, kernel, iterations=1)
    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        return None
    lc = np.array([cx - x0, cy - y0], dtype=float)
    best = None
    for c in contours:
        if len(c) < 5:
            continue
        M = cv2.moments(c)
        if abs(M.get("m00", 0.0)) > 1e-9:
            cc = np.array([M["m10"] / M["m00"], M["m01"] / M["m00"]], dtype=float)
        else:
            pts = c.reshape(-1, 2).astype(float); cc = pts.mean(axis=0)
        dist = float(np.linalg.norm(cc - lc))
        area = float(cv2.contourArea(c))
        if area <= 0:
            continue
        eqd = math.sqrt(4.0 * area / math.pi)
        if not (0.55 * D <= eqd <= 1.65 * D):
            continue
        q = dist / max(D, 1e-9) + 0.08 * abs(eqd / D - 1.0)
        if best is None or q < best[0]:
            best = (q, c)
    if best is None:
        return None
    c = best[1].astype(np.int32).copy()
    c[:, 0, 0] += x0; c[:, 0, 1] += y0
    return c


def _ellipse_contour_v21(center: np.ndarray, x_axis: np.ndarray, y_axis: np.ndarray, w: float, h: float) -> np.ndarray:
    c = np.asarray(center, dtype=float)
    xa = np.asarray(x_axis, dtype=float); xa /= np.linalg.norm(xa) + 1e-12
    ya = np.asarray(y_axis, dtype=float); ya /= np.linalg.norm(ya) + 1e-12
    pts = []
    for a in np.linspace(0.0, 2.0 * math.pi, 64, endpoint=False):
        p = c + 0.5 * w * math.cos(a) * xa + 0.5 * h * math.sin(a) * ya
        pts.append(p)
    return np.asarray(np.rint(pts), dtype=np.int32).reshape(-1, 1, 2)


def _measure_selected_via_centers_v21(
    gray: np.ndarray,
    selected: List[ShapeDesc],
    axis: Dict[str, object],
    nominal_diam_px: float,
    nearest_pitch: float,
) -> List[Tuple[ShapeDesc, float, float, Dict[str, float], Dict[str, float], str, str]]:
    xa = np.asarray(axis["x_axis"], dtype=float)
    ya = np.asarray(axis["y_axis_up"], dtype=float)
    out = []
    for d in selected:
        wx, mx = _measure_via_axis_width_v21(gray, d.center, xa, nominal_diam_px, nearest_pitch)
        hy, my = _measure_via_axis_width_v21(gray, d.center, ya, nominal_diam_px, nearest_pitch)
        xmethod = "size_prior_gradient" if np.isfinite(wx) else ""
        ymethod = "size_prior_gradient" if np.isfinite(hy) else ""
        local_contour = None
        if not np.isfinite(wx) or not np.isfinite(hy):
            local_contour = _local_halflevel_contour_v21(gray, d.center, nominal_diam_px)
        if local_contour is not None:
            if not np.isfinite(wx):
                w2 = _projection_span_v19(local_contour, xa)
                if V21_VIA_MEASURE_MIN_DIAM_FACTOR * nominal_diam_px <= w2 <= V21_VIA_MEASURE_MAX_DIAM_FACTOR * nominal_diam_px:
                    wx = float(w2); xmethod = "local_halflevel_fallback"
            if not np.isfinite(hy):
                h2 = _projection_span_v19(local_contour, ya)
                if V21_VIA_MEASURE_MIN_DIAM_FACTOR * nominal_diam_px <= h2 <= V21_VIA_MEASURE_MAX_DIAM_FACTOR * nominal_diam_px:
                    hy = float(h2); ymethod = "local_halflevel_fallback"
        if not np.isfinite(wx) or not np.isfinite(hy):
            continue
        out.append((d, float(wx), float(hy), mx, my, xmethod, ymethod))
    return out


def _select_and_measure_vias_v19(
    gray: np.ndarray,
    pool: List[ShapeDesc],
    px_nm: float,
    image_name: str,
    stage: str,
    condition: str,
    pattern_key: str,
) -> Tuple[List[dict], List[ShapeDesc], List[dict], Dict[str, object]]:
    nominal_nm, nominal_diam_px = _via_nominal_diam_px_v21(pattern_key, px_nm)
    score, raw, best_diam, resp_meta = _via_scale_response_v21(gray, nominal_diam_px)
    raw_floor = float(resp_meta["raw_contrast_floor"])

    scale_peaks = _via_scale_peaks_v21(gray, score, raw, best_diam, nominal_diam_px, raw_floor)
    legacy = _via_legacy_size_candidates_v21(pool, score, raw, nominal_diam_px, raw_floor, gray.shape)
    candidates = _merge_via_center_candidates_v21(scale_peaks + legacy, nominal_diam_px)
    model, model_seeds = _choose_via_lattice_v21(candidates, nominal_diam_px)

    selected: List[ShapeDesc] = []
    nearest_pitch = np.nan
    recovered = 0
    expected_positions = 0
    if model.get("used", False):
        vals = [float(model.get("basis_a_len", np.nan)), float(model.get("basis_b_len", np.nan))]
        vals = [v for v in vals if np.isfinite(v) and v > 0]
        nearest_pitch = min(vals) if vals else np.nan
        sites = lattice_sites_v14(model, gray.shape, 1, "via")
        expected_positions = len(sites)
        margin = max(3.0, V21_VIA_BORDER_DIAM_FRAC * nominal_diam_px)
        for site in sites:
            xy0 = np.array([float(site["x"]), float(site["y"])], dtype=float)
            if xy0[0] <= margin or xy0[1] <= margin or xy0[0] >= gray.shape[1] - 1 - margin or xy0[1] >= gray.shape[0] - 1 - margin:
                continue
            xy, sc, rw = _via_refine_site_center_v21(score, raw, xy0, nominal_diam_px, nearest_pitch)
            internal = bool(site.get("is_internal", False))
            sc_min = V21_VIA_INTERNAL_SITE_SCORE_MIN if internal else V21_VIA_EXTERNAL_SITE_SCORE_MIN
            rw_min = raw_floor if internal else V21_VIA_EXTERNAL_RAW_FACTOR * raw_floor
            if sc < sc_min or rw < rw_min:
                continue
            d = _make_circle_desc_v21(xy, nominal_diam_px, "v21_lattice_site")
            if d is None:
                continue
            setattr(d, "_v21_score", sc); setattr(d, "_v21_raw", rw)
            setattr(d, "_v21_detect_diam_px", float(best_diam[int(round(xy[1])), int(round(xy[0]))]))
            selected.append(d)
            if not any(float(np.linalg.norm(xy - s.center)) <= 0.30 * nominal_diam_px for s in model_seeds):
                recovered += 1
    else:
        # Conservative fallback: use strong scale-matched peaks only. The physical-size
        # edge measurement below will still reject false tiny dark regions.
        strong = [d for d in candidates if float(getattr(d, "_v21_score", 0.0)) >= 0.42 and float(getattr(d, "_v21_raw", 0.0)) >= raw_floor]
        strong = sorted(strong, key=lambda d: float(getattr(d, "_v21_score", 0.0)), reverse=True)
        selected = strong[:80]
        if len(selected) >= 4:
            comp, _rr, p = _largest_center_component_v19(selected, "via_v21_fallback_array")
            if len(comp) >= 4:
                selected = comp
                nearest_pitch = p

    # Center deduplication, independent of segmentation contour shape.
    selected = _merge_via_center_candidates_v21(selected, nominal_diam_px)

    # Expanded-ring recovery is allowed only when an entire weak row/column has support.
    # A single dark noise spot just outside an otherwise complete array must NOT create
    # a fake extra via. Internal model sites are retained; external lattice levels need
    # at least three supporting objects.
    if model.get("used", False) and len(selected) >= 6:
        try:
            B = np.column_stack([np.asarray(model["a"], dtype=float), np.asarray(model["b"], dtype=float)])
            invB = np.linalg.inv(B)
            o = np.asarray(model["origin"], dtype=float)
            ij = [tuple(np.rint(invB @ (d.center - o)).astype(int)) for d in selected]
            ci = Counter(q[0] for q in ij); cj = Counter(q[1] for q in ij)
            imin, imax = int(model["i_min"]), int(model["i_max"])
            jmin, jmax = int(model["j_min"]), int(model["j_max"])
            keep = []
            for d, q in zip(selected, ij):
                i, j = q
                ext_i_ok = (imin <= i <= imax) or (ci[i] >= 3)
                ext_j_ok = (jmin <= j <= jmax) or (cj[j] >= 3)
                if ext_i_ok and ext_j_ok:
                    keep.append(d)
            if len(keep) >= max(4, int(0.70 * len(selected))):
                selected = keep
        except Exception:
            pass

    selected.sort(key=lambda d: (float(d.center[1]), float(d.center[0])))

    # First fit of user-defined via axes from the bottom row.
    axis = _fit_via_bottom_row_axes_v19(selected, gray.shape)
    measured = _measure_selected_via_centers_v21(gray, selected, axis, nominal_diam_px, nearest_pitch)

    # Bad centers usually fail the physical-size edge test. Remove them, refit the bottom
    # row, and measure once more so the final X/Y axes are based only on real vias.
    if len(measured) >= V21_VIA_MIN_FINAL_OBJECTS:
        selected2 = [m[0] for m in measured]
        axis2 = _fit_via_bottom_row_axes_v19(selected2, gray.shape)
        measured2 = _measure_selected_via_centers_v21(gray, selected2, axis2, nominal_diam_px, nearest_pitch)
        if len(measured2) >= V21_VIA_MIN_FINAL_OBJECTS:
            measured = measured2
            axis = axis2

    selected_final = [m[0] for m in measured]
    selected_final.sort(key=lambda d: (float(d.center[1]), float(d.center[0])))
    if selected_final:
        # Refit once after ordering/cleanup, then final measurement in exactly these axes.
        axis = _fit_via_bottom_row_axes_v19(selected_final, gray.shape)
        measured_final = _measure_selected_via_centers_v21(gray, selected_final, axis, nominal_diam_px, nearest_pitch)
        if len(measured_final) >= V21_VIA_MIN_FINAL_OBJECTS:
            measured = measured_final
            selected_final = [m[0] for m in measured]

    # measured and selected_final now have identical order from the last pass.
    axis = _fit_via_bottom_row_axes_v19(selected_final, gray.shape)
    xa = np.asarray(axis["x_axis"], dtype=float)
    ya = np.asarray(axis["y_axis_up"], dtype=float)
    origin = np.asarray(axis["origin"], dtype=float)
    bottom_centers = [np.asarray(selected_final[i].center, dtype=float) for i in axis.get("bottom_indices", []) if 0 <= i < len(selected_final)]

    rows: List[dict] = []
    final_selected: List[ShapeDesc] = []
    for d, wx, hy, mx, my, xmethod, ymethod in measured:
        # Final hard physical sanity check. This is the key anti-tiny-core guard.
        w_nm = wx * px_nm; h_nm = hy * px_nm
        if not (V21_VIA_MEASURE_MIN_DIAM_FACTOR * nominal_nm <= w_nm <= V21_VIA_MEASURE_MAX_DIAM_FACTOR * nominal_nm):
            continue
        if not (V21_VIA_MEASURE_MIN_DIAM_FACTOR * nominal_nm <= h_nm <= V21_VIA_MEASURE_MAX_DIAM_FACTOR * nominal_nm):
            continue
        if not bool(getattr(d, "_v22_actual_boundary", False)):
            # Emergency fallback only. Normal V22 measurements keep the actual irregular contour.
            d.contour = _ellipse_contour_v21(d.center, xa, ya, wx, hy)
        is_bottom = any(float(np.linalg.norm(d.center - c)) <= 0.22 * nominal_diam_px for c in bottom_centers)
        oid = len(final_selected) + 1
        final_selected.append(d)
        rows.append({
            "stage": stage,
            "stage_zh": STAGE_LABELS_ZH_V19[stage],
            "condition": condition,
            "condition_label": CONDITION_LABELS_V19.get(str(condition), f"Condition {condition}"),
            "image": image_name,
            "pattern": pattern_key,
            "pattern_zh": PATTERN_LABELS_ZH_V19[pattern_key],
            "object_id": oid,
            "used_for_statistics": True,
            "is_axis_bottom_row": bool(is_bottom),
            "x_width_nm": float(w_nm),
            "y_height_nm": float(h_nm),
            "tip_gap_y_nm": np.nan,
            "x_width_px": float(wx),
            "y_height_px": float(hy),
            "center_image_x_px": float(d.center[0]),
            "center_image_y_px": float(d.center[1]),
            "center_axis_x_px": float((d.center - origin) @ xa),
            "center_axis_y_up_px": float((d.center - origin) @ ya),
            "circularity": d.circularity,
            "solidity": d.solidity,
            "axis_ratio": d.aspect,
            "selection_source": _source_of(d),
            "pixel_size_nm": px_nm,
            "via_nominal_diameter_nm": nominal_nm,
            "via_nominal_diameter_px": nominal_diam_px,
            "via_scale_match_score": float(getattr(d, "_v21_score", np.nan)),
            "via_scale_raw_contrast_gray": float(getattr(d, "_v21_raw", np.nan)),
            "via_axis_angle_from_image_x_deg": axis["row_angle_deg"],
            "via_axis_fit_rms_px": axis["row_fit_rms_px"],
            "estimated_center_pitch_px": nearest_pitch,
            "x_edge_method": xmethod,
            "y_edge_method": ymethod,
            "x_edge_contrast_gray": mx.get("contrast", np.nan),
            "y_edge_contrast_gray": my.get("contrast", np.nan),
        })

    # Refit axis/bottom flags one last time to final retained objects so annotations and
    # statistics cannot refer to an object removed by the physical-diameter QC.
    if final_selected:
        axis = _fit_via_bottom_row_axes_v19(final_selected, gray.shape)
        origin = np.asarray(axis["origin"], dtype=float)
        xa = np.asarray(axis["x_axis"], dtype=float); ya = np.asarray(axis["y_axis_up"], dtype=float)
        bottom_set = set(axis["bottom_indices"])
        for i, (d, r) in enumerate(zip(final_selected, rows)):
            r["object_id"] = i + 1
            r["is_axis_bottom_row"] = i in bottom_set
            r["center_axis_x_px"] = float((d.center - origin) @ xa)
            r["center_axis_y_up_px"] = float((d.center - origin) @ ya)
            r["via_axis_angle_from_image_x_deg"] = axis["row_angle_deg"]
            r["via_axis_fit_rms_px"] = axis["row_fit_rms_px"]

    rejects: List[dict] = []
    final_centers = [d.center for d in final_selected]
    # Only log scale-plausible center peaks; do not flood the output with every SEM texture.
    for d in candidates:
        if any(float(np.linalg.norm(d.center - c)) <= 0.35 * nominal_diam_px for c in final_centers):
            continue
        if float(getattr(d, "_v21_score", 0.0)) < 0.35:
            continue
        rejects.append(_make_reject_from_contour(
            d.contour,
            "via_v21_physical_size_grid",
            "via_candidate_not_final",
            "scale-matched dark-hole candidate was not retained by rectangular lattice and physical-diameter edge QC",
            desc=d,
            source=_source_of(d),
            scale_match_score=float(getattr(d, "_v21_score", np.nan)),
            raw_annulus_inner_contrast=float(getattr(d, "_v21_raw", np.nan)),
        ))
        if len(rejects) >= 80:
            break

    meta = {
        "array_assist_used": bool(model.get("used", False)),
        "array_lattice_type": model.get("lattice_type", "none"),
        "array_expected_positions": int(expected_positions),
        "array_recovered_count": int(recovered),
        "array_occupancy_seed": float(model.get("occupancy", np.nan)),
        "array_basis_a_length_px": float(model.get("basis_a_len", np.nan)),
        "array_basis_b_length_px": float(model.get("basis_b_len", np.nan)),
        "array_basis_angle_deg": float(model.get("basis_angle_deg", np.nan)),
        "via_bottom_row_count": int(len(axis.get("bottom_indices", []))),
        "via_axis_angle_deg": float(axis.get("row_angle_deg", np.nan)),
        "via_axis_fit_rms_px": float(axis.get("row_fit_rms_px", np.nan)),
        "via_axis_fallback": bool(axis.get("fallback", False)),
        "estimated_center_pitch_px": float(nearest_pitch) if np.isfinite(nearest_pitch) else np.nan,
        "via_v21_nominal_diameter_nm": nominal_nm,
        "via_v21_nominal_diameter_px": nominal_diam_px,
        "via_v21_scale_peak_count": len(scale_peaks),
        "via_v21_legacy_size_candidate_count": len(legacy),
        "via_v21_merged_candidate_count": len(candidates),
        "via_v21_final_count": len(final_selected),
        **resp_meta,
    }
    return rows, final_selected, rejects, {"axis": axis, "meta": meta, "lattice_model": model}


def _annotate_vias_v19(
    gray: np.ndarray,
    selected: List[ShapeDesc],
    rows: List[dict],
    axis_bundle: Dict[str, object],
) -> np.ndarray:
    out = to_color(gray)
    axis = axis_bundle["axis"]
    bottom = set(axis.get("bottom_indices", []))
    for i, d in enumerate(selected):
        color = (0, 255, 255) if i in bottom else (0, 210, 0)
        # V22 contour is normally the ACTUAL locally segmented irregular via boundary.
        cv2.drawContours(out, [d.contour.astype(np.int32)], -1, color, 2 if i in bottom else 1)
        c = tuple(np.round(d.center).astype(int))
        cv2.circle(out, c, 2, color, -1)
        cv2.putText(out, f"V{i+1}", (c[0] + 3, c[1] - 3), cv2.FONT_HERSHEY_SIMPLEX, 0.36, color, 1, cv2.LINE_AA)

    if selected:
        origin = np.asarray(axis["origin"], dtype=float)
        xa = np.asarray(axis["x_axis"], dtype=float)
        ya = np.asarray(axis["y_axis_up"], dtype=float)
        L = max(24.0, 0.16 * min(gray.shape[:2]))
        cv2.arrowedLine(out, p2i(origin), p2i(origin + L * xa), (255, 255, 0), 2, cv2.LINE_AA, tipLength=0.10)
        cv2.arrowedLine(out, p2i(origin), p2i(origin + L * ya), (255, 255, 0), 2, cv2.LINE_AA, tipLength=0.10)
        cv2.putText(out, "X", p2i(origin + L * xa + np.array([4, 0])), cv2.FONT_HERSHEY_SIMPLEX, 0.52, (255, 255, 255), 2, cv2.LINE_AA)
        cv2.putText(out, "Y", p2i(origin + L * ya + np.array([4, 0])), cv2.FONT_HERSHEY_SIMPLEX, 0.52, (255, 255, 255), 2, cv2.LINE_AA)

        mean_w = float(np.mean([r["x_width_nm"] for r in rows])) if rows else np.nan
        mean_h = float(np.mean([r["y_height_nm"] for r in rows])) if rows else np.nan
        nominal = float(rows[0].get("via_nominal_diameter_nm", np.nan)) if rows else np.nan
        text = (
            f"used={len(rows)} bottom-row={len(bottom)} nominal={nominal:.0f}nm  "
            f"Xmean={mean_w:.2f}nm Ymean={mean_h:.2f}nm"
        )
        cv2.rectangle(out, (4, 4), (min(gray.shape[1] - 4, 760), 29), (0, 0, 0), -1)
        cv2.putText(out, text, (8, 22), cv2.FONT_HERSHEY_SIMPLEX, 0.47, (255, 255, 255), 1, cv2.LINE_AA)
    return out



# -----------------------------
# V20 TRENCH: PROFILE DETECTOR FOR CLIPPED OPPOSING DARK BARS
# -----------------------------
# The upper and lower trenches are expected to be cut off by the image top/bottom.
# Therefore touching the image edge is NOT an error for trench160. Recognition is based
# on repeated dark vertical columns first; only then are the two facing tips measured.

V20_TRENCH_EXPECTED_GROUPS_HINT = 9
V20_TRENCH_MIN_GROUPS = 4
V20_TRENCH_MAX_GROUPS = 15
V20_TRENCH_TOP_BODY_RANGE = (0.04, 0.39)
V20_TRENCH_BOTTOM_BODY_RANGE = (0.61, 0.96)
V20_TRENCH_TOP_TIP_SEARCH = (0.18, 0.62)
V20_TRENCH_BOTTOM_TIP_SEARCH = (0.38, 0.82)
V20_TRENCH_X_NMS_MIN_FRAC = 0.035
V20_TRENCH_X_NMS_MAX_FRAC = 0.075
V20_TRENCH_X_PROFILE_SMOOTH_SIGMA_FRAC = 0.003
V20_TRENCH_MIN_PAIR_DARK_CONTRAST = 0.7
V20_TRENCH_TIP_REFINE_HALF_FRAC = 0.10
V20_TRENCH_TIP_BODY_WINDOW_FRAC = 0.055
V20_TRENCH_GAP_MAD_FACTOR = 7.0
V20_TRENCH_GAP_FRAC_TOL = 0.80
V20_TRENCH_WIDTH_MAD_FACTOR = 7.0
V20_TRENCH_WIDTH_FRAC_TOL = 0.85
V20_TRENCH_DRAW_BODY_LIMIT_FRAC = 0.03


def _trench_vertical_dark_score_v20(gray: np.ndarray) -> np.ndarray:
    """Dark vertical-line score robust to a globally dark background."""
    g = gray.astype(np.float32)
    ge = cv2.createCLAHE(clipLimit=1.8, tileGridSize=(8, 8)).apply(gray).astype(np.float32)
    H, W = gray.shape
    local = []
    # Strongly blur in X, weakly in Y: vertical dark bars remain dark relative to their
    # horizontal neighborhood while slow SEM background shading is removed.
    for fx in (0.012, 0.022, 0.038):
        sx = max(4.0, fx * W)
        sy = max(1.2, 0.003 * H)
        bg = cv2.GaussianBlur(ge, (0, 0), sigmaX=sx, sigmaY=sy)
        local.append(np.maximum(bg - ge, 0.0))
    local_dark = _robust01_v20(np.maximum.reduce(local), 3.0, 99.0)
    q65 = float(np.percentile(g, 65.0)); q12 = float(np.percentile(g, 12.0))
    abs_dark = _robust01_v20(np.clip((q65 - g) / max(4.0, q65 - q12), 0.0, 2.5), 2.0, 99.0)
    score = 0.76 * local_dark + 0.24 * abs_dark
    return np.clip(score, 0.0, 1.0)


def _smooth_1d_v20(x: np.ndarray, sigma: float) -> np.ndarray:
    a = np.asarray(x, dtype=np.float32).reshape(1, -1)
    if a.shape[1] < 3:
        return a.ravel().astype(np.float64)
    sigma = max(0.6, float(sigma))
    return cv2.GaussianBlur(a, (0, 0), sigmaX=sigma).ravel().astype(np.float64)


def _nms_peaks_1d_v20(profile: np.ndarray, min_dist: int, max_peaks: int, threshold: float) -> List[int]:
    p = np.asarray(profile, dtype=float)
    if len(p) < 3:
        return []
    local = np.where((p[1:-1] >= p[:-2]) & (p[1:-1] >= p[2:]))[0] + 1
    local = [int(i) for i in local if p[i] >= threshold]
    local.sort(key=lambda i: p[i], reverse=True)
    keep: List[int] = []
    for i in local:
        if all(abs(i - j) >= min_dist for j in keep):
            keep.append(i)
            if len(keep) >= max_peaks:
                break
    keep.sort()
    return keep


def _regularize_trench_peaks_v20(peaks: List[int], profile: np.ndarray, image_w: int) -> List[int]:
    """Keep/recover a roughly equally spaced trench series without forcing exactly 9.

    Border texture is removed first. Existing NMS peaks are trusted more than a
    theoretical grid; predicted sites are only added when the X profile itself is dark.
    """
    if not peaks:
        return []
    p = np.asarray(profile, dtype=float)
    border = max(3, int(round(0.015 * image_w)))
    xs0 = [int(x) for x in sorted(peaks) if border <= x <= image_w - 1 - border]
    if len(xs0) <= 3:
        return xs0

    diffs0 = np.diff(np.asarray(xs0, dtype=float))
    plausible = diffs0[(diffs0 >= 0.035 * image_w) & (diffs0 <= 0.22 * image_w)]
    if len(plausible) == 0:
        return xs0[:V20_TRENCH_MAX_GROUPS]

    # Missing one bar generates ~2*pitch. Estimate the fundamental pitch from the
    # lower/middle spacing population and then refine from inlier adjacent gaps.
    q45 = float(np.percentile(plausible, 45.0))
    fundamental = plausible[plausible <= 1.45 * q45]
    pitch = float(np.median(fundamental)) if len(fundamental) else q45
    if not np.isfinite(pitch) or pitch <= 2:
        return xs0[:V20_TRENCH_MAX_GROUPS]

    # Pick phase by maximum support from real peaks.
    best = None
    xs = np.asarray(xs0, dtype=float)
    for anchor in xs:
        phase = float(anchor % pitch)
        kk = np.rint((xs - phase) / pitch)
        pred = phase + kk * pitch
        res = np.abs(xs - pred)
        good = res <= 0.28 * pitch
        support = int(np.count_nonzero(good))
        medres = float(np.median(res[good])) if support else 1e9
        score = support - 0.20 * medres / max(pitch, 1e-9)
        if best is None or score > best[0]:
            best = (score, phase, good, kk)
    if best is None:
        return xs0[:V20_TRENCH_MAX_GROUPS]

    _, phase, good, kk = best
    if int(np.count_nonzero(good)) >= 3:
        # Refine phase and pitch by least-squares on matched integer indices.
        kval = kk[good].astype(float)
        xval = xs[good]
        A = np.column_stack([np.ones_like(kval), kval])
        coef, *_ = np.linalg.lstsq(A, xval, rcond=None)
        phase_fit, pitch_fit = float(coef[0]), float(coef[1])
        if 0.70 * pitch <= abs(pitch_fit) <= 1.30 * pitch:
            phase, pitch = phase_fit, abs(pitch_fit)

    # Match theoretical sites only over the observed X span plus at most one site on
    # either side. Existing real peaks always win. Missing sites require strong evidence.
    kmin = int(math.floor((min(xs0) - phase) / pitch)) - 1
    kmax = int(math.ceil((max(xs0) - phase) / pitch)) + 1
    site_threshold = max(float(np.percentile(p, 52.0)), float(np.median(p) + 0.18 * np.std(p)))
    out: List[int] = []
    for k in range(kmin, kmax + 1):
        s = phase + k * pitch
        if s < border or s > image_w - 1 - border:
            continue
        near = [(abs(x - s), x) for x in xs0 if abs(x - s) <= 0.30 * pitch]
        if near:
            x_use = int(min(near)[1])
        else:
            c = int(round(s)); r = max(2, int(round(0.20 * pitch)))
            a = max(border, c-r); b = min(image_w-1-border, c+r)
            if b <= a:
                continue
            j = int(a + np.argmax(p[a:b+1]))
            if p[j] < site_threshold:
                continue
            x_use = j
        if not out or abs(x_use - out[-1]) >= max(3, int(0.48 * pitch)):
            out.append(x_use)
    if len(out) > V20_TRENCH_MAX_GROUPS:
        # Prefer a contiguous high-score window rather than scattered noise peaks.
        best_window = None
        for i in range(len(out) - V20_TRENCH_MAX_GROUPS + 1):
            win = out[i:i+V20_TRENCH_MAX_GROUPS]
            sc = float(np.sum([p[x] for x in win]))
            if best_window is None or sc > best_window[0]:
                best_window = (sc, win)
        out = best_window[1] if best_window else out[:V20_TRENCH_MAX_GROUPS]
    return out


def _detect_trench_x_centers_v20(gray: np.ndarray, score: np.ndarray) -> Tuple[List[int], np.ndarray, Dict[str, object]]:
    H, W = gray.shape
    ta, tb = [int(round(v * H)) for v in V20_TRENCH_TOP_BODY_RANGE]
    ba, bb = [int(round(v * H)) for v in V20_TRENCH_BOTTOM_BODY_RANGE]
    ta, tb = max(0, ta), min(H, max(ta + 2, tb))
    ba, bb = max(0, ba), min(H, max(ba + 2, bb))

    top = np.mean(score[ta:tb], axis=0)
    bottom = np.mean(score[ba:bb], axis=0)
    # Geometric mean demands evidence in BOTH truncated bodies at the same X.
    pair = np.sqrt(np.maximum(top, 0.0) * np.maximum(bottom, 0.0))
    pair += 0.12 * np.minimum(top, bottom)
    sigma = max(0.8, V20_TRENCH_X_PROFILE_SMOOTH_SIGMA_FRAC * W)
    pair = _smooth_1d_v20(pair, sigma)

    med = float(np.median(pair)); mad = float(np.median(np.abs(pair - med)))
    robust_sigma = 1.4826 * mad
    thresholds = [
        max(float(np.percentile(pair, 74.0)), med + 0.75 * robust_sigma),
        max(float(np.percentile(pair, 66.0)), med + 0.40 * robust_sigma),
        float(np.percentile(pair, 58.0)),
    ]
    min_dist = int(np.clip(
        0.48 * W / max(V20_TRENCH_EXPECTED_GROUPS_HINT, 1),
        V20_TRENCH_X_NMS_MIN_FRAC * W,
        V20_TRENCH_X_NMS_MAX_FRAC * W,
    ))
    peaks: List[int] = []
    used_thr = thresholds[-1]
    for thr in thresholds:
        pks = _nms_peaks_1d_v20(pair, max(4, min_dist), V20_TRENCH_MAX_GROUPS + 5, thr)
        if len(pks) >= V20_TRENCH_MIN_GROUPS:
            peaks = pks
            used_thr = thr
            break
        if len(pks) > len(peaks):
            peaks = pks
            used_thr = thr
    peaks = _regularize_trench_peaks_v20(peaks, pair, W)
    return peaks, pair, {
        "trench_v20_x_peak_threshold": used_thr,
        "trench_v20_x_candidates": len(peaks),
        "trench_v20_top_body_y0": ta,
        "trench_v20_top_body_y1": tb,
        "trench_v20_bottom_body_y0": ba,
        "trench_v20_bottom_body_y1": bb,
    }


def _trench_x_edges_v20(gray: np.ndarray, x: int, approx_pitch: float) -> Tuple[float, float, float, Dict[str, float]]:
    """Find left/right trench edges from dark body + strongest horizontal transition."""
    H, W = gray.shape
    ta, tb = [int(round(v * H)) for v in V20_TRENCH_TOP_BODY_RANGE]
    ba, bb = [int(round(v * H)) for v in V20_TRENCH_BOTTOM_BODY_RANGE]
    ta, tb = max(0, ta), min(H, max(ta + 2, tb))
    ba, bb = max(0, ba), min(H, max(ba + 2, bb))
    # Median intensity over the two long bodies suppresses isolated SEM texture.
    ptop = np.median(gray[ta:tb].astype(np.float32), axis=0)
    pbot = np.median(gray[ba:bb].astype(np.float32), axis=0)
    prof = 0.5 * (ptop + pbot)
    prof = _smooth_1d_v20(prof, max(0.7, 0.0015 * W))

    half = int(max(6, min(0.42 * approx_pitch if np.isfinite(approx_pitch) else 0.08 * W, 0.10 * W)))
    a = max(1, x - half); b = min(W - 2, x + half)
    if b - a < 6:
        return float(x - 2), float(x + 2), 4.0, {}
    local = prof[a:b+1]
    ci = int(np.clip(x - a, 1, len(local) - 2))
    outer = np.r_[local[:max(2, len(local)//5)], local[-max(2, len(local)//5):]]
    core_r = max(1, int(round(0.10 * len(local))))
    core = local[max(0, ci-core_r):min(len(local), ci+core_r+1)]
    dark = float(np.percentile(core, 30.0)); bg = float(np.median(outer)); contrast = bg - dark
    level = dark + 0.50 * max(contrast, 1.0)
    grad = np.gradient(local)

    left_idx = np.arange(1, ci)
    right_idx = np.arange(ci+1, len(local)-1)
    if len(left_idx) == 0 or len(right_idx) == 0:
        return float(x - 2), float(x + 2), 4.0, {"contrast": contrast}
    lc = left_idx[local[left_idx] <= level]
    rc = right_idx[local[right_idx] >= level]
    lh = int(lc[0]) if len(lc) else int(left_idx[np.argmin(grad[left_idx])])
    rh = int(rc[0]) if len(rc) else int(right_idx[np.argmax(grad[right_idx])])
    sr = max(2, int(round(0.12 * (b - a))))
    ls = np.arange(max(1, lh-sr), min(ci-1, lh+sr)+1)
    rs = np.arange(max(ci+1, rh-sr), min(len(local)-2, rh+sr)+1)
    lg = int(ls[np.argmin(grad[ls])]) if len(ls) else lh
    rg = int(rs[np.argmax(grad[rs])]) if len(rs) else rh
    left = a + 0.62 * lg + 0.38 * lh
    right = a + 0.62 * rg + 0.38 * rh
    if right <= left + 1:
        left, right = float(a + lh), float(a + rh)
    return float(left), float(right), float(max(1.0, right-left)), {
        "contrast": contrast,
        "left_gradient": float(-grad[lg]),
        "right_gradient": float(grad[rg]),
    }


def _trench_vertical_profile_v20(gray: np.ndarray, x0: float, x1: float) -> np.ndarray:
    H, W = gray.shape
    # Use the central 50% of the dark bar so jagged sidewalls do not contaminate tip Y.
    cx = 0.5 * (x0 + x1)
    half = max(1.0, 0.25 * max(2.0, x1 - x0))
    xa = int(np.clip(math.floor(cx - half), 0, W - 1))
    xb = int(np.clip(math.ceil(cx + half), xa + 1, W))
    roi = gray[:, xa:xb].astype(np.float32)
    # 35th percentile retains the dark trench body if one side is slightly blurred.
    prof = np.percentile(roi, 35.0, axis=1)
    return _smooth_1d_v20(prof, max(0.8, 0.0025 * H))


def _best_tip_from_profile_v20(
    prof: np.ndarray,
    facing: str,
    search_range: Tuple[float, float],
    expected_y: Optional[float] = None,
    refine_half: Optional[float] = None,
) -> Tuple[float, Dict[str, float]]:
    H = len(prof)
    lo = max(2, int(round(search_range[0] * H)))
    hi = min(H - 3, int(round(search_range[1] * H)))
    if expected_y is not None and refine_half is not None:
        lo = max(lo, int(round(expected_y - refine_half)))
        hi = min(hi, int(round(expected_y + refine_half)))
    if hi <= lo:
        return np.nan, {}
    grad = np.gradient(prof)
    noise = 1.4826 * float(np.median(np.abs(grad - np.median(grad)))) + 1e-6
    win = max(3, int(round(V20_TRENCH_TIP_BODY_WINDOW_FRAC * H)))
    best = None
    for y in range(lo, hi + 1):
        a0 = max(0, y - win); a1 = y
        b0 = y + 1; b1 = min(H, y + 1 + win)
        if a1 - a0 < 2 or b1 - b0 < 2:
            continue
        above = float(np.median(prof[a0:a1]))
        below = float(np.median(prof[b0:b1]))
        if facing == "down":
            # Upper trench: dark above, bright below -> positive intensity transition.
            contrast = below - above
            edge = float(grad[y])
        else:
            # Lower trench: bright above, dark below -> negative transition.
            contrast = above - below
            edge = float(-grad[y])
        # Both absolute contrast across the boundary and sharp local derivative matter.
        score = 0.58 * max(0.0, edge) / noise + 0.42 * max(0.0, contrast) / (noise * max(2.0, math.sqrt(win)))
        # Mild preference for tips near the image middle, never enough to beat a real edge.
        score -= 0.10 * abs(y - 0.5 * H) / max(1.0, H)
        rec = (score, contrast, edge, y, above, below)
        if best is None or rec[0] > best[0]:
            best = rec
    if best is None:
        return np.nan, {}
    return float(best[3]), {
        "tip_score": float(best[0]),
        "tip_contrast": float(best[1]),
        "tip_gradient": float(best[2]),
        "body_above_gray": float(best[4]),
        "body_below_gray": float(best[5]),
    }


def _detect_trench_pairs_profile_v20(gray: np.ndarray) -> Tuple[List[dict], Dict[str, object]]:
    H, W = gray.shape
    score_map = _trench_vertical_dark_score_v20(gray)
    x_centers, x_profile, meta = _detect_trench_x_centers_v20(gray, score_map)
    if not x_centers:
        meta.update({"trench_v20_pairs_before_qc": 0, "trench_v20_pairs_after_qc": 0})
        return [], meta

    diffs = np.diff(sorted(x_centers))
    approx_pitch = float(np.median(diffs)) if len(diffs) else float(W / max(4, V20_TRENCH_EXPECTED_GROUPS_HINT))
    items = []
    for x in x_centers:
        xl, xr, width, xmeta = _trench_x_edges_v20(gray, int(x), approx_pitch)
        prof = _trench_vertical_profile_v20(gray, xl, xr)
        ty, tm = _best_tip_from_profile_v20(prof, "down", V20_TRENCH_TOP_TIP_SEARCH)
        by, bm = _best_tip_from_profile_v20(prof, "up", V20_TRENCH_BOTTOM_TIP_SEARCH)
        if not np.isfinite(ty) or not np.isfinite(by) or by <= ty + 1:
            continue
        items.append({
            "common_x": float(x),
            "x_left_px": float(xl),
            "x_right_px": float(xr),
            "width_px": float(width),
            "top_tip_y_px": float(ty),
            "bottom_tip_y_px": float(by),
            "gap_px": float(by - ty),
            "x_profile_score": float(x_profile[int(np.clip(x, 0, W-1))]),
            "top_tip_score": float(tm.get("tip_score", np.nan)),
            "bottom_tip_score": float(bm.get("tip_score", np.nan)),
            "top_tip_contrast": float(tm.get("tip_contrast", np.nan)),
            "bottom_tip_contrast": float(bm.get("tip_contrast", np.nan)),
            "top_tip_gradient": float(tm.get("tip_gradient", np.nan)),
            "bottom_tip_gradient": float(bm.get("tip_gradient", np.nan)),
            "x_edge_contrast": float(xmeta.get("contrast", np.nan)),
        })

    # First-pass global tip rows make weak/fuzzy individual pairs much more stable.
    if len(items) >= 3:
        global_top = float(np.median([z["top_tip_y_px"] for z in items]))
        global_bottom = float(np.median([z["bottom_tip_y_px"] for z in items]))
        refine_half = V20_TRENCH_TIP_REFINE_HALF_FRAC * H
        for z in items:
            prof = _trench_vertical_profile_v20(gray, z["x_left_px"], z["x_right_px"])
            ty, tm = _best_tip_from_profile_v20(
                prof, "down", V20_TRENCH_TOP_TIP_SEARCH,
                expected_y=global_top, refine_half=refine_half,
            )
            by, bm = _best_tip_from_profile_v20(
                prof, "up", V20_TRENCH_BOTTOM_TIP_SEARCH,
                expected_y=global_bottom, refine_half=refine_half,
            )
            if np.isfinite(ty) and np.isfinite(by) and by > ty + 1:
                z["top_tip_y_px"] = float(ty)
                z["bottom_tip_y_px"] = float(by)
                z["gap_px"] = float(by - ty)
                z["top_tip_score"] = float(tm.get("tip_score", np.nan))
                z["bottom_tip_score"] = float(bm.get("tip_score", np.nan))
                z["top_tip_contrast"] = float(tm.get("tip_contrast", np.nan))
                z["bottom_tip_contrast"] = float(bm.get("tip_contrast", np.nan))
                z["top_tip_gradient"] = float(tm.get("tip_gradient", np.nan))
                z["bottom_tip_gradient"] = float(bm.get("tip_gradient", np.nan))

    n_before = len(items)
    # Keep weak but repeated pairs; only reject gross population outliers. This is
    # intentionally much looser than the old connected-component geometry filters.
    if len(items) >= 4:
        widths = np.asarray([z["width_px"] for z in items], dtype=float)
        gaps = np.asarray([z["gap_px"] for z in items], dtype=float)
        mw = float(np.median(widths)); mg = float(np.median(gaps))
        madw = 1.4826 * float(np.median(np.abs(widths - mw)))
        madg = 1.4826 * float(np.median(np.abs(gaps - mg)))
        tw = max(V20_TRENCH_WIDTH_FRAC_TOL * mw, V20_TRENCH_WIDTH_MAD_FACTOR * madw, 2.5)
        tg = max(V20_TRENCH_GAP_FRAC_TOL * mg, V20_TRENCH_GAP_MAD_FACTOR * madg, 2.5)
        kept = []
        for z in items:
            if abs(z["width_px"] - mw) > tw:
                continue
            if abs(z["gap_px"] - mg) > tg:
                continue
            # At least one of the two tip transitions must show a positive contrast;
            # the other may be fuzzy in a dark after image.
            tc = max(z.get("top_tip_contrast", 0.0), z.get("bottom_tip_contrast", 0.0))
            if np.isfinite(tc) and tc < V20_TRENCH_MIN_PAIR_DARK_CONTRAST:
                continue
            kept.append(z)
        if len(kept) >= max(V20_TRENCH_MIN_GROUPS, int(math.ceil(0.55 * len(items)))):
            items = kept

    items.sort(key=lambda z: z["common_x"])
    # Fill compatibility fields expected by CSV/annotation code.
    for z in items:
        z["top_width_px"] = z["width_px"]
        z["bottom_width_px"] = z["width_px"]
        z["x_offset_px"] = 0.0
        z["x_overlap_frac"] = 1.0
        z["width_mismatch_frac"] = 0.0
        z["score"] = -float(z.get("x_profile_score", 0.0))

    meta.update({
        "trench_v20_detector": "dark_vertical_profile_plus_tip_gradient",
        "trench_v20_pairs_before_qc": n_before,
        "trench_v20_pairs_after_qc": len(items),
        "trench_v20_estimated_x_pitch_px": approx_pitch,
        "trench_v20_global_top_tip_y_px": float(np.median([z["top_tip_y_px"] for z in items])) if items else np.nan,
        "trench_v20_global_bottom_tip_y_px": float(np.median([z["bottom_tip_y_px"] for z in items])) if items else np.nan,
    })
    return items, meta


def _select_and_measure_trenches_v19(
    gray: np.ndarray,
    pool: List[ShapeDesc],
    px_nm: float,
    image_name: str,
    stage: str,
    condition: str,
    pattern_key: str,
) -> Tuple[List[dict], List[dict], List[dict], Dict[str, object]]:
    # IMPORTANT: pool is intentionally ignored. Opposing trenches are often clipped at
    # the top/bottom image boundaries, so connected-component completeness is the wrong
    # recognition primitive. Detect the repeated dark bars directly from gray profiles.
    pairs, meta = _detect_trench_pairs_profile_v20(gray)
    rows: List[dict] = []
    for oid, p in enumerate(pairs, 1):
        rows.append({
            "stage": stage,
            "stage_zh": STAGE_LABELS_ZH_V19[stage],
            "condition": condition,
            "condition_label": CONDITION_LABELS_V19.get(str(condition), f"Condition {condition}"),
            "image": image_name,
            "pattern": pattern_key,
            "pattern_zh": PATTERN_LABELS_ZH_V19[pattern_key],
            "object_id": oid,
            "used_for_statistics": True,
            "is_axis_bottom_row": False,
            "x_width_nm": p["width_px"] * px_nm,
            "y_height_nm": np.nan,
            "tip_gap_y_nm": p["gap_px"] * px_nm,
            "x_width_px": p["width_px"],
            "y_height_px": np.nan,
            "center_image_x_px": p["common_x"],
            "center_image_y_px": 0.5 * (p["top_tip_y_px"] + p["bottom_tip_y_px"]),
            "center_axis_x_px": p["common_x"],
            "center_axis_y_up_px": np.nan,
            "top_trench_width_nm": p["top_width_px"] * px_nm,
            "bottom_trench_width_nm": p["bottom_width_px"] * px_nm,
            "pair_x_offset_nm": 0.0,
            "pair_x_overlap_frac": 1.0,
            "width_mismatch_frac": 0.0,
            "pair_score": p["score"],
            "pixel_size_nm": px_nm,
            "trench_detector": "dark_profile+gradient",
            "top_tip_contrast_gray": p.get("top_tip_contrast", np.nan),
            "bottom_tip_contrast_gray": p.get("bottom_tip_contrast", np.nan),
            "top_tip_gradient": p.get("top_tip_gradient", np.nan),
            "bottom_tip_gradient": p.get("bottom_tip_gradient", np.nan),
        })
    # No contour-shape rejects are generated for trench160; drawing old Otsu rejects
    # would obscure the useful profile result in low-contrast images.
    rejects: List[dict] = []
    return rows, pairs, rejects, meta


def _annotate_trenches_v19(
    gray: np.ndarray,
    pairs: List[dict],
    rows: List[dict],
) -> np.ndarray:
    out = to_color(gray)
    H, W = gray.shape
    body_margin = int(round(V20_TRENCH_DRAW_BODY_LIMIT_FRAC * H))
    for i, (p, r) in enumerate(zip(pairs, rows), 1):
        xl = int(round(p["x_left_px"])); xr = int(round(p["x_right_px"]))
        ty = int(round(p["top_tip_y_px"])); by = int(round(p["bottom_tip_y_px"]))
        x = int(round(p["common_x"]))
        # Upper/lower bodies are intentionally allowed to touch the image edge.
        cv2.rectangle(out, (max(0, xl), max(0, body_margin)), (min(W-1, xr), max(body_margin, ty)), (0, 210, 0), 1)
        cv2.rectangle(out, (max(0, xl), min(H-1, by)), (min(W-1, xr), min(H-1, H-1-body_margin)), (0, 210, 0), 1)
        p1 = (x, ty); p2 = (x, by)
        cv2.arrowedLine(out, p1, p2, (255, 255, 255), 1, cv2.LINE_AA, tipLength=0.16)
        cv2.arrowedLine(out, p2, p1, (255, 255, 255), 1, cv2.LINE_AA, tipLength=0.16)
        cv2.circle(out, p1, 2, (255, 255, 255), -1)
        cv2.circle(out, p2, 2, (255, 255, 255), -1)
        label = f"T{i}"
        org = (min(W-35, x+3), int(round(0.5*(ty+by))))
        cv2.putText(out, label, org, cv2.FONT_HERSHEY_SIMPLEX, 0.40, (0, 0, 0), 3, cv2.LINE_AA)
        cv2.putText(out, label, org, cv2.FONT_HERSHEY_SIMPLEX, 0.40, (255, 255, 255), 1, cv2.LINE_AA)
    if rows:
        mean_gap = float(np.mean([r["tip_gap_y_nm"] for r in rows]))
        text = f"PROFILE trench pairs={len(rows)}  mean Y tip-gap={mean_gap:.2f} nm"
    else:
        text = "PROFILE trench: no valid opposing dark-bar pairs"
    cv2.rectangle(out, (4, 4), (min(W - 4, 570), 29), (0, 0, 0), -1)
    cv2.putText(out, text, (8, 22), cv2.FONT_HERSHEY_SIMPLEX, 0.48, (255, 255, 255), 1, cv2.LINE_AA)
    return out




# -----------------------------
# V19 SLOT BOTTOM-ROW WRAPPER
# -----------------------------
def _generic_bottom_slot_row_v19(items: List[ShapeDesc]) -> List[ShapeDesc]:
    if not items:
        return []
    heights = [max(1.0, _contour_xy_bounds(d)[3] - _contour_xy_bounds(d)[2]) for d in items]
    tol = max(4.0, 0.62 * float(np.median(heights)))
    ordered = sorted(items, key=lambda d: float(d.center[1]))
    clusters: List[List[ShapeDesc]] = []
    for d in ordered:
        if not clusters:
            clusters.append([d])
            continue
        med_y = float(np.median([z.center[1] for z in clusters[-1]]))
        if abs(float(d.center[1]) - med_y) <= tol:
            clusters[-1].append(d)
        else:
            clusters.append([d])
    return max(clusters, key=lambda c: float(np.median([d.center[1] for d in c])))


def _select_and_measure_bottom_slots_v19(
    gray: np.ndarray,
    pool: List[ShapeDesc],
    px_nm: float,
    image_name: str,
    stage: str,
    condition: str,
    pattern_key: str,
) -> Tuple[List[dict], List[ShapeDesc], List[ShapeDesc], List[dict], Dict[str, object]]:
    old_rows, selected_raw, _, rejects, meta = measure_slot_v18(
        gray, pool, px_nm, image_name, stage
    )
    paired = []
    for d, row in zip(selected_raw, old_rows):
        if _bbox_complete_v19(d, gray.shape, 0.10):
            paired.append((d, row))
        else:
            rejects.append(_make_reject_from_contour(
                d.contour,
                "slot_v19_border_qc",
                "slot_near_or_crossing_image_edge",
                "V18-selected slot is too close to the image edge and may be incomplete",
                desc=d,
            ))
    selected_all = [d for d, _ in paired]

    if paired:
        row_ids = sorted(set(int(row.get("slot_row", 0)) for _, row in paired if int(row.get("slot_row", 0)) > 0))
        if row_ids:
            med_y = {
                rid: float(np.median([d.center[1] for d, row in paired if int(row.get("slot_row", 0)) == rid]))
                for rid in row_ids
            }
            bottom_row_id = max(med_y, key=med_y.get)
            bottom = [d for d, row in paired if int(row.get("slot_row", 0)) == bottom_row_id]
        else:
            bottom_row_id = np.nan
            bottom = _generic_bottom_slot_row_v19([d for d, _ in paired])
    else:
        # Generic strict fallback if the hard four-row model cannot be fitted.
        strict, initial_rejects, _, mean_area, _, _ = _slot_initial_filter_v12(pool)
        rejects.extend(initial_rejects)
        strict = [d for d in _dedup_v19(strict) if _bbox_complete_v19(d, gray.shape, 0.10)]
        selected_all = strict
        bottom = _generic_bottom_slot_row_v19(strict)
        bottom_row_id = np.nan
        paired = []

    # Bottom-row-specific area consistency avoids averaging a merged or tiny outlier.
    bottom, rr = _area_consistency_filter_v19(
        bottom,
        "slot_v19_bottom_row_area_qc",
        lo_factor=0.40,
        hi_factor=2.50,
    )
    rejects.extend(rr)
    bottom = sorted(bottom, key=lambda d: float(d.center[0]))
    bottom_ids = {id(d) for d in bottom}

    rows: List[dict] = []
    for oid, d in enumerate(bottom, 1):
        x_axis = np.array([1.0, 0.0])
        y_axis = np.array([0.0, 1.0])
        w_px = _projection_span_v19(d.contour, x_axis)
        h_px = _projection_span_v19(d.contour, y_axis)
        reg = slot_regularity_metrics(d)
        rows.append({
            "stage": stage,
            "stage_zh": STAGE_LABELS_ZH_V19[stage],
            "condition": condition,
            "condition_label": CONDITION_LABELS_V19.get(str(condition), f"Condition {condition}"),
            "image": image_name,
            "pattern": pattern_key,
            "pattern_zh": PATTERN_LABELS_ZH_V19[pattern_key],
            "object_id": oid,
            "used_for_statistics": True,
            "is_axis_bottom_row": True,
            "x_width_nm": w_px * px_nm,
            "y_height_nm": h_px * px_nm,
            "tip_gap_y_nm": np.nan,
            "x_width_px": w_px,
            "y_height_px": h_px,
            "center_image_x_px": float(d.center[0]),
            "center_image_y_px": float(d.center[1]),
            "center_axis_x_px": float(d.center[0]),
            "center_axis_y_up_px": -float(d.center[1]),
            "slot_row_id_from_v18": bottom_row_id,
            "slot_pixel_area_px": contour_pixel_count(d.contour),
            "slot_regularity_solidity": reg["regularity_solidity"],
            "slot_regularity_convexity": reg["regularity_convexity"],
            "slot_best_template_iou": reg["regularity_template_iou"],
            "selection_source": _source_of(d),
            "pixel_size_nm": px_nm,
        })

    # Explicitly log complete V18 objects in the lowest row that were removed by V19 QC.
    for d, old in paired:
        if int(old.get("slot_row", -999)) == bottom_row_id and id(d) not in bottom_ids:
            if not any(r.get("_desc") is d for r in rejects):
                rejects.append(_make_reject_from_contour(
                    d.contour,
                    "slot_v19_bottom_row_final_qc",
                    "bottom_row_slot_not_used",
                    "V18 bottom-row slot was removed by final border/area consistency QC",
                    desc=d,
                ))
    meta = dict(meta)
    meta.update({
        "slot_bottom_row_id": bottom_row_id,
        "slot_bottom_row_used_count": len(bottom),
        "slot_all_valid_count": len(selected_all),
    })
    return rows, selected_all, bottom, rejects, meta


def _annotate_slots_v19(
    gray: np.ndarray,
    selected_all: List[ShapeDesc],
    bottom: List[ShapeDesc],
    rows: List[dict],
) -> np.ndarray:
    out = to_color(gray)
    bottom_ids = {id(d) for d in bottom}
    for d in selected_all:
        if id(d) in bottom_ids:
            continue
        cv2.drawContours(out, [d.contour.astype(np.int32)], -1, (255, 180, 0), 1)
    for i, d in enumerate(bottom, 1):
        cv2.drawContours(out, [d.contour.astype(np.int32)], -1, (0, 220, 0), 2)
        cc = tuple(np.round(d.center).astype(int))
        cv2.putText(out, f"S{i}", (cc[0] + 3, cc[1] - 3),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.44, (0, 255, 0), 1, cv2.LINE_AA)
    if rows:
        mw = float(np.mean([r["x_width_nm"] for r in rows]))
        mh = float(np.mean([r["y_height_nm"] for r in rows]))
        text = f"bottom-row used={len(rows)} Xmean={mw:.2f}nm Ymean={mh:.2f}nm"
        cv2.rectangle(out, (4, 4), (min(gray.shape[1] - 4, 520), 29), (0, 0, 0), -1)
        cv2.putText(out, text, (8, 22), cv2.FONT_HERSHEY_SIMPLEX,
                    0.50, (255, 255, 255), 1, cv2.LINE_AA)
    return out



# ============================================================
# V22 OVERRIDES
#   1) Via center detection stays V21, but CD is measured from the ACTUAL irregular
#      local boundary rather than an ideal ellipse or two 1-D edge positions.
#   2) Trench/slot Y axis is obtained from the average inclination of the lines that
#      connect corresponding upper/lower pattern centers. X is perpendicular to Y.
#   3) Trench/slot object IDs are left-to-right in the fitted X-axis coordinate.
# ============================================================
V22_VIA_BOUNDARY_ALPHAS = (0.34, 0.40, 0.46, 0.52, 0.58, 0.64, 0.70)
V22_VIA_BOUNDARY_SMOOTH_SIGMA_DIAM_FRAC = 0.018
V22_VIA_BOUNDARY_BG_SIGMA_DIAM_FRAC = 0.42
V22_VIA_BOUNDARY_ROI_DIAM_FACTOR = 1.02
V22_VIA_BOUNDARY_MIN_DIAM_FACTOR = 0.67
V22_VIA_BOUNDARY_MAX_DIAM_FACTOR = 1.58
V22_VIA_BOUNDARY_CENTER_TOL_DIAM_FRAC = 0.25
V22_VIA_BOUNDARY_CLOSE_DIAM_FRAC = 0.025
V22_AXIS_MIN_PAIR_COUNT = 2
V22_AXIS_OUTLIER_DEG = 10.0


def _unit_v22(v: np.ndarray, fallback=(0.0, 1.0)) -> np.ndarray:
    a = np.asarray(v, dtype=float).reshape(2)
    n = float(np.linalg.norm(a))
    if not np.isfinite(n) or n < 1e-9:
        a = np.asarray(fallback, dtype=float)
        n = float(np.linalg.norm(a))
    return a / max(n, 1e-12)


def _axes_from_top_bottom_points_v22(
    point_pairs: List[Tuple[np.ndarray, np.ndarray]],
    image_shape: Tuple[int, int],
) -> Dict[str, object]:
    """Average the inclination angle of top->bottom center lines.

    Image coordinates have +Y downward. Internally y_down follows top->bottom; the
    reported/drawn y_axis_up is its opposite. x_axis is perpendicular and forced to
    point generally to image-right.
    """
    H, W = image_shape[:2]
    angles = []
    mids = []
    lengths = []
    for top, bottom in point_pairs:
        t = np.asarray(top, dtype=float); b = np.asarray(bottom, dtype=float)
        v = b - t
        if v[1] < 0:
            v = -v
        L = float(np.linalg.norm(v))
        if not np.isfinite(L) or L < 3.0:
            continue
        ang = math.degrees(math.atan2(v[1], v[0]))
        # With dy >= 0, angle should normally lie [0,180].
        if ang < 0:
            ang += 180.0
        angles.append(float(ang)); mids.append(0.5 * (t + b)); lengths.append(L)

    fallback = False
    inlier_angles: List[float] = []
    if len(angles) >= V22_AXIS_MIN_PAIR_COUNT:
        aa = np.asarray(angles, dtype=float)
        med = float(np.median(aa))
        dev = np.abs(aa - med)
        mad = 1.4826 * float(np.median(np.abs(dev - np.median(dev))))
        tol = max(V22_AXIS_OUTLIER_DEG, 4.0 * mad)
        good = dev <= tol
        if int(np.count_nonzero(good)) >= V22_AXIS_MIN_PAIR_COUNT:
            inlier_angles = aa[good].tolist()
        else:
            inlier_angles = aa.tolist()
        mean_angle = float(np.mean(inlier_angles))
        y_down = np.array([math.cos(math.radians(mean_angle)), math.sin(math.radians(mean_angle))], dtype=float)
        if y_down[1] < 0:
            y_down = -y_down
        origin = np.median(np.asarray(mids, dtype=float), axis=0) if mids else np.array([0.5*W, 0.5*H])
    else:
        fallback = True
        mean_angle = 90.0
        y_down = np.array([0.0, 1.0], dtype=float)
        origin = np.array([0.5 * (W - 1), 0.5 * (H - 1)], dtype=float)

    y_down = _unit_v22(y_down)
    x_axis = np.array([y_down[1], -y_down[0]], dtype=float)
    if x_axis[0] < 0:
        x_axis = -x_axis
    x_axis = _unit_v22(x_axis, (1.0, 0.0))
    y_up = -y_down
    x_angle = math.degrees(math.atan2(x_axis[1], x_axis[0]))
    y_angle = math.degrees(math.atan2(y_down[1], y_down[0]))
    return {
        'origin': np.asarray(origin, dtype=float),
        'x_axis': x_axis,
        'y_axis_down': y_down,
        'y_axis_up': y_up,
        'x_axis_angle_from_image_x_deg': float(x_angle),
        'y_axis_down_angle_from_image_x_deg': float(y_angle),
        'pair_line_angles_deg': angles,
        'inlier_pair_line_angles_deg': inlier_angles,
        'pair_count': len(angles),
        'fallback': bool(fallback),
    }


def _draw_axes_v22(out: np.ndarray, axis: Dict[str, object], label_prefix: str = '') -> None:
    H, W = out.shape[:2]
    origin = np.asarray(axis.get('origin', [0.5*W, 0.5*H]), dtype=float)
    xa = _unit_v22(np.asarray(axis.get('x_axis', [1,0]), dtype=float), (1,0))
    yu = _unit_v22(np.asarray(axis.get('y_axis_up', [0,-1]), dtype=float), (0,-1))
    L = max(28.0, 0.15 * min(H, W))
    cv2.arrowedLine(out, p2i(origin), p2i(origin + L*xa), (255,255,0), 2, cv2.LINE_AA, tipLength=0.12)
    cv2.arrowedLine(out, p2i(origin), p2i(origin + L*yu), (255,0,255), 2, cv2.LINE_AA, tipLength=0.12)
    cv2.putText(out, f'{label_prefix}X', p2i(origin + 1.05*L*xa), cv2.FONT_HERSHEY_SIMPLEX, 0.45, (255,255,0), 1, cv2.LINE_AA)
    cv2.putText(out, f'{label_prefix}Y', p2i(origin + 1.05*L*yu), cv2.FONT_HERSHEY_SIMPLEX, 0.45, (255,0,255), 1, cv2.LINE_AA)



def _via_radial_gradient_contour_v22(
    gray: np.ndarray,
    center: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    nominal_diam_px: float,
    nearest_pitch: float,
) -> Optional[np.ndarray]:
    """Irregular polygon from the strongest outward dark->bright edge on many rays."""
    H,W=gray.shape; D=max(6.0,float(nominal_diam_px)); R=0.5*D; c=np.asarray(center,float)
    g=cv2.GaussianBlur(gray.astype(np.float32),(0,0),sigmaX=max(0.7,0.018*D),sigmaY=max(0.7,0.018*D))
    rmin=max(2.0,0.58*R); rmax=1.48*R
    if np.isfinite(nearest_pitch) and nearest_pitch>0:
        rmax=min(rmax,0.45*nearest_pitch)
    if rmax<=rmin+2: return None
    nr=max(28,int(math.ceil((rmax-rmin)*3.0))); rr=np.linspace(rmin,rmax,nr,dtype=np.float32)
    radii=[]; angles=np.linspace(0,2*math.pi,96,endpoint=False)
    for a in angles:
        xs=(c[0]+rr*math.cos(a)).astype(np.float32); ys=(c[1]+rr*math.sin(a)).astype(np.float32)
        if xs.min()<1 or ys.min()<1 or xs.max()>W-2 or ys.max()>H-2:
            radii.append(np.nan); continue
        prof=cv2.remap(g,xs.reshape(1,-1),ys.reshape(1,-1),cv2.INTER_LINEAR,borderMode=cv2.BORDER_REFLECT101).ravel().astype(float)
        grad=np.gradient(prof,rr.astype(float))
        gp=np.maximum(grad,0.0); scale=float(np.percentile(np.abs(grad),88))+1e-6
        prior=np.exp(-0.5*((rr-R)/max(1.0,0.42*R))**2)
        score=0.84*np.clip(gp/scale,0,3)+0.16*prior
        j=int(np.argmax(score)); strength=float(score[j])
        radii.append(float(rr[j]) if strength>=0.22 else np.nan)
    rad=np.asarray(radii,float)
    good=np.isfinite(rad)
    if np.count_nonzero(good)<0.55*len(rad): return None
    med=float(np.nanmedian(rad)); rad[~good]=med
    # Circular median-like smoothing removes one-ray spikes but preserves real edge waviness.
    sm=[]; n=len(rad)
    for i in range(n):
        vals=[rad[(i+k)%n] for k in (-2,-1,0,1,2)]
        sm.append(float(np.median(vals)))
    rad=np.asarray(sm,float)
    pts=np.column_stack([c[0]+rad*np.cos(angles),c[1]+rad*np.sin(angles)])
    cont=np.asarray(np.rint(pts),np.int32).reshape(-1,1,2)
    xa=_unit_v22(x_axis,(1,0)); ya=_unit_v22(y_axis,(0,-1))
    wx=_projection_span_v19(cont,xa); hy=_projection_span_v19(cont,ya)
    if not (V22_VIA_BOUNDARY_MIN_DIAM_FACTOR*D<=wx<=V22_VIA_BOUNDARY_MAX_DIAM_FACTOR*D): return None
    if not (V22_VIA_BOUNDARY_MIN_DIAM_FACTOR*D<=hy<=V22_VIA_BOUNDARY_MAX_DIAM_FACTOR*D): return None
    return cont


def _via_actual_boundary_contour_v22(
    gray: np.ndarray,
    center: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    nominal_diam_px: float,
    nearest_pitch: float,
) -> Tuple[Optional[np.ndarray], Dict[str, float]]:
    """Recover the actual, possibly jagged via boundary around a known lattice center.

    V21 is retained for center finding. Here we use local background correction and a
    family of half-level-like thresholds, then choose the closed contour that best
    matches both the expected physical scale and the actual grayscale edge strength.
    No ellipse fitting is used for the normal V22 measurement path.
    """
    H, W = gray.shape
    D = max(6.0, float(nominal_diam_px))
    c = np.asarray(center, dtype=float)
    half = max(6, int(math.ceil(V22_VIA_BOUNDARY_ROI_DIAM_FACTOR * D)))
    if np.isfinite(nearest_pitch) and nearest_pitch > 0:
        half = int(min(half, max(6.0, 0.47 * nearest_pitch)))
    x0 = max(0, int(math.floor(c[0] - half))); x1 = min(W, int(math.ceil(c[0] + half + 1)))
    y0 = max(0, int(math.floor(c[1] - half))); y1 = min(H, int(math.ceil(c[1] + half + 1)))
    if x1-x0 < 9 or y1-y0 < 9:
        return None, {'boundary_score': np.nan}

    roi = gray[y0:y1, x0:x1].astype(np.float32)
    sigma_s = max(0.65, V22_VIA_BOUNDARY_SMOOTH_SIGMA_DIAM_FRAC * D)
    sm = cv2.GaussianBlur(roi, (0,0), sigmaX=sigma_s, sigmaY=sigma_s)
    sigma_bg = max(4.0, V22_VIA_BOUNDARY_BG_SIGMA_DIAM_FRAC * D)
    bg = cv2.GaussianBlur(sm, (0,0), sigmaX=sigma_bg, sigmaY=sigma_bg)
    corr = sm - bg  # dark via -> negative

    yy, xx = np.mgrid[y0:y1, x0:x1]
    rr = np.sqrt((xx-c[0])**2 + (yy-c[1])**2)
    core = rr <= 0.20 * D
    ring_lo = 0.72 * D
    ring_hi = min(0.95 * D, 0.44 * nearest_pitch) if np.isfinite(nearest_pitch) and nearest_pitch > 0 else 0.95 * D
    if ring_hi <= ring_lo + 1:
        ring_lo, ring_hi = 0.64 * D, 0.88 * D
    ring = (rr >= ring_lo) & (rr <= ring_hi)
    if np.count_nonzero(core) < 5 or np.count_nonzero(ring) < 10:
        return None, {'boundary_score': np.nan}
    dark = float(np.median(corr[core])); outer = float(np.median(corr[ring]))
    contrast = outer - dark
    if not np.isfinite(contrast) or contrast <= 0.08:
        # Raw gray fallback in very flat local-background-corrected images.
        dark_raw = float(np.percentile(sm[core], 60.0)); outer_raw = float(np.median(sm[ring]))
        contrast_raw = outer_raw - dark_raw
        if contrast_raw <= 0.15:
            return None, {'boundary_score': np.nan, 'local_contrast': contrast_raw}
        field = sm; dark = dark_raw; outer = outer_raw; contrast = contrast_raw
    else:
        field = corr

    gx = cv2.Sobel(sm, cv2.CV_32F, 1, 0, ksize=3)
    gy = cv2.Sobel(sm, cv2.CV_32F, 0, 1, ksize=3)
    gmag = cv2.magnitude(gx, gy)
    gscale = float(np.percentile(gmag, 92.0)) + 1e-6
    lc = np.array([c[0]-x0, c[1]-y0], dtype=float)
    k = max(1, int(round(V22_VIA_BOUNDARY_CLOSE_DIAM_FRAC * D)))
    if k % 2 == 0: k += 1
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k,k))
    xa = _unit_v22(x_axis, (1,0)); ya = _unit_v22(y_axis, (0,-1))
    best = None
    for alpha in V22_VIA_BOUNDARY_ALPHAS:
        thr = dark + float(alpha) * contrast
        mask = (field <= thr).astype(np.uint8) * 255
        if k > 1:
            mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, kernel, iterations=1)
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        for cc in contours:
            if len(cc) < 8:
                continue
            inside = cv2.pointPolygonTest(cc.astype(np.float32), (float(lc[0]), float(lc[1])), False)
            desc = describe_contour(cc)
            if desc is None:
                continue
            dist = float(np.linalg.norm(desc.center - lc))
            if inside < 0 and dist > V22_VIA_BOUNDARY_CENTER_TOL_DIAM_FRAC * D:
                continue
            cg = cc.astype(np.int32).copy()
            cg[:,0,0] += x0; cg[:,0,1] += y0
            wx = _projection_span_v19(cg, xa); hy = _projection_span_v19(cg, ya)
            if not (V22_VIA_BOUNDARY_MIN_DIAM_FACTOR*D <= wx <= V22_VIA_BOUNDARY_MAX_DIAM_FACTOR*D):
                continue
            if not (V22_VIA_BOUNDARY_MIN_DIAM_FACTOR*D <= hy <= V22_VIA_BOUNDARY_MAX_DIAM_FACTOR*D):
                continue
            pts = cc.reshape(-1,2).astype(int)
            pxs = np.clip(pts[:,0], 0, gmag.shape[1]-1); pys = np.clip(pts[:,1], 0, gmag.shape[0]-1)
            edge_strength = float(np.median(gmag[pys,pxs]) / gscale) if len(pts) else 0.0
            size_geom = math.sqrt(max(wx*hy, 1e-9)) / D
            size_pen = abs(size_geom - 1.0)
            center_pen = dist / D
            # Edge evidence dominates; nominal size is only a guard/prior, not the final shape.
            q = 1.15*edge_strength - 0.55*size_pen - 0.35*center_pen
            rec = (q, cg, wx, hy, edge_strength, alpha, contrast, dist)
            if best is None or rec[0] > best[0]:
                best = rec
    if best is None:
        c0 = _via_radial_gradient_contour_v22(gray, center, xa, ya, D, nearest_pitch)
        if c0 is not None:
            return c0, {'boundary_score': 0.0, 'boundary_alpha': np.nan, 'local_contrast': contrast,
                        'edge_strength': np.nan, 'center_offset_px': np.nan, 'boundary_fallback': 'radial_gradient'}
        c0 = _local_halflevel_contour_v21(gray, center, D)
        if c0 is None:
            return None, {'boundary_score': np.nan, 'local_contrast': contrast}
        wx = _projection_span_v19(c0, xa); hy = _projection_span_v19(c0, ya)
        if not (V22_VIA_BOUNDARY_MIN_DIAM_FACTOR*D <= wx <= V22_VIA_BOUNDARY_MAX_DIAM_FACTOR*D and
                V22_VIA_BOUNDARY_MIN_DIAM_FACTOR*D <= hy <= V22_VIA_BOUNDARY_MAX_DIAM_FACTOR*D):
            return None, {'boundary_score': np.nan, 'local_contrast': contrast}
        return c0, {'boundary_score': 0.0, 'boundary_alpha': 0.52, 'local_contrast': contrast,
                    'edge_strength': np.nan, 'center_offset_px': np.nan, 'boundary_fallback': 'halflevel_contour'}
    return best[1], {'boundary_score': float(best[0]), 'boundary_alpha': float(best[5]),
                     'local_contrast': float(best[6]), 'edge_strength': float(best[4]),
                     'center_offset_px': float(best[7])}


def _measure_selected_via_centers_v21(
    gray: np.ndarray,
    selected: List[ShapeDesc],
    axis: Dict[str, object],
    nominal_diam_px: float,
    nearest_pitch: float,
) -> List[Tuple[ShapeDesc, float, float, Dict[str, float], Dict[str, float], str, str]]:
    """V22 override: actual irregular contour projection, with V21 edge fallback."""
    xa = np.asarray(axis['x_axis'], dtype=float)
    ya = np.asarray(axis['y_axis_up'], dtype=float)
    out = []
    for d in selected:
        contour, cm = _via_actual_boundary_contour_v22(gray, d.center, xa, ya, nominal_diam_px, nearest_pitch)
        if contour is None:
            contour = _via_radial_gradient_contour_v22(gray, d.center, xa, ya, nominal_diam_px, nearest_pitch)
            if contour is not None:
                cm = {'boundary_score': 0.0, 'boundary_alpha': np.nan, 'local_contrast': np.nan,
                      'edge_strength': np.nan, 'center_offset_px': np.nan, 'boundary_fallback': 'radial_gradient_direct'}
        if contour is not None:
            wx = float(_projection_span_v19(contour, xa)); hy = float(_projection_span_v19(contour, ya))
            if (V21_VIA_MEASURE_MIN_DIAM_FACTOR*nominal_diam_px <= wx <= V21_VIA_MEASURE_MAX_DIAM_FACTOR*nominal_diam_px and
                V21_VIA_MEASURE_MIN_DIAM_FACTOR*nominal_diam_px <= hy <= V21_VIA_MEASURE_MAX_DIAM_FACTOR*nominal_diam_px):
                d.contour = contour.astype(np.int32)
                setattr(d, '_v22_actual_boundary', True)
                ad = describe_contour(d.contour)
                if ad is not None:
                    # Keep the lattice center for array geometry, but use actual-boundary morphology metrics.
                    d.area=ad.area; d.perimeter=ad.perimeter; d.circularity=ad.circularity; d.solidity=ad.solidity
                    d.u=ad.u; d.v=ad.v; d.major_min=ad.major_min; d.major_max=ad.major_max
                    d.minor_min=ad.minor_min; d.minor_max=ad.minor_max; d.length_px=ad.length_px; d.width_px=ad.width_px
                    d.angle_from_x_deg=ad.angle_from_x_deg
                mx = {'contrast': cm.get('local_contrast', np.nan), 'boundary_score': cm.get('boundary_score', np.nan),
                      'edge_strength': cm.get('edge_strength', np.nan), 'boundary_alpha': cm.get('boundary_alpha', np.nan)}
                my = dict(mx)
                out.append((d, wx, hy, mx, my, 'actual_irregular_contour', 'actual_irregular_contour'))
                continue
        # Rare fallback: retain V21 physical-radius edge detector rather than losing the via.
        wx, mx = _measure_via_axis_width_v21(gray, d.center, xa, nominal_diam_px, nearest_pitch)
        hy, my = _measure_via_axis_width_v21(gray, d.center, ya, nominal_diam_px, nearest_pitch)
        if not np.isfinite(wx) or not np.isfinite(hy):
            continue
        setattr(d, '_v22_actual_boundary', False)
        out.append((d, float(wx), float(hy), mx, my, 'v21_gradient_fallback', 'v21_gradient_fallback'))
    return out


def _single_trench_body_x_edges_v22(
    gray: np.ndarray, x: int, approx_pitch: float, yrange: Tuple[float,float]
) -> Tuple[float,float,float,Dict[str,float]]:
    H, W = gray.shape
    ya=max(0,int(round(yrange[0]*H))); yb=min(H,int(round(yrange[1]*H)))
    if yb-ya < 3:
        return float(x-2),float(x+2),4.0,{}
    prof=np.median(gray[ya:yb].astype(np.float32),axis=0)
    prof=_smooth_1d_v20(prof,max(0.7,0.0015*W))
    half=int(max(6,min(0.42*approx_pitch if np.isfinite(approx_pitch) else 0.08*W,0.10*W)))
    a=max(1,x-half); b=min(W-2,x+half)
    if b-a<6:
        return float(x-2),float(x+2),4.0,{}
    local=prof[a:b+1]; ci=int(np.clip(x-a,1,len(local)-2))
    outer=np.r_[local[:max(2,len(local)//5)],local[-max(2,len(local)//5):]]
    core_r=max(1,int(round(0.10*len(local))))
    core=local[max(0,ci-core_r):min(len(local),ci+core_r+1)]
    dark=float(np.percentile(core,30)); bg=float(np.median(outer)); contrast=bg-dark
    level=dark+0.50*max(contrast,1.0); grad=np.gradient(local)
    left_idx=np.arange(1,ci); right_idx=np.arange(ci+1,len(local)-1)
    if len(left_idx)==0 or len(right_idx)==0:
        return float(x-2),float(x+2),4.0,{'contrast':contrast}
    # First threshold crossing outward from center, then sharpen around it with gradient.
    lvalid=left_idx[local[left_idx] <= level]
    lh=int(lvalid[0]) if len(lvalid) else int(left_idx[np.argmin(grad[left_idx])])
    rvalid=right_idx[local[right_idx] >= level]
    rh=int(rvalid[0]) if len(rvalid) else int(right_idx[np.argmax(grad[right_idx])])
    sr=max(2,int(round(0.12*(b-a))))
    ls=np.arange(max(1,lh-sr),min(ci-1,lh+sr)+1); rs=np.arange(max(ci+1,rh-sr),min(len(local)-2,rh+sr)+1)
    lg=int(ls[np.argmin(grad[ls])]) if len(ls) else lh; rg=int(rs[np.argmax(grad[rs])]) if len(rs) else rh
    left=a+0.62*lg+0.38*lh; right=a+0.62*rg+0.38*rh
    if right<=left+1: left,right=float(a+lh),float(a+rh)
    return float(left),float(right),float(max(1,right-left)),{'contrast':contrast}


def _select_and_measure_trenches_v19(
    gray: np.ndarray, pool: List[ShapeDesc], px_nm: float, image_name: str,
    stage: str, condition: str, pattern_key: str,
) -> Tuple[List[dict], List[dict], List[dict], Dict[str, object]]:
    pairs, meta = _detect_trench_pairs_profile_v20(gray)
    H,W=gray.shape
    if pairs:
        xs=sorted([p['common_x'] for p in pairs]); diffs=np.diff(xs)
        approx_pitch=float(np.median(diffs)) if len(diffs) else float(W/max(4,V20_TRENCH_EXPECTED_GROUPS_HINT))
    else:
        approx_pitch=np.nan
    point_pairs=[]
    for p in pairs:
        txl,txr,tw,_=_single_trench_body_x_edges_v22(gray,int(round(p['common_x'])),approx_pitch,V20_TRENCH_TOP_BODY_RANGE)
        bxl,bxr,bw,_=_single_trench_body_x_edges_v22(gray,int(round(p['common_x'])),approx_pitch,V20_TRENCH_BOTTOM_BODY_RANGE)
        tx=0.5*(txl+txr); bx=0.5*(bxl+bxr)
        # Visible-segment centers are used for the requested center-to-center inclination.
        top_cy=0.5*(V20_TRENCH_TOP_BODY_RANGE[0]*H + p['top_tip_y_px'])
        bot_cy=0.5*(p['bottom_tip_y_px'] + V20_TRENCH_BOTTOM_BODY_RANGE[1]*H)
        tp=np.array([tx,top_cy],dtype=float); bp=np.array([bx,bot_cy],dtype=float)
        p.update({'top_x_left_px':txl,'top_x_right_px':txr,'bottom_x_left_px':bxl,'bottom_x_right_px':bxr,
                  'top_center_x_px':tx,'bottom_center_x_px':bx,'top_center_y_px':top_cy,'bottom_center_y_px':bot_cy,
                  'top_width_px_v22':tw,'bottom_width_px_v22':bw})
        point_pairs.append((tp,bp))
    axis=_axes_from_top_bottom_points_v22(point_pairs,gray.shape)
    xa=np.asarray(axis['x_axis'],float); yd=np.asarray(axis['y_axis_down'],float); yu=np.asarray(axis['y_axis_up'],float); origin=np.asarray(axis['origin'],float)
    for p in pairs:
        top_tip=np.array([p.get('top_center_x_px',p['common_x']),p['top_tip_y_px']],float)
        bot_tip=np.array([p.get('bottom_center_x_px',p['common_x']),p['bottom_tip_y_px']],float)
        mid=0.5*(top_tip+bot_tip)
        p['tip_gap_axis_px']=float((bot_tip-top_tip)@yd)
        p['axis_center_x_px']=float((mid-origin)@xa)
        p['axis_center_y_up_px']=float((mid-origin)@yu)
        # Width transverse to the fitted Y direction. For these narrow, near-vertical bars,
        # the independently measured top/bottom side-edge spans are projected onto X.
        proj_factor=max(0.2,abs(float(xa[0])))
        p['x_width_axis_px']=0.5*(p.get('top_width_px_v22',p['width_px'])+p.get('bottom_width_px_v22',p['width_px']))*proj_factor
    pairs=sorted(pairs,key=lambda p:p.get('axis_center_x_px',p['common_x']))
    rows=[]
    for oid,p in enumerate(pairs,1):
        rows.append({
            'stage':stage,'stage_zh':STAGE_LABELS_ZH_V19[stage],'condition':condition,
            'condition_label':CONDITION_LABELS_V19.get(str(condition),f'Condition {condition}'),
            'image':image_name,'pattern':pattern_key,'pattern_zh':PATTERN_LABELS_ZH_V19[pattern_key],
            'object_id':oid,'used_for_statistics':True,'is_axis_bottom_row':False,
            'x_width_nm':p['x_width_axis_px']*px_nm,'y_height_nm':np.nan,'tip_gap_y_nm':p['tip_gap_axis_px']*px_nm,
            'x_width_px':p['x_width_axis_px'],'y_height_px':np.nan,
            'center_image_x_px':0.5*(p.get('top_center_x_px',p['common_x'])+p.get('bottom_center_x_px',p['common_x'])),
            'center_image_y_px':0.5*(p['top_tip_y_px']+p['bottom_tip_y_px']),
            'center_axis_x_px':p['axis_center_x_px'],'center_axis_y_up_px':p['axis_center_y_up_px'],
            'top_trench_width_nm':p.get('top_width_px_v22',p['width_px'])*px_nm,
            'bottom_trench_width_nm':p.get('bottom_width_px_v22',p['width_px'])*px_nm,
            'pair_x_offset_nm':abs(p.get('bottom_center_x_px',p['common_x'])-p.get('top_center_x_px',p['common_x']))*px_nm,
            'pair_x_overlap_frac':1.0,'width_mismatch_frac':abs(p.get('top_width_px_v22',p['width_px'])-p.get('bottom_width_px_v22',p['width_px']))/max(1e-9,0.5*(p.get('top_width_px_v22',p['width_px'])+p.get('bottom_width_px_v22',p['width_px']))),
            'pair_score':p['score'],'pixel_size_nm':px_nm,'trench_detector':'V20_gray_profile+V22_rotated_axis',
            'top_tip_contrast_gray':p.get('top_tip_contrast',np.nan),'bottom_tip_contrast_gray':p.get('bottom_tip_contrast',np.nan),
            'top_tip_gradient':p.get('top_tip_gradient',np.nan),'bottom_tip_gradient':p.get('bottom_tip_gradient',np.nan),
            'x_axis_angle_from_image_x_deg':axis['x_axis_angle_from_image_x_deg'],
            'y_axis_down_angle_from_image_x_deg':axis['y_axis_down_angle_from_image_x_deg'],
            'axis_pair_count':axis['pair_count'],'axis_fallback':axis['fallback'],
        })
    meta=dict(meta); meta.update({
        'v22_axis_pair_count':axis['pair_count'],'v22_axis_fallback':axis['fallback'],
        'v22_x_axis_angle_from_image_x_deg':axis['x_axis_angle_from_image_x_deg'],
        'v22_y_axis_down_angle_from_image_x_deg':axis['y_axis_down_angle_from_image_x_deg'],
    })
    for p in pairs: p['_v22_axis']=axis
    return rows,pairs,[],meta


def _annotate_trenches_v19(gray: np.ndarray, pairs: List[dict], rows: List[dict]) -> np.ndarray:
    out=to_color(gray); H,W=gray.shape; body_margin=int(round(V20_TRENCH_DRAW_BODY_LIMIT_FRAC*H))
    for i,(p,r) in enumerate(zip(pairs,rows),1):
        txl=int(round(p.get('top_x_left_px',p['x_left_px']))); txr=int(round(p.get('top_x_right_px',p['x_right_px'])))
        bxl=int(round(p.get('bottom_x_left_px',p['x_left_px']))); bxr=int(round(p.get('bottom_x_right_px',p['x_right_px'])))
        ty=int(round(p['top_tip_y_px'])); by=int(round(p['bottom_tip_y_px']))
        cv2.rectangle(out,(max(0,txl),max(0,body_margin)),(min(W-1,txr),max(body_margin,ty)),(0,210,0),1)
        cv2.rectangle(out,(max(0,bxl),min(H-1,by)),(min(W-1,bxr),min(H-1,H-1-body_margin)),(0,210,0),1)
        t=(int(round(p.get('top_center_x_px',p['common_x']))),ty); b=(int(round(p.get('bottom_center_x_px',p['common_x']))),by)
        cv2.arrowedLine(out,t,b,(255,255,255),1,cv2.LINE_AA,tipLength=0.12); cv2.arrowedLine(out,b,t,(255,255,255),1,cv2.LINE_AA,tipLength=0.12)
        cv2.putText(out,f'T{i}',(min(W-35,int(round(0.5*(t[0]+b[0])))+3),int(round(0.5*(ty+by)))),cv2.FONT_HERSHEY_SIMPLEX,0.40,(255,255,255),1,cv2.LINE_AA)
    if pairs and '_v22_axis' in pairs[0]: _draw_axes_v22(out,pairs[0]['_v22_axis'])
    if rows:
        text=f"trench n={len(rows)} Xmean={np.mean([r['x_width_nm'] for r in rows]):.2f}nm  Ygap mean={np.mean([r['tip_gap_y_nm'] for r in rows]):.2f}nm"
    else: text='trench: no valid opposing pairs'
    cv2.rectangle(out,(4,4),(min(W-4,650),29),(0,0,0),-1); cv2.putText(out,text,(8,22),cv2.FONT_HERSHEY_SIMPLEX,0.46,(255,255,255),1,cv2.LINE_AA)
    return out


def _slot_axis_pairs_v22(paired: List[Tuple[ShapeDesc,dict]], image_shape: Tuple[int,int]) -> Tuple[Dict[str,object], List[Tuple[np.ndarray,np.ndarray]]]:
    if not paired:
        return _axes_from_top_bottom_points_v22([],image_shape),[]
    rows=sorted(set(int(r.get('slot_row',0)) for _,r in paired if int(r.get('slot_row',0))>0))
    if len(rows)>=2:
        med={rid:float(np.median([d.center[1] for d,r in paired if int(r.get('slot_row',0))==rid])) for rid in rows}
        top_id=min(med,key=med.get); bot_id=max(med,key=med.get)
        top=sorted([d for d,r in paired if int(r.get('slot_row',0))==top_id],key=lambda d:d.center[0])
        bot=sorted([d for d,r in paired if int(r.get('slot_row',0))==bot_id],key=lambda d:d.center[0])
    else:
        all_d=[d for d,_ in paired]; ys=np.asarray([d.center[1] for d in all_d],float); medy=float(np.median(ys))
        top=sorted([d for d in all_d if d.center[1]<=medy],key=lambda d:d.center[0]); bot=sorted([d for d in all_d if d.center[1]>medy],key=lambda d:d.center[0])
    if not top or not bot:
        return _axes_from_top_bottom_points_v22([],image_shape),[]
    # Pair by left-to-right rank when counts match; otherwise greedily use nearest image-X.
    pp=[]
    if len(top)==len(bot):
        pp=[(np.asarray(a.center,float),np.asarray(b.center,float)) for a,b in zip(top,bot)]
    else:
        unused=set(range(len(bot)))
        for a in top:
            if not unused: break
            j=min(unused,key=lambda k:abs(float(bot[k].center[0]-a.center[0])))
            pp.append((np.asarray(a.center,float),np.asarray(bot[j].center,float))); unused.remove(j)
    return _axes_from_top_bottom_points_v22(pp,image_shape),pp


def _select_and_measure_bottom_slots_v19(
    gray: np.ndarray, pool: List[ShapeDesc], px_nm: float, image_name: str,
    stage: str, condition: str, pattern_key: str,
) -> Tuple[List[dict], List[ShapeDesc], List[ShapeDesc], List[dict], Dict[str, object]]:
    old_rows,selected_raw,_,rejects,meta=measure_slot_v18(gray,pool,px_nm,image_name,stage)
    paired=[]
    for d,row in zip(selected_raw,old_rows):
        if _bbox_complete_v19(d,gray.shape,0.10): paired.append((d,row))
        else: rejects.append(_make_reject_from_contour(d.contour,'slot_v22_border_qc','slot_near_or_crossing_image_edge','V18-selected slot may be incomplete at image edge',desc=d))
    selected_all=[d for d,_ in paired]
    if paired:
        row_ids=sorted(set(int(row.get('slot_row',0)) for _,row in paired if int(row.get('slot_row',0))>0))
        if row_ids:
            med_y={rid:float(np.median([d.center[1] for d,row in paired if int(row.get('slot_row',0))==rid])) for rid in row_ids}
            bottom_row_id=max(med_y,key=med_y.get); bottom=[d for d,row in paired if int(row.get('slot_row',0))==bottom_row_id]
        else:
            bottom_row_id=np.nan; bottom=_generic_bottom_slot_row_v19(selected_all)
    else:
        strict,initial_rejects,_,_,_,_=_slot_initial_filter_v12(pool); rejects.extend(initial_rejects)
        strict=[d for d in _dedup_v19(strict) if _bbox_complete_v19(d,gray.shape,0.10)]
        selected_all=strict; bottom=_generic_bottom_slot_row_v19(strict); bottom_row_id=np.nan; paired=[]
    bottom,rr=_area_consistency_filter_v19(bottom,'slot_v22_bottom_row_area_qc',lo_factor=0.40,hi_factor=2.50); rejects.extend(rr)
    # Fit axes from top/bottom corresponding slot centers using ALL valid V18 rows.
    axis,axis_pairs=_slot_axis_pairs_v22(paired,gray.shape)
    xa=np.asarray(axis['x_axis'],float); yu=np.asarray(axis['y_axis_up'],float); origin=np.asarray(axis['origin'],float)
    bottom=sorted(bottom,key=lambda d:float((d.center-origin)@xa)); bottom_ids={id(d) for d in bottom}
    rows=[]
    for oid,d in enumerate(bottom,1):
        w_px=_projection_span_v19(d.contour,xa); h_px=_projection_span_v19(d.contour,yu); reg=slot_regularity_metrics(d)
        rows.append({
            'stage':stage,'stage_zh':STAGE_LABELS_ZH_V19[stage],'condition':condition,'condition_label':CONDITION_LABELS_V19.get(str(condition),f'Condition {condition}'),
            'image':image_name,'pattern':pattern_key,'pattern_zh':PATTERN_LABELS_ZH_V19[pattern_key],'object_id':oid,'used_for_statistics':True,'is_axis_bottom_row':True,
            'x_width_nm':w_px*px_nm,'y_height_nm':h_px*px_nm,'tip_gap_y_nm':np.nan,'x_width_px':w_px,'y_height_px':h_px,
            'center_image_x_px':float(d.center[0]),'center_image_y_px':float(d.center[1]),'center_axis_x_px':float((d.center-origin)@xa),'center_axis_y_up_px':float((d.center-origin)@yu),
            'slot_row_id_from_v18':bottom_row_id,'slot_pixel_area_px':contour_pixel_count(d.contour),'slot_regularity_solidity':reg['regularity_solidity'],
            'slot_regularity_convexity':reg['regularity_convexity'],'slot_best_template_iou':reg['regularity_template_iou'],'selection_source':_source_of(d),'pixel_size_nm':px_nm,
            'x_axis_angle_from_image_x_deg':axis['x_axis_angle_from_image_x_deg'],'y_axis_down_angle_from_image_x_deg':axis['y_axis_down_angle_from_image_x_deg'],
            'axis_pair_count':axis['pair_count'],'axis_fallback':axis['fallback'],
        })
    for d,old in paired:
        if int(old.get('slot_row',-999))==bottom_row_id and id(d) not in bottom_ids and not any(r.get('_desc') is d for r in rejects):
            rejects.append(_make_reject_from_contour(d.contour,'slot_v22_bottom_row_final_qc','bottom_row_slot_not_used','bottom-row slot removed by final border/area consistency QC',desc=d))
    meta=dict(meta); meta.update({'slot_bottom_row_id':bottom_row_id,'slot_bottom_row_used_count':len(bottom),'slot_all_valid_count':len(selected_all),
                                 'v22_axis_pair_count':axis['pair_count'],'v22_axis_fallback':axis['fallback'],
                                 'v22_x_axis_angle_from_image_x_deg':axis['x_axis_angle_from_image_x_deg'],'v22_y_axis_down_angle_from_image_x_deg':axis['y_axis_down_angle_from_image_x_deg']})
    for d in selected_all: setattr(d,'_v22_slot_axis',axis)
    return rows,selected_all,bottom,rejects,meta


def _annotate_slots_v19(gray: np.ndarray, selected_all: List[ShapeDesc], bottom: List[ShapeDesc], rows: List[dict]) -> np.ndarray:
    out=to_color(gray); bottom_ids={id(d) for d in bottom}
    for d in selected_all:
        if id(d) not in bottom_ids: cv2.drawContours(out,[d.contour.astype(np.int32)],-1,(255,180,0),1)
    for i,d in enumerate(bottom,1):
        cv2.drawContours(out,[d.contour.astype(np.int32)],-1,(0,220,0),2); cc=tuple(np.round(d.center).astype(int))
        cv2.putText(out,f'S{i}',(cc[0]+3,cc[1]-3),cv2.FONT_HERSHEY_SIMPLEX,0.44,(0,255,0),1,cv2.LINE_AA)
    if selected_all and hasattr(selected_all[0],'_v22_slot_axis'): _draw_axes_v22(out,getattr(selected_all[0],'_v22_slot_axis'))
    if rows:
        mw=float(np.mean([r['x_width_nm'] for r in rows])); mh=float(np.mean([r['y_height_nm'] for r in rows])); text=f'bottom-row n={len(rows)} Xmean={mw:.2f}nm Ymean={mh:.2f}nm'
    else: text='slot: no valid bottom-row object'
    cv2.rectangle(out,(4,4),(min(gray.shape[1]-4,560),29),(0,0,0),-1); cv2.putText(out,text,(8,22),cv2.FONT_HERSHEY_SIMPLEX,0.48,(255,255,255),1,cv2.LINE_AA)
    return out

# -----------------------------
# V19 SUMMARIES AND PLOTS
# -----------------------------
METRICS_BY_PATTERN_V19: Dict[str, List[Tuple[str, str]]] = {
    "via40": [("x_width_nm", "X轴宽度"), ("y_height_nm", "Y轴高度")],
    "via60": [("x_width_nm", "X轴宽度"), ("y_height_nm", "Y轴高度")],
    "trench160": [("x_width_nm", "上下矩形平均宽度"), ("tip_gap_y_nm", "tip间Y轴距离")],
    "slot210": [("x_width_nm", "最下排X轴宽度"), ("y_height_nm", "最下排Y轴高度")],
}


def _metric_summary_v19(
    objects_df: pd.DataFrame,
    group_cols: List[str],
) -> pd.DataFrame:
    records: List[dict] = []
    if objects_df.empty:
        return pd.DataFrame()
    for keys, g in objects_df.groupby(group_cols, dropna=False, sort=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        base = dict(zip(group_cols, keys))
        pattern = str(g["pattern"].iloc[0])
        for metric, metric_zh in METRICS_BY_PATTERN_V19.get(pattern, []):
            vals = pd.to_numeric(g[metric], errors="coerce").dropna().to_numpy(dtype=float)
            if len(vals) == 0:
                continue
            rec = dict(base)
            rec.update({
                "pattern": pattern,
                "pattern_zh": PATTERN_LABELS_ZH_V19.get(pattern, pattern),
                "metric": metric,
                "metric_zh": metric_zh,
                "n": int(len(vals)),
                "mean_nm": float(np.mean(vals)),
                "std_nm": float(np.std(vals, ddof=1)) if len(vals) >= 2 else np.nan,
                "sem_nm": float(np.std(vals, ddof=1) / math.sqrt(len(vals))) if len(vals) >= 2 else np.nan,
                "median_nm": float(np.median(vals)),
                "min_nm": float(np.min(vals)),
                "max_nm": float(np.max(vals)),
            })
            records.append(rec)
    return pd.DataFrame(records)


def _before_after_wide_v19(condition_metric_df: pd.DataFrame) -> pd.DataFrame:
    if condition_metric_df.empty:
        return pd.DataFrame()
    index_cols = ["condition", "condition_label", "pattern", "pattern_zh", "metric", "metric_zh"]
    out = condition_metric_df[index_cols].drop_duplicates().copy()
    for stage in STAGE_ORDER_V19:
        s = condition_metric_df[condition_metric_df["stage"] == stage]
        cols = index_cols + ["n", "mean_nm", "std_nm", "sem_nm", "median_nm"]
        s = s[cols].rename(columns={
            "n": f"{stage}_n",
            "mean_nm": f"{stage}_mean_nm",
            "std_nm": f"{stage}_std_nm",
            "sem_nm": f"{stage}_sem_nm",
            "median_nm": f"{stage}_median_nm",
        })
        out = out.merge(s, on=index_cols, how="left")
    for col in ("before_mean_nm", "after_mean_nm"):
        if col not in out.columns:
            out[col] = np.nan
    out["delta_after_minus_before_nm"] = out["after_mean_nm"] - out["before_mean_nm"]
    before = pd.to_numeric(out["before_mean_nm"], errors="coerce")
    out["change_percent"] = np.where(
        np.isfinite(before) & (np.abs(before) > 1e-12),
        100.0 * out["delta_after_minus_before_nm"] / before,
        np.nan,
    )
    condition_order = {s: i for i, s in enumerate(CONDITION_IDS_V19)}
    pattern_order = {p: i for i, p in enumerate(PATTERN_KEYS_V19)}
    out["_condition_order"] = out["condition"].astype(str).map(condition_order).fillna(999)
    out["_pattern_order"] = out["pattern"].map(pattern_order).fillna(999)
    out = out.sort_values(["_condition_order", "_pattern_order", "metric"]).drop(columns=["_condition_order", "_pattern_order"])
    return out.reset_index(drop=True)


def _save_grouped_metric_bar_v19(
    condition_metric_df: pd.DataFrame,
    pattern: str,
    metric: str,
    metric_zh: str,
    out_path: Path,
) -> None:
    sub = condition_metric_df[
        (condition_metric_df["pattern"] == pattern)
        & (condition_metric_df["metric"] == metric)
    ].copy()
    if sub.empty:
        return
    _configure_chinese_matplotlib()
    fig, ax = plt.subplots(figsize=(12.5, 6.4))
    x = np.arange(len(CONDITION_IDS_V19), dtype=float)
    width = 0.36
    for j, stage in enumerate(STAGE_ORDER_V19):
        means = []
        errs = []
        ns = []
        for condition in CONDITION_IDS_V19:
            q = sub[(sub["stage"] == stage) & (sub["condition"].astype(str) == str(condition))]
            if q.empty:
                means.append(np.nan); errs.append(0.0); ns.append(0)
            else:
                row = q.iloc[0]
                means.append(float(row["mean_nm"]))
                if PLOT_ERROR_BAR_V19 == "std":
                    e = row["std_nm"]
                elif PLOT_ERROR_BAR_V19 == "sem":
                    e = row["sem_nm"]
                else:
                    e = 0.0
                errs.append(float(e) if np.isfinite(e) else 0.0)
                ns.append(int(row["n"]))
        offset = (j - 0.5) * width
        bars = ax.bar(
            x + offset,
            means,
            width,
            yerr=errs if PLOT_ERROR_BAR_V19 != "none" else None,
            capsize=3,
            label=STAGE_LABELS_ZH_V19[stage],
        )
        for bar, value, nobj in zip(bars, means, ns):
            if not np.isfinite(value):
                continue
            ax.annotate(
                f"{value:.2f}\nn={nobj}",
                xy=(bar.get_x() + bar.get_width() / 2.0, bar.get_height()),
                xytext=(0, 4),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=8,
            )
    err_name = {"std": "标准差", "sem": "均值标准误", "none": "无误差棒"}[PLOT_ERROR_BAR_V19]
    ax.set_xticks(x)
    ax.set_xticklabels([f"Condition {s}" for s in CONDITION_IDS_V19])
    ax.set_ylabel("尺寸 / nm")
    ax.set_title(f"{PATTERN_LABELS_ZH_V19[pattern]}：{metric_zh}处理前后对比（误差棒：{err_name}）")
    ax.legend()
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=PLOT_DPI_V19, bbox_inches="tight")
    plt.close(fig)


def _make_all_plots_v19(condition_metric_df: pd.DataFrame, plot_dir: Path) -> None:
    plot_dir.mkdir(parents=True, exist_ok=True)
    for pattern, metrics in METRICS_BY_PATTERN_V19.items():
        for metric, metric_zh in metrics:
            _save_grouped_metric_bar_v19(
                condition_metric_df,
                pattern,
                metric,
                metric_zh,
                plot_dir / f"{pattern}__{metric}__before_after.png",
            )



V22_PER_OBJECT_METRICS: Dict[str, List[Tuple[str, str]]] = {
    'trench160': [('x_width_nm', '上下矩形平均宽度'), ('tip_gap_y_nm', 'tip间Y轴距离')],
    'slot210': [('x_width_nm', 'X轴宽度'), ('y_height_nm', 'Y轴高度')],
}


def _object_before_after_wide_v22(objects_df: pd.DataFrame) -> pd.DataFrame:
    """One row per Condition/pattern/left-to-right object number/metric."""
    if objects_df.empty:
        return pd.DataFrame()
    records=[]
    for pattern, metrics in V22_PER_OBJECT_METRICS.items():
        gpat=objects_df[objects_df['pattern']==pattern].copy()
        if gpat.empty: continue
        for condition in CONDITION_IDS_V19:
            gc=gpat[gpat['condition'].astype(str)==str(condition)]
            if gc.empty: continue
            ids=sorted(set(pd.to_numeric(gc['object_id'],errors='coerce').dropna().astype(int).tolist()))
            for oid in ids:
                go=gc[pd.to_numeric(gc['object_id'],errors='coerce')==oid]
                for metric,metric_zh in metrics:
                    rec={'condition':str(condition),'condition_label':CONDITION_LABELS_V19.get(str(condition),f'Condition {condition}'),
                         'pattern':pattern,'pattern_zh':PATTERN_LABELS_ZH_V19.get(pattern,pattern),'object_id':int(oid),
                         'object_label':('T' if pattern=='trench160' else 'S')+str(int(oid)),
                         'metric':metric,'metric_zh':metric_zh}
                    for stage in STAGE_ORDER_V19:
                        gs=go[go['stage']==stage]
                        vals=pd.to_numeric(gs[metric],errors='coerce').dropna() if metric in gs.columns else pd.Series(dtype=float)
                        rec[f'{stage}_value_nm']=float(vals.iloc[0]) if len(vals) else np.nan
                        rec[f'{stage}_image']=str(gs['image'].iloc[0]) if len(gs) and 'image' in gs.columns else ''
                    bv=rec.get('before_value_nm',np.nan); av=rec.get('after_value_nm',np.nan)
                    rec['delta_after_minus_before_nm']=float(av-bv) if np.isfinite(av) and np.isfinite(bv) else np.nan
                    rec['change_percent']=float(100.0*(av-bv)/bv) if np.isfinite(av) and np.isfinite(bv) and abs(bv)>1e-12 else np.nan
                    records.append(rec)
    if not records: return pd.DataFrame()
    out=pd.DataFrame(records)
    order={s:i for i,s in enumerate(CONDITION_IDS_V19)}; po={'trench160':0,'slot210':1}
    out['_co']=out['condition'].map(order).fillna(999); out['_po']=out['pattern'].map(po).fillna(999)
    out=out.sort_values(['_co','_po','metric','object_id']).drop(columns=['_co','_po']).reset_index(drop=True)
    return out


def _save_per_object_bar_v22(
    objects_df: pd.DataFrame, condition: str, pattern: str, metric: str, metric_zh: str, out_path: Path
) -> None:
    sub=objects_df[(objects_df['condition'].astype(str)==str(condition)) & (objects_df['pattern']==pattern)].copy()
    if sub.empty or metric not in sub.columns: return
    ids=sorted(set(pd.to_numeric(sub['object_id'],errors='coerce').dropna().astype(int).tolist()))
    if not ids: return
    _configure_chinese_matplotlib(); fig_w=max(9.0,1.05*len(ids)+4.0); fig,ax=plt.subplots(figsize=(fig_w,6.0))
    x=np.arange(len(ids),dtype=float); width=0.36
    prefix='T' if pattern=='trench160' else 'S'
    all_values=[]
    for j,stage in enumerate(STAGE_ORDER_V19):
        vals=[]
        for oid in ids:
            q=sub[(sub['stage']==stage) & (pd.to_numeric(sub['object_id'],errors='coerce')==oid)]
            v=pd.to_numeric(q[metric],errors='coerce').dropna()
            vals.append(float(v.iloc[0]) if len(v) else np.nan)
        all_values.extend([v for v in vals if np.isfinite(v)])
        xpos=x+(j-0.5)*width
        # Draw each finite bar separately: if its before/after counterpart is absent,
        # the existing bar remains visible by itself as requested.
        first=True
        for xx,v in zip(xpos,vals):
            if not np.isfinite(v): continue
            bar=ax.bar([xx],[v],width,label=STAGE_LABELS_ZH_V19[stage] if first else None)[0]
            first=False
            ax.annotate(f'{v:.2f}',xy=(bar.get_x()+bar.get_width()/2,bar.get_height()),xytext=(0,4),textcoords='offset points',ha='center',va='bottom',fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels([f'{prefix}{i}' for i in ids])
    ax.set_xlabel('从左到右图案序号'); ax.set_ylabel('尺寸 / nm')
    ax.set_title(f"Condition {condition}  {PATTERN_LABELS_ZH_V19[pattern]}：逐图案{metric_zh}处理前后对比")
    handles,labels=ax.get_legend_handles_labels()
    if handles: ax.legend()
    ax.grid(axis='y',alpha=0.25); fig.tight_layout(); out_path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out_path,dpi=PLOT_DPI_V19,bbox_inches='tight'); plt.close(fig)


def _make_per_object_plots_v22(objects_df: pd.DataFrame, plot_dir: Path) -> None:
    if objects_df.empty: return
    root=plot_dir/'per_object'
    for condition in CONDITION_IDS_V19:
        for pattern,metrics in V22_PER_OBJECT_METRICS.items():
            for metric,metric_zh in metrics:
                _save_per_object_bar_v22(objects_df,condition,pattern,metric,metric_zh,
                    root/f'Condition_{condition}'/f'{pattern}__{metric}__objects_before_after.png')


def _autofit_excel_v19(writer: pd.ExcelWriter) -> None:
    for ws in writer.book.worksheets:
        ws.freeze_panes = "A2"
        if ws.max_row >= 1 and ws.max_column >= 1:
            ws.auto_filter.ref = ws.dimensions
        for col_cells in ws.columns:
            values = [str(c.value) if c.value is not None else "" for c in col_cells[:300]]
            max_len = max([len(v) for v in values] + [8])
            ws.column_dimensions[col_cells[0].column_letter].width = min(46, max_len + 2)


# -----------------------------
# V19 IMAGE/FOLDER PIPELINE
# -----------------------------
def _empty_array_meta_v19() -> Dict[str, object]:
    return {
        "array_assist_used": False,
        "array_expected_positions": 0,
        "array_recovered_count": 0,
        "array_lattice_type": "none",
    }


def _process_one_image_v19(
    image_path: Path,
    stage: str,
    condition: str,
    ann_dir: Path,
) -> Tuple[List[dict], dict, List[dict]]:
    pattern_key = _classify_pattern_key_v19(image_path)
    if pattern_key is None:
        return [], {
            "stage": stage,
            "stage_zh": STAGE_LABELS_ZH_V19[stage],
            "condition": condition,
            "condition_label": CONDITION_LABELS_V19.get(str(condition), f"Condition {condition}"),
            "image": image_path.name,
            "path": str(image_path),
            "pattern": "unknown",
            "status": "SKIPPED_UNKNOWN_PREFIX",
            "note": "Filename must begin with 40, 60, 160, or 210",
        }, []

    base_pattern = _base_pattern_v19(pattern_key)
    print(f"    {image_path.name} -> {pattern_key}", flush=True)
    px_info = read_pixel_size_nm(image_path)
    px_nm = px_info.value_nm
    gray = _read_gray_image_unicode_v19(image_path)
    mask = segment_features(gray)
    primary_descs, segmentation_rejects = extract_descriptors_with_rejections(mask)
    search_meta: Dict[str, object] = {}
    measure_meta: Dict[str, object] = _empty_array_meta_v19()

    if base_pattern in {"via", "slot"}:
        candidate_pool, search_meta = build_candidate_pool_v12(
            gray, primary_descs, base_pattern, image_path.name
        )
    elif base_pattern == "trench":
        # V20 trench detector works directly from grayscale profiles because the upper
        # and lower dark bars are intentionally clipped by the image top/bottom.
        # Old contour pools/completeness filters are therefore bypassed.
        candidate_pool = []
        search_meta = {"trench_search_mode": "V20_gray_profile_direct"}
        segmentation_rejects = []
    else:
        candidate_pool = list(primary_descs)

    rows: List[dict] = []
    shape_rejects: List[dict] = []
    if base_pattern == "via":
        rows, selected, shape_rejects, bundle = _select_and_measure_vias_v19(
            gray, candidate_pool, px_nm, image_path.name, stage, condition, pattern_key
        )
        measure_meta = bundle["meta"]
        annotated = _annotate_vias_v19(gray, selected, rows, bundle)
        n_selected = len(selected)
    elif base_pattern == "trench":
        rows, pairs, shape_rejects, measure_meta = _select_and_measure_trenches_v19(
            gray, candidate_pool, px_nm, image_path.name, stage, condition, pattern_key
        )
        annotated = _annotate_trenches_v19(gray, pairs, rows)
        n_selected = 2 * len(pairs)
    elif base_pattern == "slot":
        rows, selected_all, bottom, shape_rejects, measure_meta = _select_and_measure_bottom_slots_v19(
            gray, candidate_pool, px_nm, image_path.name, stage, condition, pattern_key
        )
        annotated = _annotate_slots_v19(gray, selected_all, bottom, rows)
        n_selected = len(bottom)
    else:
        annotated = to_color(gray)
        n_selected = 0

    all_rejects = segmentation_rejects + shape_rejects
    exported_rejects: List[dict] = []
    for rid, r in enumerate(all_rejects, 1):
        r["reject_id"] = rid
        r["stage"] = stage
        r["condition"] = condition
        r["image"] = image_path.name
        r["path"] = str(image_path)
        r["pattern"] = pattern_key
        r["pixel_size_nm"] = px_nm
        exported_rejects.append(exportable_rejection_record(r))
    if DRAW_REJECTED_OBJECTS_V19:
        annotated = annotate_rejections(annotated, all_rejects)

    safe = _safe_stem_v19(image_path)
    ann_path = ann_dir / f"{safe}__{pattern_key}__annotated.png"
    _imwrite_unicode_v19(ann_path, annotated)
    if SAVE_BINARY_MASK_V19:
        _imwrite_unicode_v19(ann_dir / f"{safe}__primary_mask.png", mask)

    status = {
        "stage": stage,
        "stage_zh": STAGE_LABELS_ZH_V19[stage],
        "condition": condition,
        "condition_label": CONDITION_LABELS_V19.get(str(condition), f"Condition {condition}"),
        "image": image_path.name,
        "path": str(image_path),
        "pattern": pattern_key,
        "pattern_zh": PATTERN_LABELS_ZH_V19[pattern_key],
        "pixel_size_nm": px_nm,
        "pixelsize_txt": str(px_info.txt_path),
        "pixelsize_encoding": px_info.encoding,
        "pixelsize_matched_line": px_info.matched_line,
        "image_height_px": int(gray.shape[0]),
        "image_width_px": int(gray.shape[1]),
        "n_primary_candidates": len(primary_descs),
        "n_candidate_pool": len(candidate_pool),
        "n_selected_structures": n_selected,
        "n_measurements": len(rows),
        "n_rejected_total": len(all_rejects),
        "rejected_reason_summary": format_rejection_summary(all_rejects),
        "annotated_path": str(ann_path),
        "status": "OK" if rows else "CHECK_NO_VALID_MEASUREMENT",
        "note": "",
    }
    status.update({f"search_{k}": v for k, v in search_meta.items()})
    status.update({f"measure_{k}": v for k, v in measure_meta.items() if np.isscalar(v) or isinstance(v, (str, bool))})
    return rows, status, exported_rejects



# ============================================================
# V23 FINAL OVERRIDES
#   * VIA: actual outer edge + light smoothing + morphology QC + strict lattice membership.
#   * SLOT: only complete 4-object columns survive; Y axis uses a TLS line through all
#           four centers in each column and averages the column-line inclinations.
#   * TRENCH: Y-axis points are the midpoints of the OUTER short edges of the two
#             edge-detected rectangular bars (the short edges farther from image center).
#   * Annotation: every measured object is labeled with its X/Y measurement.
# ============================================================
V23_VIA_MIN_SPAN_FACTOR = 0.50
V23_VIA_MAX_SPAN_FACTOR = 2.00
V23_VIA_MIN_EQ_DIAM_FACTOR = 0.74
V23_VIA_MAX_EQ_DIAM_FACTOR = 2.00
V23_VIA_MAX_ASPECT = 1.62
V23_VIA_MIN_SOLIDITY = 0.80
V23_VIA_MIN_CIRCULARITY = 0.50
V23_VIA_MIN_ELLIPSE_IOU = 0.55
V23_VIA_MAX_RADIAL_CV = 0.25
V23_VIA_GRID_SITE_TOL_DIAM_FRAC = 0.34
V23_VIA_OUTER_EDGE_MIN_RADIUS_FACTOR = 0.72
V23_VIA_OUTER_EDGE_MAX_RADIUS_FACTOR = 1.38
V23_VIA40_AFTER_SMOOTH_PASSES = 3
V23_VIA_OTHER_SMOOTH_PASSES = 2
V23_VIA_MIN_GOOD_RAYS_FRAC = 0.70
V23_SLOT_COLUMN_TOL_PITCH_FRAC = 0.34
V23_SLOT_COLUMN_TOL_WIDTH_FACTOR = 0.72
V23_SLOT_REQUIRE_ROWS = (1, 2, 3, 4)
V23_TRENCH_OUTER_DARK_FRAC = 0.28
V23_ANNOTATION_FONT = 0.31


def _circular_smooth_points_v23(pts: np.ndarray, passes: int) -> np.ndarray:
    p = np.asarray(pts, dtype=np.float64).reshape(-1, 2).copy()
    if len(p) < 8:
        return p
    for _ in range(max(0, int(passes))):
        p = 0.25 * np.roll(p, 1, axis=0) + 0.50 * p + 0.25 * np.roll(p, -1, axis=0)
    return p


def _ellipse_iou_v23(contour: np.ndarray) -> float:
    c = np.asarray(contour, np.int32)
    if len(c) < 5:
        return 0.0
    try:
        ell = cv2.fitEllipse(c.astype(np.float32))
    except Exception:
        return 0.0
    pts = c.reshape(-1, 2)
    x0, y0 = np.floor(pts.min(axis=0) - 3).astype(int)
    x1, y1 = np.ceil(pts.max(axis=0) + 3).astype(int)
    w = max(8, int(x1 - x0 + 1)); h = max(8, int(y1 - y0 + 1))
    shift = np.array([[[x0, y0]]], dtype=np.int32)
    cc = c - shift
    m1 = np.zeros((h, w), np.uint8); m2 = np.zeros_like(m1)
    cv2.drawContours(m1, [cc], -1, 255, -1)
    (ecx, ecy), axes, ang = ell
    center = (int(round(ecx - x0)), int(round(ecy - y0)))
    ax = (max(1, int(round(axes[0] / 2.0))), max(1, int(round(axes[1] / 2.0))))
    cv2.ellipse(m2, center, ax, float(ang), 0, 360, 255, -1)
    inter = int(np.count_nonzero((m1 > 0) & (m2 > 0)))
    union = int(np.count_nonzero((m1 > 0) | (m2 > 0)))
    return float(inter / union) if union else 0.0


def _via_shape_metrics_v23(contour: np.ndarray, nominal_diam_px: float) -> Dict[str, float]:
    d = describe_contour(np.asarray(contour, np.int32))
    if d is None:
        return {
            'circularity': 0.0, 'solidity': 0.0, 'aspect': np.inf,
            'ellipse_iou': 0.0, 'radial_cv': np.inf, 'eq_diam_px': 0.0,
        }
    pts = contour.reshape(-1, 2).astype(np.float64)
    center = np.mean(pts, axis=0)
    rr = np.linalg.norm(pts - center, axis=1)
    radial_cv = float(np.std(rr) / max(np.mean(rr), 1e-9)) if len(rr) else np.inf
    eqd = float(math.sqrt(max(0.0, 4.0 * d.area / math.pi)))
    return {
        'circularity': float(d.circularity),
        'solidity': float(d.solidity),
        'aspect': float(d.aspect),
        'ellipse_iou': float(_ellipse_iou_v23(contour)),
        'radial_cv': radial_cv,
        'eq_diam_px': eqd,
        'area_px2': float(d.area),
    }


def _via_shape_qc_v23(
    contour: np.ndarray,
    nominal_diam_px: float,
    stage: str,
    pattern_key: str,
) -> Tuple[bool, Dict[str, float], List[str]]:
    m = _via_shape_metrics_v23(contour, nominal_diam_px)
    reasons: List[str] = []
    # via40-after receives MORE smoothing, not a blanket shape relaxation.  A tiny
    # irregular object must still fail the same physical/array logic as any other via.
    circ_min = V23_VIA_MIN_CIRCULARITY
    sol_min = V23_VIA_MIN_SOLIDITY
    iou_min = V23_VIA_MIN_ELLIPSE_IOU
    rcv_max = V23_VIA_MAX_RADIAL_CV
    if stage == 'after' and pattern_key == 'via40':
        # Only a very small tolerance for genuine etch roughness after smoothing.
        circ_min = 0.48
        sol_min = 0.79
        iou_min = 0.53
        rcv_max = 0.27
    D = max(1e-9, float(nominal_diam_px))
    if m['eq_diam_px'] < V23_VIA_MIN_EQ_DIAM_FACTOR * D:
        reasons.append(f"equivalent diameter {m['eq_diam_px']:.2f}px too small")
    if m['eq_diam_px'] > V23_VIA_MAX_EQ_DIAM_FACTOR * D:
        reasons.append(f"equivalent diameter {m['eq_diam_px']:.2f}px too large")
    if m['circularity'] < circ_min:
        reasons.append(f"circularity {m['circularity']:.3f} < {circ_min:.3f}")
    if m['solidity'] < sol_min:
        reasons.append(f"solidity {m['solidity']:.3f} < {sol_min:.3f}")
    if m['aspect'] > V23_VIA_MAX_ASPECT:
        reasons.append(f"axis ratio {m['aspect']:.3f} > {V23_VIA_MAX_ASPECT:.3f}")
    if m['ellipse_iou'] < iou_min:
        reasons.append(f"ellipse IoU {m['ellipse_iou']:.3f} < {iou_min:.3f}")
    if m['radial_cv'] > rcv_max:
        reasons.append(f"radial CV {m['radial_cv']:.3f} > {rcv_max:.3f}")
    return len(reasons) == 0, m, reasons


def _via_outer_edge_contour_v23(
    gray: np.ndarray,
    center: np.ndarray,
    nominal_diam_px: float,
    nearest_pitch: float,
    stage: str,
    pattern_key: str,
) -> Tuple[Optional[np.ndarray], Dict[str, float]]:
    """Trace the real OUTER via boundary on radial profiles.

    The search begins outside 0.72*nominal radius, so a very dark inner core cannot
    become the measured via.  via40-after is smoothed a little more before edge search.
    """
    H, W = gray.shape
    D = max(7.0, float(nominal_diam_px)); R = 0.5 * D
    c = np.asarray(center, dtype=np.float64)
    sigma = 0.65
    if stage == 'after' and pattern_key == 'via40':
        sigma = max(0.90, 0.030 * D)
    elif stage == 'after':
        sigma = max(0.75, 0.020 * D)
    else:
        sigma = max(0.60, 0.015 * D)
    sm = cv2.GaussianBlur(gray.astype(np.float32), (0, 0), sigmaX=sigma, sigmaY=sigma)

    # Anti-hallucination support check: a real ~40/60 nm via must remain darker than
    # its surroundings through a substantial fraction of the nominal-radius annulus.
    # A tiny black core or a thin star/line may be very dark at the center, but it does
    # not provide this broad dark-disk support and is rejected before radial tracing.
    support_half = int(math.ceil(1.38 * R))
    sx0=max(0,int(math.floor(c[0]-support_half))); sx1=min(W,int(math.ceil(c[0]+support_half+1)))
    sy0=max(0,int(math.floor(c[1]-support_half))); sy1=min(H,int(math.ceil(c[1]+support_half+1)))
    yy0,xx0=np.mgrid[sy0:sy1,sx0:sx1]; rr0=np.sqrt((xx0-c[0])**2+(yy0-c[1])**2)
    core0=rr0<=0.32*R; shoulder0=(rr0>=0.56*R)&(rr0<=0.92*R); outer0=(rr0>=1.14*R)&(rr0<=1.36*R)
    if np.count_nonzero(core0)<5 or np.count_nonzero(shoulder0)<12 or np.count_nonzero(outer0)<12:
        return None, {'good_ray_fraction':0.0,'smooth_sigma_px':sigma,'broad_disk_support':0.0}
    loc=sm[sy0:sy1,sx0:sx1]
    core_g=float(np.median(loc[core0])); shoulder_g=float(np.median(loc[shoulder0])); outer_g=float(np.median(loc[outer0]))
    core_contrast=outer_g-core_g; shoulder_contrast=outer_g-shoulder_g
    support_ratio=shoulder_contrast/max(core_contrast,1e-6)
    support_min=max(0.25,0.060*max(core_contrast,0.0))
    if core_contrast<=0.20 or shoulder_contrast<support_min:
        return None, {'good_ray_fraction':0.0,'smooth_sigma_px':sigma,'broad_disk_support':support_ratio,
                      'core_contrast_gray':core_contrast,'shoulder_contrast_gray':shoulder_contrast}

    r0 = max(1.0, 0.12 * R)
    rmax = V23_VIA_OUTER_EDGE_MAX_RADIUS_FACTOR * R
    if np.isfinite(nearest_pitch) and nearest_pitch > 0:
        rmax = min(rmax, 0.44 * nearest_pitch)
    rmin_edge = V23_VIA_OUTER_EDGE_MIN_RADIUS_FACTOR * R
    if rmax <= rmin_edge + 2.0:
        return None, {'good_ray_fraction': 0.0, 'smooth_sigma_px': sigma}

    nr = max(54, int(math.ceil((rmax - r0) * 4.0)))
    rr = np.linspace(r0, rmax, nr, dtype=np.float32)
    angles = np.linspace(0.0, 2.0 * math.pi, 128, endpoint=False)
    radii = np.full(len(angles), np.nan, dtype=np.float64)
    strengths = np.zeros(len(angles), dtype=np.float64)
    contrasts = np.zeros(len(angles), dtype=np.float64)

    for ia, a in enumerate(angles):
        xs = (c[0] + rr * math.cos(a)).astype(np.float32)
        ys = (c[1] + rr * math.sin(a)).astype(np.float32)
        if xs.min() < 1 or ys.min() < 1 or xs.max() > W - 2 or ys.max() > H - 2:
            continue
        prof = cv2.remap(sm, xs.reshape(1, -1), ys.reshape(1, -1), cv2.INTER_LINEAR,
                         borderMode=cv2.BORDER_REFLECT101).ravel().astype(np.float64)
        core_mask = rr <= 0.46 * R
        bg_mask = rr >= min(1.12 * R, 0.90 * rmax)
        if np.count_nonzero(core_mask) < 3 or np.count_nonzero(bg_mask) < 3:
            continue
        dark = float(np.percentile(prof[core_mask], 55.0))
        bg = float(np.median(prof[bg_mask]))
        contrast = bg - dark
        if not np.isfinite(contrast) or contrast <= 0.18:
            continue
        edge_mask = (rr >= rmin_edge) & (rr <= rmax)
        idx = np.where(edge_mask)[0]
        if len(idx) < 5:
            continue
        grad = np.gradient(prof, rr.astype(np.float64))
        gp = np.maximum(grad[idx], 0.0)
        gscale = float(np.percentile(np.abs(grad[idx]), 85.0)) + 1e-6
        gnorm = np.clip(gp / gscale, 0.0, 3.0)
        prior = np.exp(-0.5 * ((rr[idx] - R) / max(1.0, 0.28 * R)) ** 2)
        half = dark + 0.50 * contrast
        half_score = np.exp(-np.abs(prof[idx] - half) / max(0.18 * contrast, 0.35))
        score = 0.64 * gnorm + 0.21 * prior + 0.15 * half_score
        jloc = int(np.argmax(score)); j = int(idx[jloc])
        # A prior near nominal radius is not allowed to manufacture an edge.  The
        # selected location itself must carry a meaningful positive dark->bright gradient.
        grad_here=float(gp[jloc])
        grad_min=max(0.45,0.070*contrast)
        if float(score[jloc]) < 0.28 or grad_here < grad_min:
            continue
        radii[ia] = float(rr[j])
        strengths[ia] = float(score[jloc])
        contrasts[ia] = contrast

    good = np.isfinite(radii)
    good_frac = float(np.mean(good)) if len(good) else 0.0
    if good_frac < V23_VIA_MIN_GOOD_RAYS_FRAC:
        return None, {'good_ray_fraction': good_frac, 'smooth_sigma_px': sigma}

    # Circular interpolation of missing rays.
    idx_good = np.where(good)[0]
    if len(idx_good) < 3:
        return None, {'good_ray_fraction': good_frac, 'smooth_sigma_px': sigma}
    xg = np.r_[idx_good - len(radii), idx_good, idx_good + len(radii)]
    yg = np.r_[radii[idx_good], radii[idx_good], radii[idx_good]]
    missing = np.where(~good)[0]
    if len(missing):
        radii[missing] = np.interp(missing, xg, yg)

    passes = V23_VIA40_AFTER_SMOOTH_PASSES if (stage == 'after' and pattern_key == 'via40') else V23_VIA_OTHER_SMOOTH_PASSES
    # Smooth radii first; this removes tiny saw-tooth SEM edge noise but not real ovality.
    for _ in range(passes):
        radii = 0.20 * np.roll(radii, 2) + 0.20 * np.roll(radii, 1) + 0.20 * radii + 0.20 * np.roll(radii, -1) + 0.20 * np.roll(radii, -2)

    pts = np.column_stack([c[0] + radii * np.cos(angles), c[1] + radii * np.sin(angles)])
    pts = _circular_smooth_points_v23(pts, 1 if passes >= 2 else 0)
    cont = np.asarray(np.rint(pts), np.int32).reshape(-1, 1, 2)
    return cont, {
        'good_ray_fraction': good_frac,
        'median_edge_radius_px': float(np.median(radii)),
        'median_ray_contrast_gray': float(np.median(contrasts[good])) if np.any(good) else np.nan,
        'median_ray_edge_score': float(np.median(strengths[good])) if np.any(good) else np.nan,
        'smooth_sigma_px': float(sigma),
        'smooth_passes': int(passes),
        'broad_disk_support': float(support_ratio),
        'core_contrast_gray': float(core_contrast),
        'shoulder_contrast_gray': float(shoulder_contrast),
    }


def _validate_via_model_v23(model: Dict[str, object], nominal_diam_px: float) -> bool:
    if not model.get('used', False):
        return False
    la = float(model.get('basis_a_len', np.nan)); lb = float(model.get('basis_b_len', np.nan))
    ang = float(model.get('basis_angle_deg', np.nan)); occ = float(model.get('occupancy', 0.0))
    if not (np.isfinite(la) and np.isfinite(lb) and np.isfinite(ang)):
        return False
    if min(la, lb) < V21_VIA_LATTICE_MIN_PITCH_DIAM_FACTOR * nominal_diam_px:
        return False
    if max(la, lb) > V21_VIA_LATTICE_MAX_PITCH_DIAM_FACTOR * nominal_diam_px:
        return False
    if not (V21_VIA_LATTICE_MIN_BASIS_ANGLE <= ang <= V21_VIA_LATTICE_MAX_BASIS_ANGLE):
        return False
    return int(model.get('supported_sites', 0)) >= 4 and occ >= 0.20


def _select_and_measure_vias_v19(
    gray: np.ndarray,
    pool: List[ShapeDesc],
    px_nm: float,
    image_name: str,
    stage: str,
    condition: str,
    pattern_key: str,
) -> Tuple[List[dict], List[ShapeDesc], List[dict], Dict[str, object]]:
    nominal_nm, D = _via_nominal_diam_px_v21(pattern_key, px_nm)
    score, raw, best_diam, resp_meta = _via_scale_response_v21(gray, D)
    raw_floor = float(resp_meta['raw_contrast_floor'])
    scale_peaks = _via_scale_peaks_v21(gray, score, raw, best_diam, D, raw_floor)
    legacy = _via_legacy_size_candidates_v21(pool, score, raw, D, raw_floor, gray.shape)
    candidates = _merge_via_center_candidates_v21(scale_peaks + legacy, D)
    model, model_seeds = _choose_via_lattice_v21(candidates, D)

    # If the first fit fails, retry only with strong scale-matched centers.  We never
    # keep free-floating fallback vias: every final object must be explainable by a grid.
    if not _validate_via_model_v23(model, D):
        strong = [d for d in candidates if float(getattr(d, '_v21_score', 0.0)) >= 0.34 and float(getattr(d, '_v21_raw', 0.0)) >= 0.90 * raw_floor]
        model2, seeds2 = _choose_via_lattice_v21(strong, D)
        if _validate_via_model_v23(model2, D):
            model, model_seeds = model2, seeds2
    rejects: List[dict] = []
    if not _validate_via_model_v23(model, D):
        for d in candidates[:60]:
            if float(getattr(d, '_v21_score', 0.0)) >= 0.34:
                rejects.append(_make_reject_from_contour(
                    d.contour, 'via_v23_grid_qc', 'via_no_reliable_rectangular_lattice',
                    'candidate may be dark/round, but no reliable rectangular via lattice was established',
                    desc=d, scale_match_score=float(getattr(d, '_v21_score', np.nan))))
        axis = _fit_via_bottom_row_axes_v19([], gray.shape)
        meta = {'via_v23_final_count': 0, 'via_v23_grid_required': True, 'via_v23_lattice_used': False, **resp_meta}
        return [], [], rejects, {'axis': axis, 'meta': meta, 'lattice_model': model}

    lens = [float(model.get('basis_a_len', np.nan)), float(model.get('basis_b_len', np.nan))]
    lens = [x for x in lens if np.isfinite(x) and x > 0]
    nearest_pitch = min(lens) if lens else np.nan
    margin = max(3.0, V21_VIA_BORDER_DIAM_FRAC * D)
    sites = lattice_sites_v14(model, gray.shape, 1, 'via')
    pre: List[ShapeDesc] = []
    recovered = 0
    for site in sites:
        xy0 = np.array([float(site['x']), float(site['y'])], dtype=float)
        if xy0[0] <= margin or xy0[1] <= margin or xy0[0] >= gray.shape[1]-1-margin or xy0[1] >= gray.shape[0]-1-margin:
            continue
        xy, sc, rw = _via_refine_site_center_v21(score, raw, xy0, D, nearest_pitch)
        internal = bool(site.get('is_internal', False))
        # Slightly lower signal gate for blurry via40-after, but grid membership remains hard.
        sc_min = V21_VIA_INTERNAL_SITE_SCORE_MIN if internal else V21_VIA_EXTERNAL_SITE_SCORE_MIN
        rw_min = raw_floor if internal else V21_VIA_EXTERNAL_RAW_FACTOR * raw_floor
        if stage == 'after' and pattern_key == 'via40':
            sc_min *= 0.88; rw_min *= 0.86
        if sc < sc_min or rw < rw_min:
            continue
        d = _make_circle_desc_v21(xy, D, 'v23_lattice_site')
        if d is None:
            continue
        setattr(d, '_v21_score', float(sc)); setattr(d, '_v21_raw', float(rw))
        setattr(d, '_v23_lattice_i', int(site['i'])); setattr(d, '_v23_lattice_j', int(site['j']))
        setattr(d, '_v23_lattice_internal', bool(internal))
        setattr(d, '_v23_lattice_pred_xy', xy0.copy())
        pre.append(d)
        if not any(float(np.linalg.norm(xy - s.center)) <= 0.30 * D for s in model_seeds):
            recovered += 1

    pre = _merge_via_center_candidates_v21(pre, D)
    # External expanded rows/columns survive only when they have repeated support.
    if pre:
        ci = Counter(int(getattr(d, '_v23_lattice_i', 0)) for d in pre)
        cj = Counter(int(getattr(d, '_v23_lattice_j', 0)) for d in pre)
        imin, imax = int(model['i_min']), int(model['i_max']); jmin, jmax = int(model['j_min']), int(model['j_max'])
        pre2 = []
        for d in pre:
            i = int(getattr(d, '_v23_lattice_i', 0)); j = int(getattr(d, '_v23_lattice_j', 0))
            ext_i_ok = (imin <= i <= imax) or ci[i] >= 3
            ext_j_ok = (jmin <= j <= jmax) or cj[j] >= 3
            if ext_i_ok and ext_j_ok:
                pre2.append(d)
            else:
                rejects.append(_make_reject_from_contour(
                    d.contour, 'via_v23_grid_qc', 'isolated_external_lattice_site',
                    'expanded lattice site lacks repeated row/column support', desc=d,
                    lattice_i=i, lattice_j=j))
        pre = pre2

    # Trace actual outer edges, smooth only enough to suppress saw-tooth pixels, then
    # require round/elliptical morphology.  No ideal ellipse fallback is accepted.
    survivors: List[ShapeDesc] = []
    edge_meta_by_id: Dict[int, Dict[str, float]] = {}
    shape_meta_by_id: Dict[int, Dict[str, float]] = {}
    for d in pre:
        pred = np.asarray(getattr(d, '_v23_lattice_pred_xy', d.center), dtype=float)
        grid_dist = float(np.linalg.norm(d.center - pred))
        grid_tol = max(2.5, V23_VIA_GRID_SITE_TOL_DIAM_FRAC * D)
        if np.isfinite(nearest_pitch):
            grid_tol = min(grid_tol, 0.22 * nearest_pitch)
        if grid_dist > grid_tol:
            rejects.append(_make_reject_from_contour(
                d.contour, 'via_v23_grid_qc', 'via_center_off_lattice',
                f'center is {grid_dist:.2f}px from predicted grid site (> {grid_tol:.2f}px)',
                desc=d, lattice_distance_px=grid_dist))
            continue
        cont, em = _via_outer_edge_contour_v23(gray, d.center, D, nearest_pitch, stage, pattern_key)
        if cont is None:
            rejects.append(_make_reject_from_contour(
                d.contour, 'via_v23_edge_qc', 'via_outer_edge_not_reliable',
                'outer dark-to-bright edge could not be traced on enough rays', desc=d,
                good_ray_fraction=em.get('good_ray_fraction', np.nan)))
            continue
        ok, smet, reasons = _via_shape_qc_v23(cont, D, stage, pattern_key)
        if not ok:
            rejects.append(_make_reject_from_contour(
                cont, 'via_v23_shape_qc', 'via_irregular_shape', '; '.join(reasons), desc=d, **smet))
            continue
        d.contour = cont.astype(np.int32)
        ad = describe_contour(d.contour)
        if ad is not None:
            # Keep the lattice center for grid coordinates, but morphology comes from real edge.
            d.area = ad.area; d.perimeter = ad.perimeter; d.circularity = ad.circularity; d.solidity = ad.solidity
            d.u = ad.u; d.v = ad.v; d.major_min = ad.major_min; d.major_max = ad.major_max
            d.minor_min = ad.minor_min; d.minor_max = ad.minor_max; d.length_px = ad.length_px; d.width_px = ad.width_px
            d.angle_from_x_deg = ad.angle_from_x_deg
        survivors.append(d); edge_meta_by_id[id(d)] = em; shape_meta_by_id[id(d)] = smet

    # Final via X/Y axes from the bottom row of only the valid grid/shape objects.
    axis = _fit_via_bottom_row_axes_v19(survivors, gray.shape)
    xa = np.asarray(axis['x_axis'], dtype=float); ya = np.asarray(axis['y_axis_up'], dtype=float); origin = np.asarray(axis['origin'], dtype=float)
    final: List[ShapeDesc] = []; rows: List[dict] = []
    bottom_set = set(axis.get('bottom_indices', []))
    for idx, d in enumerate(survivors):
        wx = float(_projection_span_v19(d.contour, xa)); hy = float(_projection_span_v19(d.contour, ya))
        w_nm = wx * px_nm; h_nm = hy * px_nm
        if not (V23_VIA_MIN_SPAN_FACTOR * nominal_nm <= w_nm <= V23_VIA_MAX_SPAN_FACTOR * nominal_nm and
                V23_VIA_MIN_SPAN_FACTOR * nominal_nm <= h_nm <= V23_VIA_MAX_SPAN_FACTOR * nominal_nm):
            rejects.append(_make_reject_from_contour(
                d.contour, 'via_v23_physical_qc', 'via_final_span_out_of_range',
                f'actual X/Y spans {w_nm:.2f}/{h_nm:.2f}nm inconsistent with nominal {nominal_nm:.0f}nm via',
                desc=d, x_width_nm=w_nm, y_height_nm=h_nm))
            continue
        em = edge_meta_by_id.get(id(d), {}); smet = shape_meta_by_id.get(id(d), {})
        final.append(d)
        rows.append({
            'stage': stage, 'stage_zh': STAGE_LABELS_ZH_V19[stage], 'condition': condition,
            'condition_label': CONDITION_LABELS_V19.get(str(condition), f'Condition {condition}'),
            'image': image_name, 'pattern': pattern_key, 'pattern_zh': PATTERN_LABELS_ZH_V19[pattern_key],
            'object_id': len(final), 'used_for_statistics': True, 'is_axis_bottom_row': idx in bottom_set,
            'x_width_nm': w_nm, 'y_height_nm': h_nm, 'tip_gap_y_nm': np.nan,
            'x_width_px': wx, 'y_height_px': hy,
            'center_image_x_px': float(d.center[0]), 'center_image_y_px': float(d.center[1]),
            'center_axis_x_px': float((d.center-origin)@xa), 'center_axis_y_up_px': float((d.center-origin)@ya),
            'circularity': float(smet.get('circularity', d.circularity)), 'solidity': float(smet.get('solidity', d.solidity)),
            'axis_ratio': float(smet.get('aspect', d.aspect)), 'via_ellipse_iou': float(smet.get('ellipse_iou', np.nan)),
            'via_radial_cv': float(smet.get('radial_cv', np.nan)), 'via_equivalent_diameter_nm': float(smet.get('eq_diam_px', np.nan))*px_nm,
            'selection_source': _source_of(d), 'pixel_size_nm': px_nm,
            'via_nominal_diameter_nm': nominal_nm, 'via_nominal_diameter_px': D,
            'via_scale_match_score': float(getattr(d, '_v21_score', np.nan)),
            'via_scale_raw_contrast_gray': float(getattr(d, '_v21_raw', np.nan)),
            'via_lattice_i': int(getattr(d, '_v23_lattice_i', 0)), 'via_lattice_j': int(getattr(d, '_v23_lattice_j', 0)),
            'via_lattice_site_distance_px': float(np.linalg.norm(d.center - np.asarray(getattr(d, '_v23_lattice_pred_xy', d.center)))),
            'via_axis_angle_from_image_x_deg': axis.get('row_angle_deg', np.nan), 'via_axis_fit_rms_px': axis.get('row_fit_rms_px', np.nan),
            'estimated_center_pitch_px': nearest_pitch, 'x_edge_method': 'v23_actual_outer_edge', 'y_edge_method': 'v23_actual_outer_edge',
            'via_good_ray_fraction': em.get('good_ray_fraction', np.nan), 'via_edge_smooth_sigma_px': em.get('smooth_sigma_px', np.nan),
            'via_edge_smooth_passes': em.get('smooth_passes', np.nan),
        })

    # Refit axis after final physical filter, update object IDs and bottom flags.
    if final:
        axis = _fit_via_bottom_row_axes_v19(final, gray.shape)
        xa = np.asarray(axis['x_axis'], dtype=float); ya = np.asarray(axis['y_axis_up'], dtype=float); origin = np.asarray(axis['origin'], dtype=float)
        bset = set(axis.get('bottom_indices', []))
        for i, (d, r) in enumerate(zip(final, rows)):
            r['object_id'] = i + 1; r['is_axis_bottom_row'] = i in bset
            r['center_axis_x_px'] = float((d.center-origin)@xa); r['center_axis_y_up_px'] = float((d.center-origin)@ya)
            r['x_width_px'] = float(_projection_span_v19(d.contour, xa)); r['y_height_px'] = float(_projection_span_v19(d.contour, ya))
            r['x_width_nm'] = r['x_width_px'] * px_nm; r['y_height_nm'] = r['y_height_px'] * px_nm
            r['via_axis_angle_from_image_x_deg'] = axis.get('row_angle_deg', np.nan); r['via_axis_fit_rms_px'] = axis.get('row_fit_rms_px', np.nan)

    meta = {
        'array_assist_used': True, 'array_lattice_type': model.get('lattice_type', 'rectangular'),
        'array_expected_positions': len(sites), 'array_recovered_count': recovered,
        'array_occupancy_seed': float(model.get('occupancy', np.nan)),
        'array_basis_a_length_px': float(model.get('basis_a_len', np.nan)), 'array_basis_b_length_px': float(model.get('basis_b_len', np.nan)),
        'array_basis_angle_deg': float(model.get('basis_angle_deg', np.nan)), 'via_bottom_row_count': len(axis.get('bottom_indices', [])),
        'via_axis_angle_deg': float(axis.get('row_angle_deg', np.nan)), 'via_axis_fit_rms_px': float(axis.get('row_fit_rms_px', np.nan)),
        'via_axis_fallback': bool(axis.get('fallback', False)), 'estimated_center_pitch_px': nearest_pitch,
        'via_v23_nominal_diameter_nm': nominal_nm, 'via_v23_nominal_diameter_px': D,
        'via_v23_grid_required': True, 'via_v23_lattice_used': True, 'via_v23_pre_qc_count': len(pre), 'via_v23_final_count': len(final),
        **resp_meta,
    }
    return rows, final, rejects, {'axis': axis, 'meta': meta, 'lattice_model': model}


def _label_box_v23(out: np.ndarray, xy: Tuple[int, int], text: str, color=(255,255,255), font: float = V23_ANNOTATION_FONT) -> None:
    H, W = out.shape[:2]
    x = int(np.clip(xy[0], 1, W-2)); y = int(np.clip(xy[1], 10, H-2))
    (tw, th), base = cv2.getTextSize(text, cv2.FONT_HERSHEY_SIMPLEX, font, 1)
    x = min(x, max(1, W - tw - 3)); y = max(th + 2, min(y, H - base - 2))
    # Thin black halo preserves readability without covering neighboring SEM patterns.
    cv2.putText(out, text, (x, y), cv2.FONT_HERSHEY_SIMPLEX, font, (0,0,0), 2, cv2.LINE_AA)
    cv2.putText(out, text, (x, y), cv2.FONT_HERSHEY_SIMPLEX, font, color, 1, cv2.LINE_AA)


def _annotate_vias_v19(gray: np.ndarray, selected: List[ShapeDesc], rows: List[dict], axis_bundle: Dict[str, object]) -> np.ndarray:
    out = to_color(gray); axis = axis_bundle['axis']; bottom = set(axis.get('bottom_indices', []))
    for i, (d, r) in enumerate(zip(selected, rows)):
        color = (0,255,255) if i in bottom else (0,210,0)
        cv2.drawContours(out, [np.rint(d.contour).astype(np.int32)], -1, color, 2 if i in bottom else 1)
        c = tuple(np.round(d.center).astype(int)); cv2.circle(out, c, 2, color, -1)
        _label_box_v23(out, (c[0]-8, c[1]-6), f"V{i+1}", color, 0.28)
        _label_box_v23(out, (c[0]-18, c[1]+11), f"{r['x_width_nm']:.1f}x{r['y_height_nm']:.1f}", color, 0.27)
    if selected:
        _draw_axes_v22(out, axis)
        mw = float(np.mean([r['x_width_nm'] for r in rows])) if rows else np.nan
        mh = float(np.mean([r['y_height_nm'] for r in rows])) if rows else np.nan
        nominal = float(rows[0].get('via_nominal_diameter_nm', np.nan)) if rows else np.nan
        text = f"used={len(rows)} nominal={nominal:.0f}nm Xmean={mw:.2f}nm Ymean={mh:.2f}nm; shape+grid QC"
        cv2.rectangle(out,(4,4),(min(gray.shape[1]-4,790),29),(0,0,0),-1)
        cv2.putText(out,text,(8,22),cv2.FONT_HERSHEY_SIMPLEX,0.45,(255,255,255),1,cv2.LINE_AA)
    return out


def _trench_outer_short_edge_midpoint_v23(
    gray: np.ndarray, xl: float, xr: float, tip_y: float, is_top: bool
) -> Tuple[np.ndarray, Dict[str, float]]:
    """Midpoint of the detected rectangle short edge farthest from image center."""
    H, W = gray.shape
    x0 = max(1, int(math.floor(min(xl, xr)))); x1 = min(W-2, int(math.ceil(max(xl, xr))))
    if x1 <= x0 + 1:
        return np.array([0.5*(xl+xr), 0.0 if is_top else H-1.0]), {'outer_edge_fallback': 1.0}
    # Median inside the bar vs adjacent side background gives a Y-wise darkness score.
    pad = max(2, int(round(0.8 * (x1-x0+1))))
    left = gray[:, max(0,x0-pad):x0].astype(np.float32)
    right = gray[:, x1+1:min(W,x1+1+pad)].astype(np.float32)
    body = gray[:, x0:x1+1].astype(np.float32)
    body_p = np.median(body, axis=1)
    side_parts = []
    if left.size: side_parts.append(left)
    if right.size: side_parts.append(right)
    if side_parts:
        bg_p = np.median(np.concatenate(side_parts, axis=1), axis=1)
    else:
        bg_p = np.full(H, float(np.median(gray)), dtype=np.float32)
    dark = _smooth_1d_v20((bg_p-body_p).astype(np.float64), 1.0)
    ty = int(np.clip(round(tip_y), 1, H-2))
    if is_top:
        ref_slice = dark[max(0, ty-max(5,int(0.16*H))):max(1,ty-2)]
        ref = float(np.percentile(ref_slice, 55)) if len(ref_slice) else float(np.max(dark[:ty]))
        thr = max(0.18, V23_TRENCH_OUTER_DARK_FRAC * max(ref, 0.5))
        y = ty - 2; misses = 0; outer = y
        while y >= 0:
            if dark[y] >= thr:
                outer = y; misses = 0
            else:
                misses += 1
                if misses > 2: break
            y -= 1
    else:
        ref_slice = dark[min(H,ty+2):min(H,ty+max(5,int(0.16*H)))]
        ref = float(np.percentile(ref_slice, 55)) if len(ref_slice) else float(np.max(dark[ty:]))
        thr = max(0.18, V23_TRENCH_OUTER_DARK_FRAC * max(ref, 0.5))
        y = ty + 2; misses = 0; outer = y
        while y < H:
            if dark[y] >= thr:
                outer = y; misses = 0
            else:
                misses += 1
                if misses > 2: break
            y += 1
    return np.array([0.5*(xl+xr), float(np.clip(outer,0,H-1))]), {'outer_dark_ref': ref, 'outer_dark_threshold': thr}


def _select_and_measure_trenches_v19(
    gray: np.ndarray, pool: List[ShapeDesc], px_nm: float, image_name: str,
    stage: str, condition: str, pattern_key: str,
) -> Tuple[List[dict], List[dict], List[dict], Dict[str, object]]:
    pairs, meta = _detect_trench_pairs_profile_v20(gray); H,W = gray.shape
    if pairs:
        xs=sorted([p['common_x'] for p in pairs]); diffs=np.diff(xs)
        approx_pitch=float(np.median(diffs)) if len(diffs) else float(W/max(4,V20_TRENCH_EXPECTED_GROUPS_HINT))
    else:
        approx_pitch=np.nan
    axis_pairs=[]
    for p in pairs:
        txl,txr,tw,_=_single_trench_body_x_edges_v22(gray,int(round(p['common_x'])),approx_pitch,V20_TRENCH_TOP_BODY_RANGE)
        bxl,bxr,bw,_=_single_trench_body_x_edges_v22(gray,int(round(p['common_x'])),approx_pitch,V20_TRENCH_BOTTOM_BODY_RANGE)
        top_outer, tom = _trench_outer_short_edge_midpoint_v23(gray,txl,txr,p['top_tip_y_px'],True)
        bot_outer, bom = _trench_outer_short_edge_midpoint_v23(gray,bxl,bxr,p['bottom_tip_y_px'],False)
        # The Y-axis line for each column is now defined by the OUTER short-edge midpoints.
        p.update({
            'top_x_left_px':txl,'top_x_right_px':txr,'bottom_x_left_px':bxl,'bottom_x_right_px':bxr,
            'top_center_x_px':float(top_outer[0]),'bottom_center_x_px':float(bot_outer[0]),
            'top_center_y_px':float(top_outer[1]),'bottom_center_y_px':float(bot_outer[1]),
            'top_outer_short_mid_x_px':float(top_outer[0]),'top_outer_short_mid_y_px':float(top_outer[1]),
            'bottom_outer_short_mid_x_px':float(bot_outer[0]),'bottom_outer_short_mid_y_px':float(bot_outer[1]),
            'top_width_px_v22':tw,'bottom_width_px_v22':bw,
            'top_outer_dark_ref':tom.get('outer_dark_ref',np.nan),'bottom_outer_dark_ref':bom.get('outer_dark_ref',np.nan),
        })
        axis_pairs.append((top_outer,bot_outer))
    axis=_axes_from_top_bottom_points_v22(axis_pairs,gray.shape)
    xa=np.asarray(axis['x_axis'],float); yd=np.asarray(axis['y_axis_down'],float); yu=np.asarray(axis['y_axis_up'],float); origin=np.asarray(axis['origin'],float)
    for p in pairs:
        # Tip points remain the inward-facing short-edge midpoints for gap measurement.
        top_tip=np.array([0.5*(p['top_x_left_px']+p['top_x_right_px']),p['top_tip_y_px']],float)
        bot_tip=np.array([0.5*(p['bottom_x_left_px']+p['bottom_x_right_px']),p['bottom_tip_y_px']],float)
        mid=0.5*(top_tip+bot_tip)
        p['tip_gap_axis_px']=abs(float((bot_tip-top_tip)@yd)); p['axis_center_x_px']=float((mid-origin)@xa); p['axis_center_y_up_px']=float((mid-origin)@yu)
        # User definition: trench width is the arithmetic mean of the upper/lower
        # edge-detected rectangle widths.  Do NOT project this width onto the fitted X axis.
        p['x_width_axis_px']=0.5*(p.get('top_width_px_v22',p['width_px'])+p.get('bottom_width_px_v22',p['width_px']))
    pairs=sorted(pairs,key=lambda p:p.get('axis_center_x_px',p['common_x']))
    rows=[]
    for oid,p in enumerate(pairs,1):
        rows.append({
            'stage':stage,'stage_zh':STAGE_LABELS_ZH_V19[stage],'condition':condition,'condition_label':CONDITION_LABELS_V19.get(str(condition),f'Condition {condition}'),
            'image':image_name,'pattern':pattern_key,'pattern_zh':PATTERN_LABELS_ZH_V19[pattern_key],'object_id':oid,'used_for_statistics':True,'is_axis_bottom_row':False,
            'x_width_nm':p['x_width_axis_px']*px_nm,'y_height_nm':np.nan,'tip_gap_y_nm':p['tip_gap_axis_px']*px_nm,'x_width_px':p['x_width_axis_px'],'y_height_px':np.nan,
            'trench_mean_rectangle_width_nm':p['x_width_axis_px']*px_nm,'trench_mean_rectangle_width_px':p['x_width_axis_px'],
            'center_image_x_px':0.5*(p['top_x_left_px']+p['top_x_right_px']+p['bottom_x_left_px']+p['bottom_x_right_px'])/2.0,
            'center_image_y_px':0.5*(p['top_tip_y_px']+p['bottom_tip_y_px']),'center_axis_x_px':p['axis_center_x_px'],'center_axis_y_up_px':p['axis_center_y_up_px'],
            'top_trench_width_nm':p.get('top_width_px_v22',p['width_px'])*px_nm,'bottom_trench_width_nm':p.get('bottom_width_px_v22',p['width_px'])*px_nm,
            'pair_x_offset_nm':abs(p['bottom_outer_short_mid_x_px']-p['top_outer_short_mid_x_px'])*px_nm,'pair_x_overlap_frac':1.0,
            'width_mismatch_frac':abs(p.get('top_width_px_v22',p['width_px'])-p.get('bottom_width_px_v22',p['width_px']))/max(1e-9,0.5*(p.get('top_width_px_v22',p['width_px'])+p.get('bottom_width_px_v22',p['width_px']))),
            'pair_score':p['score'],'pixel_size_nm':px_nm,'trench_detector':'V20_gray_profile+V24_y_axis_gap+mean_rectangle_width',
            'top_tip_contrast_gray':p.get('top_tip_contrast',np.nan),'bottom_tip_contrast_gray':p.get('bottom_tip_contrast',np.nan),
            'top_tip_gradient':p.get('top_tip_gradient',np.nan),'bottom_tip_gradient':p.get('bottom_tip_gradient',np.nan),
            'top_axis_point_x_px':p['top_outer_short_mid_x_px'],'top_axis_point_y_px':p['top_outer_short_mid_y_px'],
            'bottom_axis_point_x_px':p['bottom_outer_short_mid_x_px'],'bottom_axis_point_y_px':p['bottom_outer_short_mid_y_px'],
            'x_axis_angle_from_image_x_deg':axis['x_axis_angle_from_image_x_deg'],'y_axis_down_angle_from_image_x_deg':axis['y_axis_down_angle_from_image_x_deg'],
            'axis_pair_count':axis['pair_count'],'axis_fallback':axis['fallback'],
        })
    meta=dict(meta); meta.update({'v23_axis_pair_count':axis['pair_count'],'v23_axis_fallback':axis['fallback'],
        'v23_x_axis_angle_from_image_x_deg':axis['x_axis_angle_from_image_x_deg'],'v23_y_axis_down_angle_from_image_x_deg':axis['y_axis_down_angle_from_image_x_deg'],
        'v23_trench_axis_point_definition':'outer_short_edge_midpoints'})
    for p in pairs: p['_v22_axis']=axis
    return rows,pairs,[],meta


def _annotate_trenches_v19(gray: np.ndarray, pairs: List[dict], rows: List[dict]) -> np.ndarray:
    out=to_color(gray); H,W=gray.shape
    for i,(p,r) in enumerate(zip(pairs,rows),1):
        txl=int(round(p['top_x_left_px'])); txr=int(round(p['top_x_right_px'])); bxl=int(round(p['bottom_x_left_px'])); bxr=int(round(p['bottom_x_right_px']))
        ty=int(round(p['top_tip_y_px'])); by=int(round(p['bottom_tip_y_px'])); toy=int(round(p['top_outer_short_mid_y_px'])); boy=int(round(p['bottom_outer_short_mid_y_px']))
        cv2.rectangle(out,(max(0,txl),max(0,min(toy,ty))),(min(W-1,txr),min(H-1,max(toy,ty))),(0,210,0),1)
        cv2.rectangle(out,(max(0,bxl),max(0,min(by,boy))),(min(W-1,bxr),min(H-1,max(by,boy))),(0,210,0),1)
        top_axis=(int(round(p['top_outer_short_mid_x_px'])),toy); bot_axis=(int(round(p['bottom_outer_short_mid_x_px'])),boy)
        cv2.circle(out,top_axis,3,(255,0,255),-1); cv2.circle(out,bot_axis,3,(255,0,255),-1)
        t=(int(round(0.5*(p['top_x_left_px']+p['top_x_right_px']))),ty); b=(int(round(0.5*(p['bottom_x_left_px']+p['bottom_x_right_px']))),by)
        cv2.arrowedLine(out,t,b,(255,255,255),1,cv2.LINE_AA,tipLength=0.12); cv2.arrowedLine(out,b,t,(255,255,255),1,cv2.LINE_AA,tipLength=0.12)
        _label_box_v23(out,(t[0]-12,ty-6),f"T{i} W{r['x_width_nm']:.1f}",(255,255,255),0.28)
        _label_box_v23(out,(b[0]-13,by+13),f"Y{r['tip_gap_y_nm']:.1f}",(255,255,255),0.28)
    if pairs and '_v22_axis' in pairs[0]: _draw_axes_v22(out,pairs[0]['_v22_axis'])
    if rows:
        text=f"trench n={len(rows)} Wmean={np.mean([r['x_width_nm'] for r in rows]):.2f}nm Ygap mean={np.mean([r['tip_gap_y_nm'] for r in rows]):.2f}nm; magenta=axis points"
    else: text='trench: no valid opposing pairs'
    cv2.rectangle(out,(4,4),(min(W-4,820),29),(0,0,0),-1); cv2.putText(out,text,(8,22),cv2.FONT_HERSHEY_SIMPLEX,0.43,(255,255,255),1,cv2.LINE_AA)
    return out


def _slot_prelim_axis_v23(paired: List[Tuple[ShapeDesc,dict]], image_shape: Tuple[int,int]) -> Dict[str,object]:
    axis,_ = _slot_axis_pairs_v22(paired,image_shape)
    return axis


def _slot_pitch_in_axis_v23(paired: List[Tuple[ShapeDesc,dict]], axis: Dict[str,object]) -> float:
    xa=np.asarray(axis['x_axis'],float); origin=np.asarray(axis['origin'],float)
    diffs=[]
    for rid in V23_SLOT_REQUIRE_ROWS:
        us=sorted(float((d.center-origin)@xa) for d,r in paired if int(r.get('slot_row',0))==rid)
        if len(us)>=2:
            dd=np.diff(us); diffs.extend(float(x) for x in dd if x>=4.0)
    if not diffs: return np.nan
    arr=np.asarray(diffs,float); med=float(np.median(arr))
    core=arr[arr<=1.55*med] if med>0 else arr
    return float(np.median(core)) if len(core) else med


def _slot_complete_columns_once_v23(
    paired: List[Tuple[ShapeDesc,dict]], axis: Dict[str,object]
) -> Tuple[List[Dict[int,Tuple[ShapeDesc,dict]]], float, float]:
    xa=np.asarray(axis['x_axis'],float); origin=np.asarray(axis['origin'],float)
    pitch=_slot_pitch_in_axis_v23(paired,axis)
    widths=[]
    for d,_ in paired:
        widths.append(float(_projection_span_v19(d.contour,xa)))
    medw=float(np.median(widths)) if widths else 5.0
    if not np.isfinite(pitch) or pitch<=0:
        pitch=max(2.8*medw,12.0)
    tol=max(3.0,min(V23_SLOT_COLUMN_TOL_PITCH_FRAC*pitch,max(4.0,V23_SLOT_COLUMN_TOL_WIDTH_FACTOR*medw)))
    rowmaps={rid:[] for rid in V23_SLOT_REQUIRE_ROWS}
    for d,r in paired:
        rid=int(r.get('slot_row',0))
        if rid in rowmaps:
            rowmaps[rid].append((float((d.center-origin)@xa),d,r))
    for rid in rowmaps: rowmaps[rid].sort(key=lambda z:z[0])
    anchors=rowmaps[1]
    used={rid:set() for rid in V23_SLOT_REQUIRE_ROWS}
    cols=[]
    for anchor_k,(ua,da,ra) in enumerate(anchors):
        if anchor_k in used[1]:
            continue
        col={1:(da,ra)}; choices=[]; valid=True
        for rid in (2,3,4):
            opts=[(abs(u-ua),k,u,d,r) for k,(u,d,r) in enumerate(rowmaps[rid]) if k not in used[rid] and abs(u-ua)<=tol]
            if not opts:
                valid=False; break
            choices.append((rid,min(opts,key=lambda z:z[0])))
        if not valid: continue
        used[1].add(anchor_k)
        for rid,choice in choices:
            _,k,_,d,r=choice; used[rid].add(k); col[rid]=(d,r)
        cols.append(col)
    return cols,float(pitch),float(tol)


def _slot_axis_from_columns_v23(columns: List[Dict[int,Tuple[ShapeDesc,dict]]], image_shape: Tuple[int,int]) -> Dict[str,object]:
    pp=[]
    for col in columns:
        pts=np.asarray([col[r][0].center for r in V23_SLOT_REQUIRE_ROWS],dtype=float)
        if len(pts)<4: continue
        o,u=_fit_tls_line_v19(pts)
        if float(u @ (pts[-1]-pts[0])) < 0: u=-u
        tt=(pts-o)@u
        pp.append((o+u*float(np.min(tt)),o+u*float(np.max(tt))))
    axis=_axes_from_top_bottom_points_v22(pp,image_shape)
    axis['slot_complete_column_count']=len(columns)
    return axis


def _select_and_measure_bottom_slots_v19(
    gray: np.ndarray, pool: List[ShapeDesc], px_nm: float, image_name: str,
    stage: str, condition: str, pattern_key: str,
) -> Tuple[List[dict], List[ShapeDesc], List[ShapeDesc], List[dict], Dict[str, object]]:
    old_rows,selected_raw,_,rejects,meta=measure_slot_v18(gray,pool,px_nm,image_name,stage)
    paired=[]
    for d,row in zip(selected_raw,old_rows):
        if _bbox_complete_v19(d,gray.shape,0.10): paired.append((d,row))
        else: rejects.append(_make_reject_from_contour(d.contour,'slot_v23_border_qc','slot_near_or_crossing_image_edge','slot may be incomplete at image edge',desc=d))
    if not paired:
        meta=dict(meta); meta.update({'v23_complete_slot_columns':0,'v23_slot_grid_required':True})
        return [],[],[],rejects,meta

    prelim=_slot_prelim_axis_v23(paired,gray.shape)
    cols,pitch,tol=_slot_complete_columns_once_v23(paired,prelim)
    if cols:
        axis=_slot_axis_from_columns_v23(cols,gray.shape)
        # One refinement pass: complete-column membership is reevaluated in the final axis.
        cols2,pitch2,tol2=_slot_complete_columns_once_v23(paired,axis)
        if len(cols2)>=len(cols) or len(cols2)>=2:
            cols,pitch,tol=cols2,pitch2,tol2; axis=_slot_axis_from_columns_v23(cols,gray.shape)
    else:
        axis=_axes_from_top_bottom_points_v22([],gray.shape)

    if not cols:
        for d,r in paired:
            rejects.append(_make_reject_from_contour(d.contour,'slot_v23_column_qc','slot_not_in_complete_4row_column','candidate is not part of a complete column containing rows 1,2,3,4',desc=d,slot_row=int(r.get('slot_row',0))))
        meta=dict(meta); meta.update({'v23_complete_slot_columns':0,'v23_slot_grid_required':True,'v23_slot_column_pitch_px':pitch,'v23_slot_column_tol_px':tol})
        return [],[],[],rejects,meta

    xa=np.asarray(axis['x_axis'],float); yu=np.asarray(axis['y_axis_up'],float); origin=np.asarray(axis['origin'],float)
    cols=sorted(cols,key=lambda col:float((np.mean([col[r][0].center for r in V23_SLOT_REQUIRE_ROWS],axis=0)-origin)@xa))
    selected_all=[]; bottom=[]; used_ids=set(); rows=[]
    for ci,col in enumerate(cols,1):
        for rid in V23_SLOT_REQUIRE_ROWS:
            d,old=col[rid]; selected_all.append(d); used_ids.add(id(d))
            setattr(d,'_v23_slot_row',rid); setattr(d,'_v23_slot_col',ci); setattr(d,'_v22_slot_axis',axis)
        d,old=col[4]; bottom.append(d)
        w_px=float(_projection_span_v19(d.contour,xa)); h_px=float(_projection_span_v19(d.contour,yu)); reg=slot_regularity_metrics(d)
        rows.append({
            'stage':stage,'stage_zh':STAGE_LABELS_ZH_V19[stage],'condition':condition,'condition_label':CONDITION_LABELS_V19.get(str(condition),f'Condition {condition}'),
            'image':image_name,'pattern':pattern_key,'pattern_zh':PATTERN_LABELS_ZH_V19[pattern_key],'object_id':ci,'used_for_statistics':True,'is_axis_bottom_row':True,
            'x_width_nm':w_px*px_nm,'y_height_nm':h_px*px_nm,'tip_gap_y_nm':np.nan,'x_width_px':w_px,'y_height_px':h_px,
            'center_image_x_px':float(d.center[0]),'center_image_y_px':float(d.center[1]),'center_axis_x_px':float((d.center-origin)@xa),'center_axis_y_up_px':float((d.center-origin)@yu),
            'slot_row_id_from_v18':4,'slot_column_id_v23':ci,'slot_pixel_area_px':contour_pixel_count(d.contour),'slot_regularity_solidity':reg['regularity_solidity'],
            'slot_regularity_convexity':reg['regularity_convexity'],'slot_best_template_iou':reg['regularity_template_iou'],'selection_source':_source_of(d),'pixel_size_nm':px_nm,
            'x_axis_angle_from_image_x_deg':axis['x_axis_angle_from_image_x_deg'],'y_axis_down_angle_from_image_x_deg':axis['y_axis_down_angle_from_image_x_deg'],
            'axis_pair_count':axis['pair_count'],'axis_fallback':axis['fallback'],'slot_complete_column_count':len(cols),
        })
    # Everything V18 selected but not belonging to a complete 4-row column is rejected.
    for d,old in paired:
        if id(d) not in used_ids:
            rejects.append(_make_reject_from_contour(d.contour,'slot_v23_column_qc','slot_not_in_complete_4row_column','candidate is not part of a complete column containing rows 1,2,3,4',desc=d,slot_row=int(old.get('slot_row',0))))
    meta=dict(meta); meta.update({'v23_complete_slot_columns':len(cols),'v23_slot_grid_required':True,'v23_slot_column_pitch_px':pitch,'v23_slot_column_tol_px':tol,
        'v23_axis_pair_count':axis['pair_count'],'v23_axis_fallback':axis['fallback'],'v23_x_axis_angle_from_image_x_deg':axis['x_axis_angle_from_image_x_deg'],
        'v23_y_axis_down_angle_from_image_x_deg':axis['y_axis_down_angle_from_image_x_deg'],'v23_slot_axis_definition':'TLS through all 4 centers of every complete column'})
    return rows,selected_all,bottom,rejects,meta


def _annotate_slots_v19(gray: np.ndarray, selected_all: List[ShapeDesc], bottom: List[ShapeDesc], rows: List[dict]) -> np.ndarray:
    out=to_color(gray); bottom_ids={id(d) for d in bottom}; row_by_col={int(r['object_id']):r for r in rows}
    for d in selected_all:
        if id(d) not in bottom_ids:
            cv2.drawContours(out,[np.rint(d.contour).astype(np.int32)],-1,(255,180,0),1)
            c=tuple(np.round(d.center).astype(int)); ci=int(getattr(d,'_v23_slot_col',0)); rid=int(getattr(d,'_v23_slot_row',0))
            cv2.putText(out,f'{ci}.{rid}',(c[0]+2,c[1]-2),cv2.FONT_HERSHEY_SIMPLEX,0.27,(255,180,0),1,cv2.LINE_AA)
    for d in bottom:
        ci=int(getattr(d,'_v23_slot_col',0)); r=row_by_col.get(ci)
        cv2.drawContours(out,[np.rint(d.contour).astype(np.int32)],-1,(0,220,0),2); c=tuple(np.round(d.center).astype(int))
        if r is not None:
            _label_box_v23(out,(c[0]-10,c[1]-7),f"S{ci}",(0,255,0),0.30)
            _label_box_v23(out,(c[0]-22,c[1]+13),f"{r['x_width_nm']:.1f}x{r['y_height_nm']:.1f}",(0,255,0),0.28)
    if selected_all and hasattr(selected_all[0],'_v22_slot_axis'): _draw_axes_v22(out,getattr(selected_all[0],'_v22_slot_axis'))
    if rows:
        mw=float(np.mean([r['x_width_nm'] for r in rows])); mh=float(np.mean([r['y_height_nm'] for r in rows])); text=f"complete columns={len(rows)} bottom row Xmean={mw:.2f}nm Ymean={mh:.2f}nm"
    else: text='slot: no complete 4-row column'
    cv2.rectangle(out,(4,4),(min(gray.shape[1]-4,650),29),(0,0,0),-1); cv2.putText(out,text,(8,22),cv2.FONT_HERSHEY_SIMPLEX,0.45,(255,255,255),1,cv2.LINE_AA)
    return out




# -----------------------------
# V23 PER-OBJECT BEFORE/AFTER MATCHING
# -----------------------------
def _normalized_object_axis_x_v23(g: pd.DataFrame) -> np.ndarray:
    if g.empty:
        return np.asarray([], dtype=float)
    x = pd.to_numeric(g.get('center_axis_x_px', pd.Series(index=g.index,dtype=float)), errors='coerce').to_numpy(dtype=float)
    if len(x) == 1:
        return np.array([0.5], dtype=float)
    if not np.all(np.isfinite(x)) or float(np.nanmax(x)-np.nanmin(x)) < 1e-9:
        return np.linspace(0.0, 1.0, len(x))
    return (x - float(np.min(x))) / max(float(np.max(x)-np.min(x)), 1e-12)


def _monotonic_object_pairs_v23(gc: pd.DataFrame) -> List[Tuple[Optional[pd.Series], Optional[pd.Series]]]:
    """Align left-to-right objects across before/after without index shifting after a miss."""
    stages={}
    for st in STAGE_ORDER_V19:
        q=gc[gc['stage']==st].copy()
        if q.empty:
            stages[st]=q; continue
        if 'center_axis_x_px' in q.columns and pd.to_numeric(q['center_axis_x_px'],errors='coerce').notna().any():
            q['_sx']=pd.to_numeric(q['center_axis_x_px'],errors='coerce')
        else:
            q['_sx']=pd.to_numeric(q['object_id'],errors='coerce')
        q=q.sort_values(['_sx','object_id']).reset_index(drop=True)
        stages[st]=q
    b=stages.get('before',pd.DataFrame()); a=stages.get('after',pd.DataFrame())
    if b.empty:
        return [(None,a.iloc[j]) for j in range(len(a))]
    if a.empty:
        return [(b.iloc[i],None) for i in range(len(b))]
    ub=_normalized_object_axis_x_v23(b); ua=_normalized_object_axis_x_v23(a)
    nb,na=len(b),len(a); gap=0.22
    dp=np.full((nb+1,na+1),np.inf,float); bt=np.full((nb+1,na+1),'',object); dp[0,0]=0.0
    for i in range(1,nb+1): dp[i,0]=dp[i-1,0]+gap; bt[i,0]='B'
    for j in range(1,na+1): dp[0,j]=dp[0,j-1]+gap; bt[0,j]='A'
    for i in range(1,nb+1):
        for j in range(1,na+1):
            d=abs(float(ub[i-1]-ua[j-1]))
            match_cost=d if d<=0.34 else 0.60+d
            opts=[(dp[i-1,j-1]+match_cost,'M'),(dp[i-1,j]+gap,'B'),(dp[i,j-1]+gap,'A')]
            val,op=min(opts,key=lambda z:z[0]); dp[i,j]=val; bt[i,j]=op
    rev=[]; i,j=nb,na
    while i>0 or j>0:
        op=bt[i,j]
        if op=='M': rev.append((b.iloc[i-1],a.iloc[j-1])); i-=1; j-=1
        elif op=='B': rev.append((b.iloc[i-1],None)); i-=1
        elif op=='A': rev.append((None,a.iloc[j-1])); j-=1
        else:
            if i>0: rev.append((b.iloc[i-1],None)); i-=1
            elif j>0: rev.append((None,a.iloc[j-1])); j-=1
    return list(reversed(rev))


def _object_before_after_wide_v22(objects_df: pd.DataFrame) -> pd.DataFrame:
    """V23: monotonic geometric correspondence; a missing object leaves a single bar."""
    if objects_df.empty:
        return pd.DataFrame()
    records=[]
    for pattern,metrics in V22_PER_OBJECT_METRICS.items():
        gp=objects_df[objects_df['pattern']==pattern].copy()
        if gp.empty: continue
        for condition in CONDITION_IDS_V19:
            gc=gp[gp['condition'].astype(str)==str(condition)].copy()
            if gc.empty: continue
            pairs=_monotonic_object_pairs_v23(gc)
            prefix='T' if pattern=='trench160' else 'S'
            for mid,(br,ar) in enumerate(pairs,1):
                for metric,metric_zh in metrics:
                    rec={'condition':str(condition),'condition_label':CONDITION_LABELS_V19.get(str(condition),f'Condition {condition}'),
                         'pattern':pattern,'pattern_zh':PATTERN_LABELS_ZH_V19.get(pattern,pattern),'object_id':mid,'object_label':f'{prefix}{mid}',
                         'metric':metric,'metric_zh':metric_zh,
                         'before_local_object_id':int(br['object_id']) if br is not None and pd.notna(br.get('object_id',np.nan)) else np.nan,
                         'after_local_object_id':int(ar['object_id']) if ar is not None and pd.notna(ar.get('object_id',np.nan)) else np.nan}
                    for st,row in [('before',br),('after',ar)]:
                        if row is None:
                            rec[f'{st}_value_nm']=np.nan; rec[f'{st}_image']=''; rec[f'{st}_center_axis_x_px']=np.nan
                        else:
                            v=pd.to_numeric(pd.Series([row.get(metric,np.nan)]),errors='coerce').iloc[0]
                            rec[f'{st}_value_nm']=float(v) if np.isfinite(v) else np.nan
                            rec[f'{st}_image']=str(row.get('image','')); rec[f'{st}_center_axis_x_px']=float(row.get('center_axis_x_px',np.nan))
                    bv=rec['before_value_nm']; av=rec['after_value_nm']
                    rec['delta_after_minus_before_nm']=float(av-bv) if np.isfinite(av) and np.isfinite(bv) else np.nan
                    rec['change_percent']=float(100.0*(av-bv)/bv) if np.isfinite(av) and np.isfinite(bv) and abs(bv)>1e-12 else np.nan
                    records.append(rec)
    if not records: return pd.DataFrame()
    out=pd.DataFrame(records); order={s:i for i,s in enumerate(CONDITION_IDS_V19)}; po={'trench160':0,'slot210':1}
    out['_co']=out['condition'].map(order).fillna(999); out['_po']=out['pattern'].map(po).fillna(999)
    return out.sort_values(['_co','_po','metric','object_id']).drop(columns=['_co','_po']).reset_index(drop=True)


def _save_per_object_bar_v23(
    objects_df: pd.DataFrame, condition: str, pattern: str, metric: str, metric_zh: str, out_path: Path
) -> None:
    sub=objects_df[(objects_df['condition'].astype(str)==str(condition))&(objects_df['pattern']==pattern)].copy()
    if sub.empty or metric not in sub.columns: return
    pairs=_monotonic_object_pairs_v23(sub)
    if not pairs: return
    _configure_chinese_matplotlib(); fig_w=max(9.0,1.05*len(pairs)+4.0); fig,ax=plt.subplots(figsize=(fig_w,6.0))
    x=np.arange(len(pairs),dtype=float); width=0.36; prefix='T' if pattern=='trench160' else 'S'; all_values=[]
    for jj,st in enumerate(STAGE_ORDER_V19):
        vals=[]
        for br,ar in pairs:
            row=br if st=='before' else ar
            if row is None: vals.append(np.nan); continue
            v=pd.to_numeric(pd.Series([row.get(metric,np.nan)]),errors='coerce').iloc[0]
            vals.append(float(v) if np.isfinite(v) else np.nan)
        all_values.extend([v for v in vals if np.isfinite(v)]); xpos=x+(jj-0.5)*width; first=True
        for xx,v in zip(xpos,vals):
            if not np.isfinite(v): continue
            bar=ax.bar([xx],[v],width,label=STAGE_LABELS_ZH_V19[st] if first else None)[0]; first=False
            ax.annotate(f'{v:.2f}',xy=(bar.get_x()+bar.get_width()/2,bar.get_height()),xytext=(0,4),textcoords='offset points',ha='center',va='bottom',fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels([f'{prefix}{i+1}' for i in range(len(pairs))]); ax.set_xlabel('从左到右匹配后的图案序号'); ax.set_ylabel('尺寸 / nm')
    ax.set_title(f"Condition {condition}  {PATTERN_LABELS_ZH_V19[pattern]}：逐图案{metric_zh}处理前后对比")
    handles,labels=ax.get_legend_handles_labels();
    if handles: ax.legend()
    ax.grid(axis='y',alpha=0.25); fig.tight_layout(); out_path.parent.mkdir(parents=True,exist_ok=True); fig.savefig(out_path,dpi=PLOT_DPI_V19,bbox_inches='tight'); plt.close(fig)


def _make_per_object_plots_v22(objects_df: pd.DataFrame, plot_dir: Path) -> None:
    if objects_df.empty: return
    root=plot_dir/'per_object'
    for condition in CONDITION_IDS_V19:
        for pattern,metrics in V22_PER_OBJECT_METRICS.items():
            for metric,metric_zh in metrics:
                _save_per_object_bar_v23(objects_df,condition,pattern,metric,metric_zh,
                    root/f'Condition_{condition}'/f'{pattern}__{metric}__objects_before_after.png')

def _expected_inventory_status_v19(
    stage: str,
    condition: str,
    images: List[Path],
) -> List[dict]:
    counts = Counter(_classify_pattern_key_v19(p) for p in images)
    rows = []
    for pattern in PATTERN_KEYS_V19:
        count = int(counts.get(pattern, 0))
        if count == 1:
            continue
        rows.append({
            "stage": stage,
            "stage_zh": STAGE_LABELS_ZH_V19[stage],
            "condition": condition,
            "condition_label": CONDITION_LABELS_V19.get(str(condition), f"Condition {condition}"),
            "image": "",
            "path": "",
            "pattern": pattern,
            "pattern_zh": PATTERN_LABELS_ZH_V19[pattern],
            "status": "MISSING_EXPECTED_IMAGE" if count == 0 else "MULTIPLE_IMAGES_WITH_SAME_PREFIX",
            "note": f"Expected exactly one image for {pattern}; found {count}",
        })
    return rows


def main_v19(
    before_root: Path,
    after_root: Path,
    output_root: Optional[Path] = None,
) -> None:
    before_root = Path(before_root)
    after_root = Path(after_root)
    if not before_root.exists():
        raise FileNotFoundError(f"Before root does not exist: {before_root}")
    if not after_root.exists():
        raise FileNotFoundError(f"After root does not exist: {after_root}")
    if output_root is None:
        output_root = before_root.parent / OUTPUT_DIR_NAME_V19
    output_root = Path(output_root)
    ann_root = output_root / "annotated"
    plot_dir = output_root / "plots"
    output_root.mkdir(parents=True, exist_ok=True)

    print("=" * 90)
    print("SEM BEFORE/AFTER V25: strict first pass + automatic zero-result fallback; structural QC retained")
    print(f"Script          : {Path(__file__).resolve()}")
    print(f"Before root     : {before_root.resolve()}")
    print(f"After root      : {after_root.resolve()}")
    print(f"Output root     : {output_root.resolve()}")
    print(f"Script version  : {SCRIPT_VERSION}")
    print("=" * 90)

    roots = {"before": before_root, "after": after_root}
    all_rows: List[dict] = []
    status_rows: List[dict] = []
    rejected_rows: List[dict] = []

    for stage in STAGE_ORDER_V19:
        root = roots[stage]
        print(f"\n[{STAGE_LABELS_ZH_V19[stage]}] {root}")
        for condition in CONDITION_IDS_V19:
            condition_dir = root / condition
            ann_dir = ann_root / stage / condition
            ann_dir.mkdir(parents=True, exist_ok=True)
            if not condition_dir.exists():
                status_rows.append({
                    "stage": stage,
                    "stage_zh": STAGE_LABELS_ZH_V19[stage],
                    "condition": condition,
                    "condition_label": CONDITION_LABELS_V19.get(str(condition), f"Condition {condition}"),
                    "image": "",
                    "path": str(condition_dir),
                    "pattern": "",
                    "status": "MISSING_CONDITION_FOLDER",
                    "note": f"Condition folder not found: {condition_dir}",
                })
                print(f"  [Condition {condition}] missing folder: {condition_dir}")
                continue

            images = _find_sem_images_v19(condition_dir)
            print(f"  [Condition {condition}] {len(images)} image(s)")
            status_rows.extend(_expected_inventory_status_v19(stage, condition, images))
            for image_path in images:
                try:
                    rows, status, rejected = _process_one_image_v19(
                        image_path, stage, condition, ann_dir
                    )
                    all_rows.extend(rows)
                    status_rows.append(status)
                    rejected_rows.extend(rejected)
                    print(
                        f"      selected measurements={len(rows)}; status={status['status']}; "
                        f"PixelSize={status.get('pixel_size_nm', np.nan):g} nm/px"
                    )
                except PixelSizeTxtNotFound as e:
                    status_rows.append({
                        "stage": stage,
                        "stage_zh": STAGE_LABELS_ZH_V19[stage],
                        "condition": condition,
                        "condition_label": CONDITION_LABELS_V19.get(str(condition), f"Condition {condition}"),
                        "image": image_path.name,
                        "path": str(image_path),
                        "pattern": _classify_pattern_key_v19(image_path) or "unknown",
                        "status": "ERROR_TXT_NOT_FOUND",
                        "note": str(e),
                    })
                    print(f"      [PixelSize TXT missing] {e}")
                except PixelSizeFieldNotFound as e:
                    status_rows.append({
                        "stage": stage,
                        "stage_zh": STAGE_LABELS_ZH_V19[stage],
                        "condition": condition,
                        "condition_label": CONDITION_LABELS_V19.get(str(condition), f"Condition {condition}"),
                        "image": image_path.name,
                        "path": str(image_path),
                        "pattern": _classify_pattern_key_v19(image_path) or "unknown",
                        "status": "ERROR_PIXELSIZE_FIELD_NOT_FOUND",
                        "note": str(e),
                    })
                    print(f"      [PixelSize field missing] {e}")
                except Exception as e:
                    status_rows.append({
                        "stage": stage,
                        "stage_zh": STAGE_LABELS_ZH_V19[stage],
                        "condition": condition,
                        "condition_label": CONDITION_LABELS_V19.get(str(condition), f"Condition {condition}"),
                        "image": image_path.name,
                        "path": str(image_path),
                        "pattern": _classify_pattern_key_v19(image_path) or "unknown",
                        "status": "ERROR_OTHER",
                        "note": f"{type(e).__name__}: {e}",
                    })
                    print(f"      [ERROR] {type(e).__name__}: {e}")

    objects_df = pd.DataFrame(all_rows)
    status_df = pd.DataFrame(status_rows)
    rejected_df = pd.DataFrame(rejected_rows)

    image_metric_df = _metric_summary_v19(
        objects_df,
        ["stage", "stage_zh", "condition", "condition_label", "image", "pattern", "pattern_zh"],
    )
    condition_metric_df = _metric_summary_v19(
        objects_df,
        ["stage", "stage_zh", "condition", "condition_label", "pattern", "pattern_zh"],
    )
    before_after_df = _before_after_wide_v19(condition_metric_df)
    object_before_after_df = _object_before_after_wide_v22(objects_df)

    objects_df.to_csv(output_root / "all_object_measurements.csv", index=False, encoding="utf-8-sig")
    status_df.to_csv(output_root / "image_status.csv", index=False, encoding="utf-8-sig")
    rejected_df.to_csv(output_root / "rejected_objects.csv", index=False, encoding="utf-8-sig")
    image_metric_df.to_csv(output_root / "image_metric_summary.csv", index=False, encoding="utf-8-sig")
    condition_metric_df.to_csv(output_root / "condition_metric_summary.csv", index=False, encoding="utf-8-sig")
    before_after_df.to_csv(output_root / "before_after_comparison.csv", index=False, encoding="utf-8-sig")
    object_before_after_df.to_csv(output_root / "object_before_after_comparison.csv", index=False, encoding="utf-8-sig")

    xlsx_path = output_root / "sem_before_after_results.xlsx"
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        status_df.to_excel(writer, sheet_name="image_status", index=False)
        objects_df.to_excel(writer, sheet_name="all_objects", index=False)
        rejected_df.to_excel(writer, sheet_name="rejected_objects", index=False)
        image_metric_df.to_excel(writer, sheet_name="image_metric", index=False)
        condition_metric_df.to_excel(writer, sheet_name="condition_metric", index=False)
        before_after_df.to_excel(writer, sheet_name="before_after", index=False)
        object_before_after_df.to_excel(writer, sheet_name="object_compare", index=False)
        if not objects_df.empty:
            for pattern in PATTERN_KEYS_V19:
                sub = objects_df[objects_df["pattern"] == pattern]
                if not sub.empty:
                    sub.to_excel(writer, sheet_name=pattern[:31], index=False)
        _autofit_excel_v19(writer)

    if not condition_metric_df.empty:
        _make_all_plots_v19(condition_metric_df, plot_dir)
    if not objects_df.empty:
        _make_per_object_plots_v22(objects_df, plot_dir)

    print("\nDone.")
    print(f"Output folder     : {output_root}")
    print(f"Excel             : {xlsx_path}")
    print(f"Annotated images  : {ann_root}")
    print(f"Bar charts        : {plot_dir} (including per_object/Condition_1..8)")
    print("Annotation legend : GREEN=used; YELLOW=via bottom row; BLUE=other valid slot rows; MAGENTA trench points define Y axis. Labels show per-object X/Y; trench uses W=mean upper/lower rectangle width and Y=Y-axis tip gap. Rejected overlays are off by default.")



# ============================================================
# V24 FINAL VIA OVERRIDES
#   * via40-after: lower centre-evidence gate only at predicted lattice sites;
#     mild local smoothing recovers blurry holes without accepting free-floating noise.
#   * via40/via60 boundary: follow the OUTERMOST sustained dark-to-bright transition,
#     then merge it with the connected local dark component.  This prevents an inner
#     dark-core edge from shrinking the measured hole.
#   * shape/grid/physical-size QC remains mandatory.
# ============================================================
V24_VIA40_AFTER_RELAXED_PEAK_SCORE = 0.19
V24_VIA40_AFTER_RELAXED_RAW_FACTOR = 0.62
V24_VIA40_AFTER_SITE_SCORE_FACTOR = 0.70
V24_VIA40_AFTER_SITE_RAW_FACTOR = 0.68
V24_VIA40_AFTER_SUPPORT_SCORE_FACTOR = 0.42
V24_VIA40_AFTER_SUPPORT_RAW_FACTOR = 0.28
V24_VIA_RADIAL_ANGLES = 160
V24_VIA_MAX_RADIUS_FACTOR = 1.58
V24_VIA_MIN_EDGE_RADIUS_FACTOR = 0.56
V24_VIA40_AFTER_DARKNESS_THRESHOLD = 0.115
V24_VIA60_DARKNESS_THRESHOLD = 0.115
V24_VIA_OTHER_DARKNESS_THRESHOLD = 0.150
V24_VIA40_AFTER_MIN_GOOD_RAYS = 0.50
V24_VIA60_MIN_GOOD_RAYS = 0.56
V24_VIA_OTHER_MIN_GOOD_RAYS = 0.62
V24_VIA_MAX_OUTER_DARK_LEAK = 0.34
V24_VIA_COMPONENT_CLOSE_DIAM_FRAC = 0.075
V24_VIA_COMPONENT_MIN_CORE_OVERLAP_FRAC = 0.12
V24_VIA_RADIUS_SPIKE_LIMIT_R = 0.24


def _circular_median_v24(values: np.ndarray, half_window: int = 2) -> np.ndarray:
    a = np.asarray(values, dtype=np.float64)
    if len(a) < 3 or half_window <= 0:
        return a.copy()
    stack = [np.roll(a, k) for k in range(-half_window, half_window + 1)]
    return np.nanmedian(np.vstack(stack), axis=0)


def _smooth_profile_v24(profile: np.ndarray, sigma_samples: float = 1.15) -> np.ndarray:
    p = np.asarray(profile, dtype=np.float32).reshape(1, -1)
    if p.shape[1] < 5:
        return p.ravel().astype(np.float64)
    return cv2.GaussianBlur(p, (0, 0), sigmaX=max(0.35, float(sigma_samples))).ravel().astype(np.float64)


def _via_processing_image_v24(
    gray: np.ndarray, nominal_diam_px: float, stage: str, pattern_key: str
) -> Tuple[np.ndarray, float]:
    """Small-object smoothing only; preserve the physical outer edge."""
    D = max(7.0, float(nominal_diam_px))
    src = gray.astype(np.float32)
    if stage == 'after' and pattern_key == 'via40':
        # A 3x3 median removes isolated SEM spikes; Gaussian smoothing then joins a
        # slightly broken/blurred etched edge.  It is deliberately not a large blur.
        med = cv2.medianBlur(np.clip(src, 0, 255).astype(np.uint8), 3).astype(np.float32)
        sigma = max(0.95, min(1.75, 0.034 * D))
        sm = cv2.GaussianBlur(med, (0, 0), sigmaX=sigma, sigmaY=sigma)
    elif pattern_key == 'via60':
        # Keep the weak outer dark shoulder of a 60-nm hole instead of emphasizing
        # only its very dark centre.
        sigma = max(0.78, min(1.45, 0.020 * D))
        sm = cv2.GaussianBlur(src, (0, 0), sigmaX=sigma, sigmaY=sigma)
    elif stage == 'after':
        sigma = max(0.72, min(1.35, 0.020 * D))
        sm = cv2.GaussianBlur(src, (0, 0), sigmaX=sigma, sigmaY=sigma)
    else:
        sigma = max(0.58, min(1.10, 0.014 * D))
        sm = cv2.GaussianBlur(src, (0, 0), sigmaX=sigma, sigmaY=sigma)
    return sm.astype(np.float32), float(sigma)


def _via_site_support_v24(
    gray: np.ndarray,
    center: np.ndarray,
    nominal_diam_px: float,
    nearest_pitch: float,
    stage: str,
    pattern_key: str,
) -> Dict[str, float]:
    """Broad dark-disk evidence at one predicted lattice site.

    This is used only to rescue weak via40-after sites.  A tiny dark speck fails because
    the shoulder (roughly 0.55R--0.96R) must also be darker than the local background.
    """
    H, W = gray.shape
    D = max(7.0, float(nominal_diam_px)); R = 0.5 * D
    c = np.asarray(center, dtype=np.float64)
    sm, sigma = _via_processing_image_v24(gray, D, stage, pattern_key)
    rlim = 1.50 * R
    if np.isfinite(nearest_pitch) and nearest_pitch > 0:
        rlim = min(rlim, 0.46 * nearest_pitch)
    half = int(math.ceil(rlim + 2.0))
    x0=max(0,int(math.floor(c[0]-half))); x1=min(W,int(math.ceil(c[0]+half+1)))
    y0=max(0,int(math.floor(c[1]-half))); y1=min(H,int(math.ceil(c[1]+half+1)))
    if x1-x0 < 7 or y1-y0 < 7:
        return {'valid':0.0,'disk_contrast_gray':0.0,'shoulder_contrast_gray':0.0,'shoulder_dark_fraction':0.0}
    yy,xx=np.mgrid[y0:y1,x0:x1]
    rr=np.sqrt((xx-c[0])**2+(yy-c[1])**2)
    loc=sm[y0:y1,x0:x1]
    core=rr<=0.40*R
    shoulder=(rr>=0.56*R)&(rr<=0.96*R)
    outer=(rr>=1.18*R)&(rr<=rlim)
    if min(np.count_nonzero(core),np.count_nonzero(shoulder),np.count_nonzero(outer)) < 8:
        return {'valid':0.0,'disk_contrast_gray':0.0,'shoulder_contrast_gray':0.0,'shoulder_dark_fraction':0.0}
    # Use a bright-side percentile outside.  If part of the outer annulus is still
    # occupied by the hole, this remains a better background estimate than its median.
    bg=float(np.percentile(loc[outer],72.0))
    core_g=float(np.median(loc[core])); shoulder_g=float(np.median(loc[shoulder]))
    contrast=bg-core_g; sh_contrast=bg-shoulder_g
    global_range=max(1.0,float(np.percentile(sm,95)-np.percentile(sm,5)))
    contrast_min=max(0.16,0.0035*global_range)
    sh_min=max(0.07,0.025*max(contrast,0.0))
    qthr = 0.10 if (stage=='after' and pattern_key=='via40') else 0.12
    dark_cut=bg-qthr*max(contrast,1e-6)
    sh_frac=float(np.mean(loc[shoulder] <= dark_cut))
    frac_min=0.23 if (stage=='after' and pattern_key=='via40') else 0.30
    valid=bool(contrast>=contrast_min and sh_contrast>=sh_min and sh_frac>=frac_min)
    return {
        'valid':1.0 if valid else 0.0,
        'disk_contrast_gray':float(contrast),
        'shoulder_contrast_gray':float(sh_contrast),
        'shoulder_dark_fraction':float(sh_frac),
        'support_smooth_sigma_px':float(sigma),
        'local_background_gray':float(bg),
        'local_core_gray':float(core_g),
    }


def _via_relaxed_scale_peaks_v24(
    score: np.ndarray,
    raw: np.ndarray,
    best_diam: np.ndarray,
    nominal_diam_px: float,
    raw_floor: float,
) -> List[ShapeDesc]:
    """Extra centre proposals for low-contrast via40-after images.

    They are proposals only.  They still need a reliable rectangular lattice, broad
    disk support, physical size, and final shape QC.
    """
    H,W=score.shape
    k=int(max(3,round(0.48*nominal_diam_px))); k += (k % 2 == 0)
    ker=cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(k,k))
    dil=cv2.dilate(score,ker)
    maxima=(score>=dil-1e-6)&(score>=V24_VIA40_AFTER_RELAXED_PEAK_SCORE)&(raw>=V24_VIA40_AFTER_RELAXED_RAW_FACTOR*raw_floor)
    margin=max(3,int(math.ceil(V21_VIA_BORDER_DIAM_FRAC*nominal_diam_px)))
    maxima[:margin,:]=False; maxima[-margin:,:]=False; maxima[:,:margin]=False; maxima[:,-margin:]=False
    peaks=_peak_components_v21(maxima,score)
    out=[]
    for sc,x,y in peaks[:V21_VIA_MAX_PEAKS]:
        dd=float(best_diam[y,x]) if np.isfinite(best_diam[y,x]) else nominal_diam_px
        d=_make_circle_desc_v21(np.array([x,y],float),dd,'v24_relaxed_via40_after_peak')
        if d is None: continue
        setattr(d,'_v21_score',float(sc)); setattr(d,'_v21_raw',float(raw[y,x])); setattr(d,'_v21_detect_diam_px',dd)
        out.append(d)
    return out


def _via_component_radii_v24(
    sm: np.ndarray,
    center: np.ndarray,
    R: float,
    rmax: float,
    bg: float,
    contrast: float,
    darkness_threshold: float,
    angles: np.ndarray,
    rr_sample: np.ndarray,
    stage: str,
    pattern_key: str,
) -> Tuple[np.ndarray, Dict[str,float]]:
    """Radii of the centre-connected local dark component at a permissive threshold."""
    H,W=sm.shape; c=np.asarray(center,float)
    half=int(math.ceil(rmax+3.0))
    x0=max(0,int(math.floor(c[0]-half))); x1=min(W,int(math.ceil(c[0]+half+1)))
    y0=max(0,int(math.floor(c[1]-half))); y1=min(H,int(math.ceil(c[1]+half+1)))
    yy,xx=np.mgrid[y0:y1,x0:x1]
    radial=np.sqrt((xx-c[0])**2+(yy-c[1])**2)
    loc=sm[y0:y1,x0:x1]
    darkness=np.clip((bg-loc)/max(contrast,1e-6),-0.5,2.0)
    # Threshold is intentionally near the background side so a weak but visibly dark
    # outer shoulder is included.  Grid/shape/size QC prevents arbitrary expansion.
    mask=((darkness>=darkness_threshold)&(radial<=rmax)).astype(np.uint8)
    k=max(3,int(round(V24_VIA_COMPONENT_CLOSE_DIAM_FRAC*(2.0*R))))
    if k%2==0: k+=1
    ker=cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(k,k))
    mask=cv2.morphologyEx(mask,cv2.MORPH_CLOSE,ker,iterations=1)
    # Remove isolated one-pixel dark noise, but do not erode a blurry via40-after edge.
    if not (stage=='after' and pattern_key=='via40'):
        mask=cv2.morphologyEx(mask,cv2.MORPH_OPEN,cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(3,3)),iterations=1)
    n,lab=cv2.connectedComponents(mask,8)
    if n<=1:
        return np.full(len(angles),np.nan,float),{'component_used':0.0,'component_core_overlap':0.0,'component_bound_fraction':0.0}
    core=radial<=0.42*R
    best_lab=0; best_score=-1.0; best_overlap=0.0
    core_n=max(1,int(np.count_nonzero(core)))
    for q in range(1,n):
        comp=lab==q
        overlap=float(np.count_nonzero(comp&core))/core_n
        area=float(np.count_nonzero(comp))
        # Prefer the component that occupies the expected centre disk, not a nearby speck.
        s=4.0*overlap+min(area/(math.pi*R*R+1e-6),1.5)
        if s>best_score:
            best_score=s; best_lab=q; best_overlap=overlap
    if best_lab<=0 or best_overlap<V24_VIA_COMPONENT_MIN_CORE_OVERLAP_FRAC:
        return np.full(len(angles),np.nan,float),{'component_used':0.0,'component_core_overlap':best_overlap,'component_bound_fraction':0.0}
    comp=(lab==best_lab).astype(np.float32)
    xs=(c[0]+np.outer(np.cos(angles),rr_sample)).astype(np.float32)
    ys=(c[1]+np.outer(np.sin(angles),rr_sample)).astype(np.float32)
    # Convert full-image coordinates to local-ROI coordinates for remap.
    vals=cv2.remap(comp,(xs-x0).astype(np.float32),(ys-y0).astype(np.float32),cv2.INTER_NEAREST,borderMode=cv2.BORDER_CONSTANT,borderValue=0)
    cr=np.full(len(angles),np.nan,float)
    hit_bound=0
    for i,row in enumerate(vals):
        idx=np.where(row>0.5)[0]
        if len(idx):
            j=int(idx[-1]); cr[i]=float(rr_sample[j])
            if cr[i]>=rmax-0.6: hit_bound+=1
    return cr,{
        'component_used':1.0,
        'component_core_overlap':float(best_overlap),
        'component_bound_fraction':float(hit_bound/max(1,len(angles))),
        'component_area_px':float(np.count_nonzero(comp)),
    }


def _via_outer_edge_contour_v24(
    gray: np.ndarray,
    center: np.ndarray,
    nominal_diam_px: float,
    nearest_pitch: float,
    stage: str,
    pattern_key: str,
) -> Tuple[Optional[np.ndarray], Dict[str,float]]:
    """Trace the complete physical via edge, choosing the outermost sustained transition.

    The old strongest-gradient rule could stop at an internal dark-core/outer-shoulder
    transition.  Here each ray is normalized by its local core/background contrast;
    the outermost sustained dark-to-bright crossing is selected and then enlarged to
    cover the centre-connected dark component when that component proves the edge lies
    farther out.
    """
    H,W=gray.shape; D=max(7.0,float(nominal_diam_px)); R=0.5*D; c=np.asarray(center,float)
    sm,sigma=_via_processing_image_v24(gray,D,stage,pattern_key)
    rmax=V24_VIA_MAX_RADIUS_FACTOR*R
    if np.isfinite(nearest_pitch) and nearest_pitch>0:
        rmax=min(rmax,0.46*nearest_pitch)
    rmin=max(1.0,0.08*R); rmin_edge=V24_VIA_MIN_EDGE_RADIUS_FACTOR*R
    if rmax<=rmin_edge+2.0:
        return None,{'good_ray_fraction':0.0,'smooth_sigma_px':sigma}

    # Common robust core/background estimate for the 2-D component.  Per-ray estimates
    # below adapt to local directional shading.
    half=int(math.ceil(rmax+3)); x0=max(0,int(c[0]-half)); x1=min(W,int(c[0]+half+1)); y0=max(0,int(c[1]-half)); y1=min(H,int(c[1]+half+1))
    yy,xx=np.mgrid[y0:y1,x0:x1]; radial=np.sqrt((xx-c[0])**2+(yy-c[1])**2); loc=sm[y0:y1,x0:x1]
    core=radial<=0.40*R; outer=(radial>=max(1.20*R,0.78*rmax))&(radial<=rmax)
    shoulder=(radial>=0.55*R)&(radial<=0.96*R)
    if min(np.count_nonzero(core),np.count_nonzero(outer),np.count_nonzero(shoulder))<8:
        return None,{'good_ray_fraction':0.0,'smooth_sigma_px':sigma}
    bg_global=float(np.percentile(loc[outer],72.0)); dark_global=float(np.median(loc[core])); shoulder_g=float(np.median(loc[shoulder]))
    contrast_global=bg_global-dark_global; shoulder_contrast=bg_global-shoulder_g
    robust_range=max(1.0,float(np.percentile(sm,95)-np.percentile(sm,5)))
    cmin=max(0.15,0.0032*robust_range)
    shoulder_ratio=shoulder_contrast/max(contrast_global,1e-6)
    # A true via must have broad dark support near its expected outer half-radius.
    # This blocks tiny specks from being inflated to a nominal-size contour.
    broad_cut=bg_global-0.08*max(contrast_global,1e-6)
    shoulder_dark_fraction=float(np.mean(loc[shoulder]<=broad_cut))
    if stage=='after' and pattern_key=='via40':
        shmin=max(0.100,0.045*max(contrast_global,0.0),0.0040*robust_range); shfrac_min=0.22
        qfloor=0.060; min_good=V24_VIA40_AFTER_MIN_GOOD_RAYS
    elif pattern_key=='via60':
        shmin=max(0.200,0.085*max(contrast_global,0.0),0.0060*robust_range); shfrac_min=0.30
        qfloor=0.072; min_good=V24_VIA60_MIN_GOOD_RAYS
    else:
        shmin=max(0.200,0.085*max(contrast_global,0.0),0.0060*robust_range); shfrac_min=0.30
        qfloor=0.090; min_good=V24_VIA_OTHER_MIN_GOOD_RAYS
    if contrast_global<cmin or shoulder_contrast<shmin or shoulder_dark_fraction<shfrac_min:
        return None,{
            'good_ray_fraction':0.0,'smooth_sigma_px':sigma,'core_contrast_gray':contrast_global,
            'shoulder_contrast_gray':shoulder_contrast,'shoulder_dark_fraction':shoulder_dark_fraction,
            'outer_edge_failure':'insufficient_broad_dark_support',
        }
    # The edge level is half-way between the OUTER shoulder and local background,
    # rather than half-way between the darkest core and background.  This includes
    # the whole hole while avoiding an over-expanded low-threshold halo.
    qthr=float(np.clip(0.56*shoulder_ratio,qfloor,0.32))

    nr=max(80,int(math.ceil((rmax-rmin)*5.0)))
    rr=np.linspace(rmin,rmax,nr,dtype=np.float32)
    angles=np.linspace(0.0,2.0*math.pi,V24_VIA_RADIAL_ANGLES,endpoint=False)
    radii=np.full(len(angles),np.nan,float); strengths=np.zeros(len(angles),float); contrasts=np.zeros(len(angles),float)

    for ia,a in enumerate(angles):
        xs=(c[0]+rr*math.cos(a)).astype(np.float32); ys=(c[1]+rr*math.sin(a)).astype(np.float32)
        if xs.min()<1 or ys.min()<1 or xs.max()>W-2 or ys.max()>H-2:
            continue
        prof=cv2.remap(sm,xs.reshape(1,-1),ys.reshape(1,-1),cv2.INTER_LINEAR,borderMode=cv2.BORDER_REFLECT101).ravel().astype(np.float64)
        prof=_smooth_profile_v24(prof,1.0 if pattern_key=='via60' else 1.15)
        cm=rr<=0.42*R; bm=rr>=max(1.18*R,0.78*rmax)
        if np.count_nonzero(cm)<4 or np.count_nonzero(bm)<4: continue
        dark=float(np.median(prof[cm])); bg=float(np.percentile(prof[bm],68.0)); contrast=bg-dark
        if contrast<max(0.10,0.40*cmin): continue
        q=np.clip((bg-prof)/max(contrast,1e-6),-0.5,2.0)
        # Per-ray threshold follows the darkness of the outer shoulder on that ray.
        # It remains below the shoulder value, so an inner-core transition is ignored.
        ray_sh=(rr>=0.55*R)&(rr<=0.96*R)
        ray_sh_q=float(np.median(q[ray_sh])) if np.count_nonzero(ray_sh)>=3 else qthr*2.0
        qthr_ray=float(np.clip(0.65*(0.56*ray_sh_q)+0.35*qthr,qfloor,0.34))
        # Smooth normalized darkness; 1 means core-dark and 0 means background.
        win=max(3,int(round(0.045*D*5.0)))
        if win%2==0: win+=1
        qsm=np.convolve(q,np.ones(win,dtype=float)/win,mode='same')
        search=np.where((rr>=rmin_edge)&(rr<=rmax))[0]
        if len(search)<8: continue
        # Find the OUTERMOST crossing that has sustained dark support inside and
        # sustained background outside.  This rejects inner-core transitions.
        in_n=max(3,int(round(0.055*R*5.0))); out_n=max(3,int(round(0.070*R*5.0)))
        chosen=None; chosen_grad=0.0
        grad=np.gradient(prof,rr.astype(np.float64))
        for j in search[-2:1:-1]:
            if j+out_n>=len(qsm) or j-in_n<0: continue
            if not (qsm[j]>=qthr_ray and qsm[j+1]<qthr_ray): continue
            in_med=float(np.median(qsm[j-in_n+1:j+1])); out_med=float(np.median(qsm[j+1:j+1+out_n]))
            g=float(np.max(grad[max(0,j-1):min(len(grad),j+3)]))
            if in_med>=0.92*qthr_ray and out_med<=1.05*qthr_ray and g>=max(0.025,0.012*contrast):
                chosen=j; chosen_grad=g; break
        if chosen is None:
            # Fallback to the end of the longest centre-connected dark run. Small
            # one-dimensional gaps are closed first, so a blurry outer shoulder survives.
            b=(qsm>=qthr_ray).astype(np.uint8)
            kk=max(3,int(round(0.035*D*5.0))); kk += (kk%2==0)
            b=cv2.morphologyEx(b.reshape(1,-1),cv2.MORPH_CLOSE,np.ones((1,kk),np.uint8),iterations=1).ravel()
            start_candidates=np.where((rr<=0.48*R)&(b>0))[0]
            if len(start_candidates):
                j=int(start_candidates[-1]); last=j; gap=0
                for k in range(j+1,len(b)):
                    if b[k]>0:
                        last=k; gap=0
                    else:
                        gap+=1
                        if gap>max(2,out_n//2): break
                if rr[last]>=rmin_edge and last<len(rr)-1:
                    chosen=last; chosen_grad=float(max(0.0,grad[min(last,len(grad)-1)]))
        if chosen is None: continue
        j=int(chosen)
        # Linear threshold interpolation gives sub-sample radius.
        if j+1<len(rr) and qsm[j]!=qsm[j+1]:
            f=float(np.clip((qsm[j]-qthr_ray)/(qsm[j]-qsm[j+1]),0.0,1.0))
            rad=float(rr[j]+f*(rr[j+1]-rr[j]))
        else:
            rad=float(rr[j])
        radii[ia]=rad; strengths[ia]=chosen_grad; contrasts[ia]=contrast

    # Merge with the centre-connected 2-D dark component.  This is what recovers a
    # visibly dark outer sector that has a weaker gradient than an inner core edge.
    comp_r,comp_meta=_via_component_radii_v24(sm,c,R,rmax,bg_global,contrast_global,qthr,angles,rr,stage,pattern_key)
    comp_ok=bool(comp_meta.get('component_used',0.0)>0.5 and comp_meta.get('component_bound_fraction',1.0)<=0.16)
    if comp_ok:
        for i in range(len(radii)):
            cr=comp_r[i]
            if not np.isfinite(cr): continue
            if not np.isfinite(radii[i]):
                radii[i]=cr
            elif cr>radii[i]:
                # Expand only by a physically reasonable amount; isolated protrusions
                # are later suppressed by circular median and shape QC.
                radii[i]=min(cr,radii[i]+0.34*R)

    good=np.isfinite(radii); good_frac=float(np.mean(good)) if len(good) else 0.0
    if good_frac<min_good or np.count_nonzero(good)<8:
        return None,{
            'good_ray_fraction':good_frac,'smooth_sigma_px':sigma,'core_contrast_gray':contrast_global,
            'shoulder_contrast_gray':shoulder_contrast,**comp_meta,
        }
    idx=np.where(good)[0]
    xg=np.r_[idx-len(radii),idx,idx+len(radii)]; yg=np.r_[radii[idx],radii[idx],radii[idx]]
    miss=np.where(~good)[0]
    if len(miss): radii[miss]=np.interp(miss,xg,yg)

    # Remove only isolated radial spikes.  A genuine oval/etched distortion remains.
    med=_circular_median_v24(radii,3 if (stage=='after' and pattern_key=='via40') else 2)
    lim=V24_VIA_RADIUS_SPIKE_LIMIT_R*R
    radii=np.clip(radii,med-lim,med+lim)
    passes=3 if (stage=='after' and pattern_key=='via40') else 2
    for _ in range(passes):
        radii=0.16*np.roll(radii,2)+0.20*np.roll(radii,1)+0.28*radii+0.20*np.roll(radii,-1)+0.16*np.roll(radii,-2)

    pts=np.column_stack([c[0]+radii*np.cos(angles),c[1]+radii*np.sin(angles)])
    cont=np.asarray(np.rint(pts),np.int32).reshape(-1,1,2)

    # Quantify whether substantial dark material still lies immediately outside the
    # selected contour.  A high value signals an inward edge and is penalized when
    # choosing among refined/predicted centres.
    leak_vals=[]
    dr=max(1.0,0.075*R)
    for ia,a in enumerate(angles):
        rtest=np.array([radii[ia]+dr,radii[ia]+2*dr],np.float32)
        xs=(c[0]+rtest*math.cos(a)).astype(np.float32); ys=(c[1]+rtest*math.sin(a)).astype(np.float32)
        if xs.min()<0 or ys.min()<0 or xs.max()>W-1 or ys.max()>H-1: continue
        vv=cv2.remap(sm,xs.reshape(1,-1),ys.reshape(1,-1),cv2.INTER_LINEAR,borderMode=cv2.BORDER_REFLECT101).ravel()
        leak_vals.append(float(np.mean((bg_global-vv)/max(contrast_global,1e-6)>=0.80*qthr)))
    leak=float(np.mean(leak_vals)) if leak_vals else 1.0
    return cont,{
        'good_ray_fraction':good_frac,'median_edge_radius_px':float(np.median(radii)),
        'median_ray_contrast_gray':float(np.median(contrasts[good])) if np.any(good) else np.nan,
        'median_ray_edge_gradient':float(np.median(strengths[good])) if np.any(good) else np.nan,
        'smooth_sigma_px':float(sigma),'smooth_passes':int(passes),
        'core_contrast_gray':float(contrast_global),'shoulder_contrast_gray':float(shoulder_contrast),
        'shoulder_dark_fraction':float(shoulder_dark_fraction),'outer_dark_leak_fraction':float(leak),'darkness_threshold':float(qthr),**comp_meta,
    }


def _via_shape_qc_v24(
    contour: np.ndarray,
    nominal_diam_px: float,
    stage: str,
    pattern_key: str,
) -> Tuple[bool,Dict[str,float],List[str]]:
    m=_via_shape_metrics_v23(contour,nominal_diam_px); reasons=[]; D=max(1e-9,float(nominal_diam_px))
    if stage=='after' and pattern_key=='via40':
        circ_min,sol_min,iou_min,rcv_max,eq_min=0.43,0.76,0.48,0.31,0.70
    elif pattern_key=='via60':
        circ_min,sol_min,iou_min,rcv_max,eq_min=0.50,0.80,0.55,0.28,0.72
    else:
        circ_min,sol_min,iou_min,rcv_max,eq_min=V23_VIA_MIN_CIRCULARITY,V23_VIA_MIN_SOLIDITY,V23_VIA_MIN_ELLIPSE_IOU,V23_VIA_MAX_RADIAL_CV,V23_VIA_MIN_EQ_DIAM_FACTOR
    if m['eq_diam_px']<eq_min*D: reasons.append(f"equivalent diameter {m['eq_diam_px']:.2f}px too small")
    if m['eq_diam_px']>V23_VIA_MAX_EQ_DIAM_FACTOR*D: reasons.append(f"equivalent diameter {m['eq_diam_px']:.2f}px too large")
    if m['circularity']<circ_min: reasons.append(f"circularity {m['circularity']:.3f} < {circ_min:.3f}")
    if m['solidity']<sol_min: reasons.append(f"solidity {m['solidity']:.3f} < {sol_min:.3f}")
    if m['aspect']>V23_VIA_MAX_ASPECT: reasons.append(f"axis ratio {m['aspect']:.3f} > {V23_VIA_MAX_ASPECT:.3f}")
    if m['ellipse_iou']<iou_min: reasons.append(f"ellipse IoU {m['ellipse_iou']:.3f} < {iou_min:.3f}")
    if m['radial_cv']>rcv_max: reasons.append(f"radial CV {m['radial_cv']:.3f} > {rcv_max:.3f}")
    return len(reasons)==0,m,reasons


def _select_and_measure_vias_v19(
    gray: np.ndarray,
    pool: List[ShapeDesc],
    px_nm: float,
    image_name: str,
    stage: str,
    condition: str,
    pattern_key: str,
    search_diam_px: Optional[float] = None,
) -> Tuple[List[dict], List[ShapeDesc], List[dict], Dict[str, object]]:
    nominal_nm,nominal_D=_via_nominal_diam_px_v21(pattern_key,px_nm)
    D = nominal_D if search_diam_px is None else max(6.0, float(search_diam_px))
    score,raw,best_diam,resp_meta=_via_scale_response_v21(gray,D); raw_floor=float(resp_meta['raw_contrast_floor'])
    scale_peaks=_via_scale_peaks_v21(gray,score,raw,best_diam,D,raw_floor)
    legacy=_via_legacy_size_candidates_v21(pool,score,raw,D,raw_floor,gray.shape)
    candidates=_merge_via_center_candidates_v21(scale_peaks+legacy,D)
    model,model_seeds=_choose_via_lattice_v21(candidates,D)

    # Low-contrast via40-after may not provide enough standard peaks to establish the
    # grid.  Add lower-threshold centre proposals, but only adopt them if the resulting
    # model still passes the same rectangular-lattice validation.
    if stage=='after' and pattern_key=='via40' and not _validate_via_model_v23(model,D):
        relaxed=_via_relaxed_scale_peaks_v24(score,raw,best_diam,D,raw_floor)
        c2=_merge_via_center_candidates_v21(candidates+relaxed,D)
        m2,s2=_choose_via_lattice_v21(c2,D)
        if _validate_via_model_v23(m2,D):
            candidates,model,model_seeds=c2,m2,s2

    if not _validate_via_model_v23(model,D):
        strong=[d for d in candidates if float(getattr(d,'_v21_score',0.0))>=0.30 and float(getattr(d,'_v21_raw',0.0))>=0.82*raw_floor]
        m2,s2=_choose_via_lattice_v21(strong,D)
        if _validate_via_model_v23(m2,D): model,model_seeds=m2,s2
    rejects=[]
    if not _validate_via_model_v23(model,D):
        for d in candidates[:60]:
            if float(getattr(d,'_v21_score',0.0))>=0.28:
                rejects.append(_make_reject_from_contour(d.contour,'via_v24_grid_qc','via_no_reliable_rectangular_lattice',
                    'candidate may be dark/round, but no reliable rectangular via lattice was established',desc=d,
                    scale_match_score=float(getattr(d,'_v21_score',np.nan))))
        axis=_fit_via_bottom_row_axes_v19([],gray.shape)
        meta={'via_v24_final_count':0,'via_v24_grid_required':True,'via_v24_lattice_used':False,**resp_meta}
        return [],[],rejects,{'axis':axis,'meta':meta,'lattice_model':model}

    lens=[float(model.get('basis_a_len',np.nan)),float(model.get('basis_b_len',np.nan))]; lens=[x for x in lens if np.isfinite(x) and x>0]
    nearest_pitch=min(lens) if lens else np.nan
    margin=max(3.0,V21_VIA_BORDER_DIAM_FRAC*D); sites=lattice_sites_v14(model,gray.shape,1,'via')
    pre=[]; recovered=0
    for site in sites:
        xy0=np.array([float(site['x']),float(site['y'])],float)
        if xy0[0]<=margin or xy0[1]<=margin or xy0[0]>=gray.shape[1]-1-margin or xy0[1]>=gray.shape[0]-1-margin: continue
        xy,sc,rw=_via_refine_site_center_v21(score,raw,xy0,D,nearest_pitch); internal=bool(site.get('is_internal',False))
        sc_min=V21_VIA_INTERNAL_SITE_SCORE_MIN if internal else V21_VIA_EXTERNAL_SITE_SCORE_MIN
        rw_min=raw_floor if internal else V21_VIA_EXTERNAL_RAW_FACTOR*raw_floor
        support=_via_site_support_v24(gray,xy,D,nearest_pitch,stage,pattern_key)
        signal_ok=bool(sc>=sc_min and rw>=rw_min)
        if stage=='after' and pattern_key=='via40':
            relaxed_ok=bool(sc>=V24_VIA40_AFTER_SITE_SCORE_FACTOR*sc_min and rw>=V24_VIA40_AFTER_SITE_RAW_FACTOR*rw_min)
            support_ok=bool(support.get('valid',0.0)>0.5 and sc>=V24_VIA40_AFTER_SUPPORT_SCORE_FACTOR*sc_min and rw>=V24_VIA40_AFTER_SUPPORT_RAW_FACTOR*rw_min)
            signal_ok=signal_ok or relaxed_ok or support_ok
        if not signal_ok: continue
        d=_make_circle_desc_v21(xy,D,'v24_lattice_site')
        if d is None: continue
        setattr(d,'_v21_score',float(sc)); setattr(d,'_v21_raw',float(rw)); setattr(d,'_v23_lattice_i',int(site['i'])); setattr(d,'_v23_lattice_j',int(site['j']))
        setattr(d,'_v23_lattice_internal',internal); setattr(d,'_v23_lattice_pred_xy',xy0.copy()); setattr(d,'_v24_site_support',support)
        pre.append(d)
        if not any(float(np.linalg.norm(xy-s.center))<=0.30*D for s in model_seeds): recovered+=1
    pre=_merge_via_center_candidates_v21(pre,D)

    # Expanded rows/columns must have repeated support; isolated off-array dark objects stay rejected.
    if pre:
        ci=Counter(int(getattr(d,'_v23_lattice_i',0)) for d in pre); cj=Counter(int(getattr(d,'_v23_lattice_j',0)) for d in pre)
        imin,imax=int(model['i_min']),int(model['i_max']); jmin,jmax=int(model['j_min']),int(model['j_max']); p2=[]
        for d in pre:
            i=int(getattr(d,'_v23_lattice_i',0)); j=int(getattr(d,'_v23_lattice_j',0))
            if ((imin<=i<=imax) or ci[i]>=3) and ((jmin<=j<=jmax) or cj[j]>=3): p2.append(d)
            else: rejects.append(_make_reject_from_contour(d.contour,'via_v24_grid_qc','isolated_external_lattice_site','expanded lattice site lacks repeated row/column support',desc=d,lattice_i=i,lattice_j=j))
        pre=p2

    survivors=[]; edge_meta_by_id={}; shape_meta_by_id={}
    for d in pre:
        pred=np.asarray(getattr(d,'_v23_lattice_pred_xy',d.center),float); grid_dist=float(np.linalg.norm(d.center-pred))
        grid_tol=max(2.5,V23_VIA_GRID_SITE_TOL_DIAM_FRAC*D)
        if np.isfinite(nearest_pitch): grid_tol=min(grid_tol,0.22*nearest_pitch)
        if grid_dist>grid_tol:
            rejects.append(_make_reject_from_contour(d.contour,'via_v24_grid_qc','via_center_off_lattice',f'center is {grid_dist:.2f}px from predicted grid site (> {grid_tol:.2f}px)',desc=d,lattice_distance_px=grid_dist)); continue

        # Try both the locally refined centre and the theoretical grid centre.  A very
        # dark off-centre core can pull the refined centre inward; the candidate with the
        # most complete plausible outer edge wins.
        centres=[np.asarray(d.center,float)]
        if float(np.linalg.norm(pred-d.center))>0.06*D: centres.append(pred)
        centres.append(0.5*(pred+np.asarray(d.center,float)))
        unique=[]
        for cc in centres:
            if not any(np.linalg.norm(cc-u)<0.5 for u in unique): unique.append(cc)
        options=[]
        for cc in unique:
            cont,em=_via_outer_edge_contour_v24(gray,cc,D,nearest_pitch,stage,pattern_key)
            if cont is None: continue
            ok,smet,reasons=_via_shape_qc_v24(cont,D,stage,pattern_key)
            if not ok: continue
            pts=cont.reshape(-1,2).astype(float); spanx=float(pts[:,0].max()-pts[:,0].min()); spany=float(pts[:,1].max()-pts[:,1].min())
            if min(spanx,spany)<V23_VIA_MIN_SPAN_FACTOR*D or max(spanx,spany)>V23_VIA_MAX_SPAN_FACTOR*D: continue
            leak=float(em.get('outer_dark_leak_fraction',1.0)); eq=float(smet.get('eq_diam_px',0.0))/max(D,1e-9)
            # Larger complete contours are preferred only while shape and leak remain good.
            q=2.0*float(em.get('good_ray_fraction',0.0))+0.65*min(eq,1.35)-1.25*max(0.0,leak-V24_VIA_MAX_OUTER_DARK_LEAK)+0.25*float(smet.get('solidity',0.0))+0.20*float(smet.get('ellipse_iou',0.0))
            options.append((q,cont,em,smet,cc))
        if not options:
            rejects.append(_make_reject_from_contour(d.contour,'via_v24_edge_shape_qc','via_outer_edge_not_reliable','complete outer edge failed sustained-transition, physical-size, or shape QC',desc=d)); continue
        options.sort(key=lambda z:z[0],reverse=True); _,cont,em,smet,used_center=options[0]
        d.contour=cont.astype(np.int32)
        ad=describe_contour(d.contour)
        if ad is not None:
            d.area=ad.area; d.perimeter=ad.perimeter; d.circularity=ad.circularity; d.solidity=ad.solidity; d.u=ad.u; d.v=ad.v
            d.major_min=ad.major_min; d.major_max=ad.major_max; d.minor_min=ad.minor_min; d.minor_max=ad.minor_max; d.length_px=ad.length_px; d.width_px=ad.width_px; d.angle_from_x_deg=ad.angle_from_x_deg
        survivors.append(d); edge_meta_by_id[id(d)]=dict(em,used_edge_center_x_px=float(used_center[0]),used_edge_center_y_px=float(used_center[1])); shape_meta_by_id[id(d)]=smet

    axis=_fit_via_bottom_row_axes_v19(survivors,gray.shape); xa=np.asarray(axis['x_axis'],float); ya=np.asarray(axis['y_axis_up'],float); origin=np.asarray(axis['origin'],float)
    final=[]; rows=[]; bottom_set=set(axis.get('bottom_indices',[]))
    for idx,d in enumerate(survivors):
        wx=float(_projection_span_v19(d.contour,xa)); hy=float(_projection_span_v19(d.contour,ya)); w_nm=wx*px_nm; h_nm=hy*px_nm
        if not (V23_VIA_MIN_SPAN_FACTOR*D*px_nm<=w_nm<=V23_VIA_MAX_SPAN_FACTOR*D*px_nm and V23_VIA_MIN_SPAN_FACTOR*D*px_nm<=h_nm<=V23_VIA_MAX_SPAN_FACTOR*D*px_nm):
            rejects.append(_make_reject_from_contour(d.contour,'via_v24_physical_qc','via_final_span_out_of_range',f'actual X/Y spans {w_nm:.2f}/{h_nm:.2f}nm inconsistent with image search diameter {D*px_nm:.2f}nm',desc=d,x_width_nm=w_nm,y_height_nm=h_nm)); continue
        em=edge_meta_by_id.get(id(d),{}); smet=shape_meta_by_id.get(id(d),{}); support=getattr(d,'_v24_site_support',{})
        final.append(d); rows.append({
            'stage':stage,'stage_zh':STAGE_LABELS_ZH_V19[stage],'condition':condition,'condition_label':CONDITION_LABELS_V19.get(str(condition),f'Condition {condition}'),
            'image':image_name,'pattern':pattern_key,'pattern_zh':PATTERN_LABELS_ZH_V19[pattern_key],'object_id':len(final),'used_for_statistics':True,'is_axis_bottom_row':idx in bottom_set,
            'x_width_nm':w_nm,'y_height_nm':h_nm,'tip_gap_y_nm':np.nan,'x_width_px':wx,'y_height_px':hy,
            'center_image_x_px':float(d.center[0]),'center_image_y_px':float(d.center[1]),'center_axis_x_px':float((d.center-origin)@xa),'center_axis_y_up_px':float((d.center-origin)@ya),
            'circularity':float(smet.get('circularity',d.circularity)),'solidity':float(smet.get('solidity',d.solidity)),'axis_ratio':float(smet.get('aspect',d.aspect)),
            'via_ellipse_iou':float(smet.get('ellipse_iou',np.nan)),'via_radial_cv':float(smet.get('radial_cv',np.nan)),'via_equivalent_diameter_nm':float(smet.get('eq_diam_px',np.nan))*px_nm,
            'selection_source':_source_of(d),'pixel_size_nm':px_nm,'via_nominal_diameter_nm':nominal_nm,'via_nominal_diameter_px':nominal_D,'via_search_diameter_px':D,'via_search_diameter_nm':D*px_nm,
            'via_scale_match_score':float(getattr(d,'_v21_score',np.nan)),'via_scale_raw_contrast_gray':float(getattr(d,'_v21_raw',np.nan)),
            'via_lattice_i':int(getattr(d,'_v23_lattice_i',0)),'via_lattice_j':int(getattr(d,'_v23_lattice_j',0)),
            'via_lattice_site_distance_px':float(np.linalg.norm(d.center-np.asarray(getattr(d,'_v23_lattice_pred_xy',d.center)))),
            'via_axis_angle_from_image_x_deg':axis.get('row_angle_deg',np.nan),'via_axis_fit_rms_px':axis.get('row_fit_rms_px',np.nan),'estimated_center_pitch_px':nearest_pitch,
            'x_edge_method':'v24_outermost_sustained_transition+connected_dark_component','y_edge_method':'v24_outermost_sustained_transition+connected_dark_component',
            'via_good_ray_fraction':em.get('good_ray_fraction',np.nan),'via_edge_smooth_sigma_px':em.get('smooth_sigma_px',np.nan),'via_edge_smooth_passes':em.get('smooth_passes',np.nan),
            'via_outer_dark_leak_fraction':em.get('outer_dark_leak_fraction',np.nan),'via_component_used':em.get('component_used',np.nan),
            'via_site_disk_contrast_gray':support.get('disk_contrast_gray',np.nan),'via_site_shoulder_contrast_gray':support.get('shoulder_contrast_gray',np.nan),
        })
    if final:
        axis=_fit_via_bottom_row_axes_v19(final,gray.shape); xa=np.asarray(axis['x_axis'],float); ya=np.asarray(axis['y_axis_up'],float); origin=np.asarray(axis['origin'],float); bset=set(axis.get('bottom_indices',[]))
        for i,(d,r) in enumerate(zip(final,rows)):
            r['object_id']=i+1; r['is_axis_bottom_row']=i in bset; r['center_axis_x_px']=float((d.center-origin)@xa); r['center_axis_y_up_px']=float((d.center-origin)@ya)
            r['x_width_px']=float(_projection_span_v19(d.contour,xa)); r['y_height_px']=float(_projection_span_v19(d.contour,ya)); r['x_width_nm']=r['x_width_px']*px_nm; r['y_height_nm']=r['y_height_px']*px_nm
            r['via_axis_angle_from_image_x_deg']=axis.get('row_angle_deg',np.nan); r['via_axis_fit_rms_px']=axis.get('row_fit_rms_px',np.nan)
    meta={
        'array_assist_used':True,'array_lattice_type':model.get('lattice_type','rectangular'),'array_expected_positions':len(sites),'array_recovered_count':recovered,
        'array_occupancy_seed':float(model.get('occupancy',np.nan)),'array_basis_a_length_px':float(model.get('basis_a_len',np.nan)),'array_basis_b_length_px':float(model.get('basis_b_len',np.nan)),
        'array_basis_angle_deg':float(model.get('basis_angle_deg',np.nan)),'via_bottom_row_count':len(axis.get('bottom_indices',[])),'via_axis_angle_deg':float(axis.get('row_angle_deg',np.nan)),
        'via_axis_fit_rms_px':float(axis.get('row_fit_rms_px',np.nan)),'via_axis_fallback':bool(axis.get('fallback',False)),'estimated_center_pitch_px':nearest_pitch,
        'via_v24_nominal_diameter_nm':nominal_nm,'via_v24_nominal_diameter_px':nominal_D,'via_search_diameter_px':D,'via_search_diameter_nm':D*px_nm,'via_v24_grid_required':True,'via_v24_lattice_used':True,
        'via_v24_pre_qc_count':len(pre),'via_v24_final_count':len(final),'via_v24_edge_definition':'outermost_sustained_dark_to_bright_plus_connected_component',**resp_meta,
    }
    return rows,final,rejects,{'axis':axis,'meta':meta,'lattice_model':model}



# ============================================================
# V25 ZERO-RESULT AUTOMATIC FALLBACK
# ============================================================
# Philosophy:
#   1) Run V24 exactly as before.
#   2) ONLY if an image returns zero measurements, run relaxed level 1.
#   3) ONLY if level 1 is also zero, run relaxed level 2.
#   4) Structural priors remain mandatory even in fallback:
#        - via: rectangular lattice + physical size + non-pathological shape
#        - slot: complete 4-object columns only
#        - trench: repeated upper/lower paired columns
# This prevents a bad image from turning random dark texture into measurements.
V25_ENABLE_ZERO_RESULT_FALLBACK = True
V25_ZERO_FALLBACK_MAX_LEVEL = 2
V25_ZERO_FALLBACK_BANNER = True
_V25_ZERO_FALLBACK_LEVEL = 0


def _temporary_global_overrides_v25(overrides: Dict[str, object], func, *args, **kwargs):
    """Temporarily change numeric tuning knobs for one emergency fallback call."""
    old = {}
    g = globals()
    for k, v in overrides.items():
        if k in g:
            old[k] = g[k]
            g[k] = v
    try:
        return func(*args, **kwargs)
    finally:
        for k, v in old.items():
            g[k] = v


def _enhance_zero_fallback_gray_v25(gray: np.ndarray, level: int, pattern_key: str) -> np.ndarray:
    """Mild contrast rescue for poorly exposed/blurred SEM images; geometry is unchanged."""
    g = np.asarray(gray, np.uint8)
    if level <= 0:
        return g.copy()
    # Keep enhancement deliberately mild at level 1. Level 2 is stronger but still local.
    clip = 1.55 if level == 1 else 2.20
    tile = 10 if min(g.shape[:2]) >= 700 else 8
    clahe = cv2.createCLAHE(clipLimit=clip, tileGridSize=(tile, tile))
    if pattern_key in ('via40', 'via60'):
        # A tiny median step helps low-SNR via40-after without erasing the actual via outline.
        base = cv2.medianBlur(g, 3) if level >= 2 else g
        eq = clahe.apply(base)
        a = 0.42 if level == 1 else 0.28
        out = cv2.addWeighted(g, a, eq, 1.0-a, 0)
        if level >= 2:
            out = cv2.GaussianBlur(out, (0, 0), 0.55)
        return out
    eq = clahe.apply(g)
    a = 0.50 if level == 1 else 0.32
    return cv2.addWeighted(g, a, eq, 1.0-a, 0)


def _validate_via_model_v23(model: Dict[str, object], nominal_diam_px: float) -> bool:
    """V25 override: strict V24 validation normally; modest relaxation only in zero-result fallback."""
    if not model.get('used', False):
        return False
    la = float(model.get('basis_a_len', np.nan)); lb = float(model.get('basis_b_len', np.nan))
    ang = float(model.get('basis_angle_deg', np.nan)); occ = float(model.get('occupancy', 0.0))
    if not (np.isfinite(la) and np.isfinite(lb) and np.isfinite(ang)):
        return False
    level = int(globals().get('_V25_ZERO_FALLBACK_LEVEL', 0))
    if level <= 0:
        pmin = V21_VIA_LATTICE_MIN_PITCH_DIAM_FACTOR
        pmax = V21_VIA_LATTICE_MAX_PITCH_DIAM_FACTOR
        amin = V21_VIA_LATTICE_MIN_BASIS_ANGLE
        amax = V21_VIA_LATTICE_MAX_BASIS_ANGLE
        min_sites, min_occ = 4, 0.20
    elif level == 1:
        pmin, pmax, amin, amax = 1.02, 8.5, 48.0, 132.0
        min_sites, min_occ = 3, 0.14
    else:
        pmin, pmax, amin, amax = 0.96, 9.0, 45.0, 135.0
        min_sites, min_occ = 3, 0.10
    if min(la, lb) < pmin * nominal_diam_px:
        return False
    if max(la, lb) > pmax * nominal_diam_px:
        return False
    if not (amin <= ang <= amax):
        return False
    return int(model.get('supported_sites', 0)) >= min_sites and occ >= min_occ


def _via_shape_qc_v24(
    contour: np.ndarray,
    nominal_diam_px: float,
    stage: str,
    pattern_key: str,
) -> Tuple[bool, Dict[str, float], List[str]]:
    """V25 override: fallback is a little tolerant to blur, never tolerant to pathological shapes."""
    m = _via_shape_metrics_v23(contour, nominal_diam_px)
    reasons: List[str] = []
    D = max(1e-9, float(nominal_diam_px))
    level = int(globals().get('_V25_ZERO_FALLBACK_LEVEL', 0))

    # Exact V24 thresholds in normal mode.
    if level <= 0:
        if stage == 'after' and pattern_key == 'via40':
            circ_min, sol_min, iou_min, rcv_max, eq_min = 0.43, 0.76, 0.48, 0.31, 0.70
        elif pattern_key == 'via60':
            circ_min, sol_min, iou_min, rcv_max, eq_min = 0.50, 0.80, 0.55, 0.28, 0.72
        else:
            circ_min, sol_min, iou_min, rcv_max, eq_min = (
                V23_VIA_MIN_CIRCULARITY, V23_VIA_MIN_SOLIDITY,
                V23_VIA_MIN_ELLIPSE_IOU, V23_VIA_MAX_RADIAL_CV,
                V23_VIA_MIN_EQ_DIAM_FACTOR,
            )
        aspect_max = V23_VIA_MAX_ASPECT
    elif level == 1:
        # Soft rescue: enough for blur/low contrast, still rejects stars, fragments and tiny blobs.
        if stage == 'after' and pattern_key == 'via40':
            circ_min, sol_min, iou_min, rcv_max, eq_min = 0.39, 0.72, 0.43, 0.35, 0.66
        elif pattern_key == 'via60':
            circ_min, sol_min, iou_min, rcv_max, eq_min = 0.44, 0.75, 0.48, 0.33, 0.68
        else:
            circ_min, sol_min, iou_min, rcv_max, eq_min = 0.43, 0.74, 0.47, 0.33, 0.68
        aspect_max = 1.72
    else:
        # Strong rescue is still bounded: grid membership + size + solidity remain mandatory.
        if stage == 'after' and pattern_key == 'via40':
            circ_min, sol_min, iou_min, rcv_max, eq_min = 0.36, 0.69, 0.39, 0.38, 0.63
        elif pattern_key == 'via60':
            circ_min, sol_min, iou_min, rcv_max, eq_min = 0.40, 0.72, 0.43, 0.36, 0.65
        else:
            circ_min, sol_min, iou_min, rcv_max, eq_min = 0.40, 0.71, 0.43, 0.36, 0.65
        aspect_max = 1.80

    eq_min = min(eq_min, V23_VIA_MIN_SPAN_FACTOR)
    if m['eq_diam_px'] < eq_min * D:
        reasons.append(f"equivalent diameter {m['eq_diam_px']:.2f}px too small")
    if m['eq_diam_px'] > V23_VIA_MAX_EQ_DIAM_FACTOR * D:
        reasons.append(f"equivalent diameter {m['eq_diam_px']:.2f}px too large")
    if m['circularity'] < circ_min:
        reasons.append(f"circularity {m['circularity']:.3f} < {circ_min:.3f}")
    if m['solidity'] < sol_min:
        reasons.append(f"solidity {m['solidity']:.3f} < {sol_min:.3f}")
    if m['aspect'] > aspect_max:
        reasons.append(f"axis ratio {m['aspect']:.3f} > {aspect_max:.3f}")
    if m['ellipse_iou'] < iou_min:
        reasons.append(f"ellipse IoU {m['ellipse_iou']:.3f} < {iou_min:.3f}")
    if m['radial_cv'] > rcv_max:
        reasons.append(f"radial CV {m['radial_cv']:.3f} > {rcv_max:.3f}")
    return len(reasons) == 0, m, reasons


def _fallback_knobs_v25(pattern_key: str, level: int) -> Dict[str, object]:
    """Pattern-specific thresholds. Structural requirements are intentionally not removed."""
    if pattern_key in ('via40', 'via60'):
        if level == 1:
            return {
                'V21_VIA_PEAK_SCORE_MIN': 0.21,
                'V21_VIA_INTERNAL_SITE_SCORE_MIN': 0.22,
                'V21_VIA_EXTERNAL_SITE_SCORE_MIN': 0.29,
                'V21_VIA_EXTERNAL_RAW_FACTOR': 0.92,
                'V23_VIA_MIN_SPAN_FACTOR': 0.45,
                'V23_VIA_MAX_SPAN_FACTOR': 2.10,
                'V24_VIA40_AFTER_SITE_SCORE_FACTOR': 0.60,
                'V24_VIA40_AFTER_SITE_RAW_FACTOR': 0.58,
                'V24_VIA40_AFTER_SUPPORT_SCORE_FACTOR': 0.34,
                'V24_VIA40_AFTER_SUPPORT_RAW_FACTOR': 0.22,
                'V24_VIA40_AFTER_MIN_GOOD_RAYS': 0.42,
                'V24_VIA60_MIN_GOOD_RAYS': 0.48,
                'V24_VIA_OTHER_MIN_GOOD_RAYS': 0.52,
                'V24_VIA_MAX_OUTER_DARK_LEAK': 0.40,
                'V24_VIA40_AFTER_DARKNESS_THRESHOLD': 0.090,
                'V24_VIA60_DARKNESS_THRESHOLD': 0.095,
                'V24_VIA_OTHER_DARKNESS_THRESHOLD': 0.120,
            }
        return {
            'V21_VIA_PEAK_SCORE_MIN': 0.17,
            'V21_VIA_INTERNAL_SITE_SCORE_MIN': 0.18,
            'V21_VIA_EXTERNAL_SITE_SCORE_MIN': 0.24,
            'V21_VIA_EXTERNAL_RAW_FACTOR': 0.84,
            'V23_VIA_MIN_SPAN_FACTOR': 0.40,
            'V23_VIA_MAX_SPAN_FACTOR': 2.20,
            'V24_VIA40_AFTER_SITE_SCORE_FACTOR': 0.52,
            'V24_VIA40_AFTER_SITE_RAW_FACTOR': 0.48,
            'V24_VIA40_AFTER_SUPPORT_SCORE_FACTOR': 0.27,
            'V24_VIA40_AFTER_SUPPORT_RAW_FACTOR': 0.17,
            'V24_VIA40_AFTER_MIN_GOOD_RAYS': 0.36,
            'V24_VIA60_MIN_GOOD_RAYS': 0.42,
            'V24_VIA_OTHER_MIN_GOOD_RAYS': 0.46,
            'V24_VIA_MAX_OUTER_DARK_LEAK': 0.46,
            'V24_VIA40_AFTER_DARKNESS_THRESHOLD': 0.070,
            'V24_VIA60_DARKNESS_THRESHOLD': 0.078,
            'V24_VIA_OTHER_DARKNESS_THRESHOLD': 0.100,
        }
    if pattern_key == 'slot210':
        if level == 1:
            return {
                'SLOT_MIN_AREA_PX': 22,
                'V17_SLOT_SEED_MIN_SOLIDITY': 0.64,
                'V17_SLOT_SEED_MIN_CONVEXITY': 0.48,
                'V17_SLOT_SEED_MIN_TEMPLATE_IOU': 0.32,
                'V17_SLOT_RECOVER_MIN_SOLIDITY': 0.62,
                'V17_SLOT_RECOVER_MIN_CONVEXITY': 0.46,
                'V17_SLOT_RECOVER_MIN_TEMPLATE_IOU': 0.30,
                'V17_SLOT_RECOVER_LOCAL_CONTRAST_GRAY': 0.065,
                'V18_SLOT_LOWGRAY_MODEL_AREA_MIN': 14,
            }
        return {
            'SLOT_MIN_AREA_PX': 16,
            'V17_SLOT_SEED_MIN_SOLIDITY': 0.58,
            'V17_SLOT_SEED_MIN_CONVEXITY': 0.42,
            'V17_SLOT_SEED_MIN_TEMPLATE_IOU': 0.27,
            'V17_SLOT_RECOVER_MIN_SOLIDITY': 0.56,
            'V17_SLOT_RECOVER_MIN_CONVEXITY': 0.40,
            'V17_SLOT_RECOVER_MIN_TEMPLATE_IOU': 0.25,
            'V17_SLOT_RECOVER_LOCAL_CONTRAST_GRAY': 0.040,
            'V18_SLOT_LOWGRAY_MODEL_AREA_MIN': 10,
        }
    if pattern_key == 'trench160':
        if level == 1:
            return {
                'V20_TRENCH_MIN_GROUPS': 3,
                'V20_TRENCH_MIN_PAIR_DARK_CONTRAST': 0.28,
                'V20_TRENCH_GAP_FRAC_TOL': 0.90,
                'V20_TRENCH_WIDTH_FRAC_TOL': 0.95,
            }
        return {
            'V20_TRENCH_MIN_GROUPS': 2,
            'V20_TRENCH_MIN_PAIR_DARK_CONTRAST': 0.10,
            'V20_TRENCH_GAP_FRAC_TOL': 1.00,
            'V20_TRENCH_WIDTH_FRAC_TOL': 1.10,
        }
    return {}


def _fallback_candidate_pool_v25(
    gray: np.ndarray, base_pattern: str, image_name: str
) -> Tuple[List[ShapeDesc], Dict[str, object]]:
    """Rebuild via/slot candidates from the enhanced image during fallback only."""
    if base_pattern not in ('via', 'slot'):
        return [], {}
    try:
        mask = segment_features(gray)
        descs, _ = extract_descriptors_with_rejections(mask)
        pool, meta = build_candidate_pool_v12(gray, descs, base_pattern, image_name)
        return pool, meta
    except Exception as e:
        return [], {'fallback_candidate_pool_error': f'{type(e).__name__}: {e}'}


def _run_zero_fallback_v25(
    gray: np.ndarray,
    original_pool: List[ShapeDesc],
    px_nm: float,
    image_name: str,
    stage: str,
    condition: str,
    pattern_key: str,
    level: int,
):
    """Run one fallback level and return a pattern-specific result bundle."""
    global _V25_ZERO_FALLBACK_LEVEL
    base_pattern = _base_pattern_v19(pattern_key)
    fgray = _enhance_zero_fallback_gray_v25(gray, level, pattern_key)
    knobs = _fallback_knobs_v25(pattern_key, level)
    old_level = _V25_ZERO_FALLBACK_LEVEL
    _V25_ZERO_FALLBACK_LEVEL = int(level)
    try:
        if base_pattern in ('via', 'slot'):
            # Level 1 keeps the original pool; level 2 also rebuilds candidates from enhanced contrast.
            pool = original_pool
            pool_meta = {}
            if level >= 2:
                rebuilt, pool_meta = _fallback_candidate_pool_v25(fgray, base_pattern, image_name)
                if rebuilt:
                    pool = rebuilt
        else:
            pool = original_pool
            pool_meta = {}

        if base_pattern == 'via':
            result = _temporary_global_overrides_v25(
                knobs, _select_and_measure_vias_v19,
                fgray, pool, px_nm, image_name, stage, condition, pattern_key,
            )
            rows, selected, rejects, bundle = result
            bundle = dict(bundle)
            bundle['meta'] = dict(bundle.get('meta', {}), **pool_meta)
            return rows, (selected, bundle), rejects, bundle['meta'], fgray
        if base_pattern == 'slot':
            result = _temporary_global_overrides_v25(
                knobs, _select_and_measure_bottom_slots_v19,
                fgray, pool, px_nm, image_name, stage, condition, pattern_key,
            )
            rows, selected_all, bottom, rejects, meta = result
            meta = dict(meta, **pool_meta)
            return rows, (selected_all, bottom), rejects, meta, fgray
        if base_pattern == 'trench':
            result = _temporary_global_overrides_v25(
                knobs, _select_and_measure_trenches_v19,
                fgray, [], px_nm, image_name, stage, condition, pattern_key,
            )
            rows, pairs, rejects, meta = result
            return rows, pairs, rejects, meta, fgray
        return [], None, [], {}, fgray
    finally:
        _V25_ZERO_FALLBACK_LEVEL = old_level


def _tag_rows_fallback_v25(rows: List[dict], level: int) -> None:
    for r in rows:
        r['auto_fallback_used'] = bool(level > 0)
        r['auto_fallback_level'] = int(level)
        r['selection_mode'] = 'strict_v24' if level <= 0 else f'zero_result_fallback_L{level}'
        r['measurement_confidence_flag'] = 'NORMAL' if level <= 0 else 'RELAXED_FALLBACK_REVIEW'


def _add_fallback_banner_v25(out: np.ndarray, level: int) -> np.ndarray:
    if not V25_ZERO_FALLBACK_BANNER or level <= 0:
        return out
    H, W = out.shape[:2]
    msg = f'AUTO FALLBACK L{level}: strict pass found 0; relaxed thresholds used'
    y0 = max(0, H - 25)
    cv2.rectangle(out, (3, y0), (min(W-4, 610), H-3), (0, 0, 0), -1)
    cv2.putText(out, msg, (7, H-8), cv2.FONT_HERSHEY_SIMPLEX, 0.38, (0, 220, 255), 1, cv2.LINE_AA)
    return out


def _process_one_image_v19(
    image_path: Path,
    stage: str,
    condition: str,
    ann_dir: Path,
) -> Tuple[List[dict], dict, List[dict]]:
    """V25 override: V24 strict first; automatic fallback only when strict result is empty."""
    pattern_key = _classify_pattern_key_v19(image_path)
    if pattern_key is None:
        return [], {
            'stage': stage, 'stage_zh': STAGE_LABELS_ZH_V19[stage], 'condition': condition,
            'condition_label': CONDITION_LABELS_V19.get(str(condition), f'Condition {condition}'),
            'image': image_path.name, 'path': str(image_path), 'pattern': 'unknown',
            'status': 'SKIPPED_UNKNOWN_PREFIX', 'note': 'Filename must begin with 40, 60, 160, or 210',
        }, []

    base_pattern = _base_pattern_v19(pattern_key)
    print(f'    {image_path.name} -> {pattern_key}', flush=True)
    px_info = read_pixel_size_nm(image_path); px_nm = px_info.value_nm
    gray = _read_gray_image_unicode_v19(image_path)
    mask = segment_features(gray)
    primary_descs, segmentation_rejects = extract_descriptors_with_rejections(mask)
    search_meta: Dict[str, object] = {}
    measure_meta: Dict[str, object] = _empty_array_meta_v19()

    if base_pattern in {'via', 'slot'}:
        candidate_pool, search_meta = build_candidate_pool_v12(gray, primary_descs, base_pattern, image_path.name)
    elif base_pattern == 'trench':
        candidate_pool = []
        search_meta = {'trench_search_mode': 'V20_gray_profile_direct'}
        segmentation_rejects = []
    else:
        candidate_pool = list(primary_descs)

    rows: List[dict] = []
    strict_rejects: List[dict] = []
    selected_payload = None

    # ---------------- strict V24 pass ----------------
    if base_pattern == 'via':
        rows, selected, strict_rejects, bundle = _select_and_measure_vias_v19(
            gray, candidate_pool, px_nm, image_path.name, stage, condition, pattern_key)
        measure_meta = bundle['meta']; selected_payload = (selected, bundle)
    elif base_pattern == 'trench':
        rows, pairs, strict_rejects, measure_meta = _select_and_measure_trenches_v19(
            gray, candidate_pool, px_nm, image_path.name, stage, condition, pattern_key)
        selected_payload = pairs
    elif base_pattern == 'slot':
        rows, selected_all, bottom, strict_rejects, measure_meta = _select_and_measure_bottom_slots_v19(
            gray, candidate_pool, px_nm, image_path.name, stage, condition, pattern_key)
        selected_payload = (selected_all, bottom)

    fallback_level = 0
    fallback_rejects: List[dict] = []
    fallback_attempts: List[str] = []

    # ---------------- zero-result fallback ----------------
    if V25_ENABLE_ZERO_RESULT_FALLBACK and not rows and base_pattern in ('via', 'slot', 'trench'):
        print('      strict pass selected 0 -> starting automatic fallback', flush=True)
        for level in range(1, int(V25_ZERO_FALLBACK_MAX_LEVEL) + 1):
            fallback_attempts.append(f'L{level}')
            fb = _run_zero_fallback_v25(
                gray, candidate_pool, px_nm, image_path.name, stage, condition, pattern_key, level)
            fb_rows, payload, fb_rej, fb_meta, _ = fb
            fallback_rejects.extend(fb_rej)
            print(f'      fallback L{level}: measurements={len(fb_rows)}', flush=True)
            if fb_rows:
                rows = fb_rows; selected_payload = payload; measure_meta = fb_meta
                fallback_level = level
                break

    _tag_rows_fallback_v25(rows, fallback_level)

    # Annotation always uses the original SEM image, never the enhanced fallback image.
    if base_pattern == 'via':
        if selected_payload is None:
            selected, bundle = [], {'axis': _fit_via_bottom_row_axes_v19([], gray.shape), 'meta': {}}
        else:
            selected, bundle = selected_payload
        annotated = _annotate_vias_v19(gray, selected, rows, bundle)
        n_selected = len(selected)
    elif base_pattern == 'trench':
        pairs = selected_payload if selected_payload is not None else []
        annotated = _annotate_trenches_v19(gray, pairs, rows)
        n_selected = 2 * len(pairs)
    elif base_pattern == 'slot':
        if selected_payload is None:
            selected_all, bottom = [], []
        else:
            selected_all, bottom = selected_payload
        annotated = _annotate_slots_v19(gray, selected_all, bottom, rows)
        n_selected = len(bottom)
    else:
        annotated = to_color(gray); n_selected = 0

    if fallback_level > 0:
        annotated = _add_fallback_banner_v25(annotated, fallback_level)

    all_rejects = segmentation_rejects + strict_rejects + fallback_rejects
    exported_rejects: List[dict] = []
    for rid, r in enumerate(all_rejects, 1):
        r['reject_id'] = rid; r['stage'] = stage; r['condition'] = condition
        r['image'] = image_path.name; r['path'] = str(image_path); r['pattern'] = pattern_key
        r['pixel_size_nm'] = px_nm
        exported_rejects.append(exportable_rejection_record(r))
    if DRAW_REJECTED_OBJECTS_V19:
        annotated = annotate_rejections(annotated, all_rejects)

    safe = _safe_stem_v19(image_path)
    ann_path = ann_dir / f'{safe}__{pattern_key}__annotated.png'
    _imwrite_unicode_v19(ann_path, annotated)
    if SAVE_BINARY_MASK_V19:
        _imwrite_unicode_v19(ann_dir / f'{safe}__primary_mask.png', mask)

    if rows and fallback_level == 0:
        status_code = 'OK'
    elif rows and fallback_level > 0:
        status_code = f'OK_AUTO_FALLBACK_L{fallback_level}'
    else:
        status_code = 'CHECK_NO_VALID_MEASUREMENT_AFTER_FALLBACK' if fallback_attempts else 'CHECK_NO_VALID_MEASUREMENT'

    note = ''
    if fallback_level > 0:
        note = f'Strict V24 pass returned zero; automatic fallback level {fallback_level} produced measurements. Review annotation.'
    elif fallback_attempts and not rows:
        note = f'Strict pass and fallback levels {",".join(fallback_attempts)} all returned zero; no measurement was forced.'

    status = {
        'stage': stage, 'stage_zh': STAGE_LABELS_ZH_V19[stage], 'condition': condition,
        'condition_label': CONDITION_LABELS_V19.get(str(condition), f'Condition {condition}'),
        'image': image_path.name, 'path': str(image_path), 'pattern': pattern_key,
        'pattern_zh': PATTERN_LABELS_ZH_V19[pattern_key], 'pixel_size_nm': px_nm,
        'pixelsize_txt': str(px_info.txt_path), 'pixelsize_encoding': px_info.encoding,
        'pixelsize_matched_line': px_info.matched_line, 'image_height_px': int(gray.shape[0]),
        'image_width_px': int(gray.shape[1]), 'n_primary_candidates': len(primary_descs),
        'n_candidate_pool': len(candidate_pool), 'n_selected_structures': n_selected,
        'n_measurements': len(rows), 'n_rejected_total': len(all_rejects),
        'rejected_reason_summary': format_rejection_summary(all_rejects), 'annotated_path': str(ann_path),
        'status': status_code, 'note': note,
        'auto_fallback_enabled': bool(V25_ENABLE_ZERO_RESULT_FALLBACK),
        'auto_fallback_triggered': bool(fallback_attempts),
        'auto_fallback_level_used': int(fallback_level),
        'auto_fallback_attempts': ','.join(fallback_attempts),
    }
    status.update({f'search_{k}': v for k, v in search_meta.items()})
    status.update({f'measure_{k}': v for k, v in measure_meta.items() if np.isscalar(v) or isinstance(v, (str, bool))})
    return rows, status, exported_rejects





# Yakun V7 primitives from the user-provided archive; see ALGORITHM_V32.md.

def _yakun_v111_moving_average(profile: np.ndarray, window: int = 3) -> np.ndarray:
    """Exact V1.11/V13 reflected moving-average convention."""
    a = np.asarray(profile, dtype=np.float64)
    w = max(1, int(window))
    if w <= 1 or a.size < 3:
        return a.copy()
    if w % 2 == 0:
        w += 1
    w = min(w, a.size if a.size % 2 else max(1, a.size - 1))
    pad = w // 2
    return np.convolve(np.pad(a, (pad, pad), mode="reflect"), np.ones(w) / w, mode="valid")


def _yakun_v111_gradient(profile: np.ndarray) -> np.ndarray:
    a = np.asarray(profile, dtype=np.float64)
    g = np.empty_like(a)
    if a.size < 3:
        g[:] = 0.0
        return g
    g[1:-1] = 0.5 * (a[2:] - a[:-2])
    g[0] = a[1] - a[0]
    g[-1] = a[-1] - a[-2]
    return g


def _yakun_v111_mean(profile: np.ndarray, start: float, end: float) -> float:
    n = len(profile)
    s = max(0, min(n, int(math.floor(start))))
    e = max(0, min(n, int(math.ceil(end))))
    if e < s:
        s, e = e, s
    return float(np.mean(profile[s:e])) if e > s else math.nan


def _yakun_v111_median(profile: np.ndarray, start: float, end: float) -> float:
    n = len(profile)
    s = max(0, min(n, int(math.floor(start))))
    e = max(0, min(n, int(math.ceil(end))))
    if e < s:
        s, e = e, s
    return float(np.median(profile[s:e])) if e > s else math.nan


def _yakun_v111_local_level(profile: np.ndarray, index: int) -> float:
    return _yakun_v111_median(profile, index - 1, index + 2)


def _yakun_v111_mad(values: np.ndarray) -> float:
    a = np.asarray(values, dtype=float)
    a = a[np.isfinite(a)]
    if not len(a):
        return 0.0
    med = float(np.median(a))
    return float(1.4826 * np.median(np.abs(a - med)))


def _yakun_v111_topology(
    profile: np.ndarray, edge: float, side: str, band: int, gap: int
) -> Tuple[float, float]:
    x = int(round(edge))
    if side == "left":
        outside = _yakun_v111_mean(profile, x - gap - band, x - gap)
        inside = _yakun_v111_mean(profile, x + gap, x + gap + band)
    else:
        inside = _yakun_v111_mean(profile, x - gap - band, x - gap)
        outside = _yakun_v111_mean(profile, x + gap, x + gap + band)
    return outside, inside


def _yakun_v111_peak(
    profile: np.ndarray,
    gradient: np.ndarray,
    start: int,
    end: int,
    side: str,
    center: float,
    anchor: float,
    max_anchor_distance: float,
    dark_level: float,
    topology_band: int,
) -> Tuple[Optional[int], float, bool]:
    start, end = max(1, int(start)), min(len(profile) - 1, int(end))
    if end - start < 3:
        return None, math.nan, False
    transformed = -gradient if side == "left" else gradient
    local_abs = np.abs(gradient[start:end])
    threshold = max(float(np.median(local_abs)) + _yakun_v111_mad(local_abs), 1e-6)

    def valid(p: int) -> bool:
        if side == "left" and p >= center:
            return False
        if side == "right" and p <= center:
            return False
        if abs(float(p) - anchor) > max_anchor_distance:
            return False
        if _yakun_v111_local_level(profile, p) < dark_level + 2.0:
            return False
        outside, inside = _yakun_v111_topology(profile, p, side, topology_band, 1)
        return bool(np.isfinite(outside) and np.isfinite(inside) and outside - inside >= 2.5)

    peaks = [p for p in range(start + 1, end - 1)
             if transformed[p] > 0 and transformed[p] >= transformed[p - 1]
             and transformed[p] > transformed[p + 1] and transformed[p] >= threshold and valid(p)]
    peaks.sort(key=lambda p: (abs(float(p) - anchor), abs(float(p) - center)))
    if peaks:
        p = int(peaks[0])
        return p, float(abs(gradient[p])), False
    eligible = [p for p in range(start, end) if transformed[p] > 0 and valid(p)]
    if not eligible:
        return None, math.nan, False
    p = max(eligible, key=lambda q: (abs(float(gradient[q])), -abs(float(q) - anchor)))
    return (int(p), float(abs(gradient[p])), True) if abs(float(gradient[p])) >= threshold else (None, math.nan, False)


def _yakun_v111_crossing(
    profile: np.ndarray, threshold: float, peak: int, center: float, side: str
) -> Optional[float]:
    if side == "left":
        stop = min(len(profile) - 1, int(math.ceil(center)))
        # Include the two samples around the discrete gradient maximum.  With a
        # symmetric 3-point smoother the sub-pixel 50% crossing can lie one
        # sample outside the central-difference peak.
        for x in range(max(0, peak - 2), stop):
            p0, p1 = float(profile[x]), float(profile[x + 1])
            if p0 >= threshold and p1 < threshold:
                den = p1 - p0
                return float(x + (threshold - p0) / den) if abs(den) > 1e-12 else None
    else:
        stop = max(0, int(math.floor(center)))
        for x in range(min(len(profile) - 2, peak + 2), stop - 1, -1):
            p0, p1 = float(profile[x]), float(profile[x + 1])
            if p0 < threshold and p1 >= threshold:
                den = p1 - p0
                return float(x + (threshold - p0) / den) if abs(den) > 1e-12 else None
    return None


def _yakun_v111_edge_pair(
    raw_profile: np.ndarray, predicted_left: float, predicted_right: float
) -> Dict[str, object]:
    """V1.11 shared-edge rule: anchored peak + topology + V1.7 side-band 50%."""
    p = _yakun_v111_moving_average(raw_profile, 3)
    width = float(predicted_right - predicted_left)
    center = 0.5 * (predicted_left + predicted_right)
    if not np.isfinite(width) or width < 4.0:
        return {"valid": False, "reason": "invalid_predicted_width"}
    search = max(4.0, 0.52 * width)
    topology_band = max(2, min(7, int(round(0.12 * width))))
    dark = _yakun_v111_median(p, predicted_left + 0.30 * width, predicted_right - 0.30 * width)
    grad = _yakun_v111_gradient(p)
    lp, lstrength, lfb = _yakun_v111_peak(
        p, grad, int(predicted_left - search), int(predicted_left + search), "left",
        center, predicted_left, search, dark, topology_band)
    rp, rstrength, rfb = _yakun_v111_peak(
        p, grad, int(predicted_right - search), int(predicted_right + search), "right",
        center, predicted_right, search, dark, topology_band)
    if lp is None or rp is None:
        return {"valid": False, "reason": "no_valid_gradient_peak_near_anchor",
                "left_peak": lp, "right_peak": rp}
    dark_ref = _yakun_v111_mean(p, predicted_left + 0.28 * width, predicted_right - 0.28 * width)
    band = max(3, int(round(0.18 * max(width, 4.0))))
    lb = _yakun_v111_mean(p, lp - band, lp - 1)
    rb = _yakun_v111_mean(p, rp + 1, rp + band)
    if not all(np.isfinite(x) for x in (dark_ref, lb, rb)) or lb <= dark_ref or rb <= dark_ref:
        return {"valid": False, "reason": "invalid_v111_side_band_reference"}
    lt, rt = dark_ref + 0.5 * (lb - dark_ref), dark_ref + 0.5 * (rb - dark_ref)
    le = _yakun_v111_crossing(p, lt, lp, center, "left")
    re = _yakun_v111_crossing(p, rt, rp, center, "right")
    if le is None or re is None or not (le < center < re):
        return {"valid": False, "reason": "v111_50pct_crossing_not_found"}
    measured = float(re - le)
    if not 0.55 * width <= measured <= 1.55 * width:
        return {"valid": False, "reason": "measured_width_outside_v111_hard_range"}
    return {
        "valid": True, "left": float(le), "right": float(re), "width": measured,
        "left_peak": float(lp), "right_peak": float(rp),
        "left_peak_strength": lstrength, "right_peak_strength": rstrength,
        "left_peak_fallback": lfb, "right_peak_fallback": rfb,
        "dark_reference": dark_ref, "left_bright_reference": lb, "right_bright_reference": rb,
        "left_threshold": lt, "right_threshold": rt,
        "threshold_reference_mode": "v111_v17_side_band_mean50",
    }


def _yakun_v111_one_sided_dark_edge(raw_profile: np.ndarray, anchor: float, side: str, scale: float) -> Dict[str, object]:
    p = _yakun_v111_moving_average(raw_profile, 3)
    g = _yakun_v111_gradient(p)
    sign_g = -g if side == "left" else g
    search = max(5, int(round(0.40 * scale)))
    start, end = max(2, int(anchor - search)), min(len(p) - 2, int(anchor + search))
    band = max(3, min(7, int(round(0.18 * scale))))
    local = np.abs(g[start:end])
    threshold_g = max(float(np.median(local)) + _yakun_v111_mad(local), 1e-6) if len(local) else math.inf

    def refs(j: int) -> Tuple[float, float]:
        if side == "right":
            dark = _yakun_v111_mean(p, j - 2 * band, j - band)
            bright = _yakun_v111_mean(p, j + 1, j + 1 + band)
        else:
            bright = _yakun_v111_mean(p, j - band, j - 1)
            dark = _yakun_v111_mean(p, j + band, j + 2 * band)
        return dark, bright

    peaks = []
    for j in range(start + 1, end - 1):
        dark, bright = refs(j)
        if sign_g[j] > 0 and sign_g[j] >= sign_g[j - 1] and sign_g[j] > sign_g[j + 1] and sign_g[j] >= threshold_g \
                and np.isfinite(dark) and np.isfinite(bright) and bright - dark >= 2.5:
            peaks.append(j)
    peaks.sort(key=lambda j: abs(j - anchor))
    if not peaks:
        eligible = []
        for j in range(start, end):
            dark, bright = refs(j)
            if sign_g[j] > 0 and np.isfinite(dark) and np.isfinite(bright) and bright - dark >= 2.5:
                eligible.append(j)
        if eligible:
            candidate = max(eligible, key=lambda j: (float(sign_g[j]), -abs(j - anchor)))
            if sign_g[candidate] >= threshold_g:
                peaks = [candidate]
        if not peaks:
            return {"valid": False, "reason": "no_one_sided_v111_peak"}
    peak = peaks[0]
    dark, bright = refs(peak)
    level = dark + 0.5 * (bright - dark)
    crossing = None
    if side == "right":
        for j in range(min(len(p) - 2, peak + 2), max(0, peak - search), -1):
            if p[j] < level <= p[j + 1]:
                crossing = float(j + (level - p[j]) / (p[j + 1] - p[j]))
                break
    else:
        for j in range(max(0, peak - 2), min(len(p) - 1, peak + search)):
            if p[j] >= level > p[j + 1]:
                crossing = float(j + (level - p[j]) / (p[j + 1] - p[j]))
                break
    return {"valid": crossing is not None, "edge": crossing, "peak": peak, "dark_reference": dark,
            "bright_reference": bright, "threshold": level,
            "reason": "" if crossing is not None else "one_sided_v111_crossing_not_found"}


def _yakun_v4_dynamic_track(
    candidate_lists: List[List[Tuple[float, float, bool, int]]],
    baseline: np.ndarray, search_half_width: float
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Circular first-order dynamic programming with an explicit jump penalty."""
    candidates = []
    for i, row in enumerate(candidate_lists):
        if row:
            candidates.append(row)
        else:
            candidates.append([(float(baseline[i]), 7.5, False, 0)])
    maximum_soft_jump = max(1.5, 0.45 * search_half_width)
    best_total = math.inf
    best_path: Optional[List[int]] = None
    for start_index in range(len(candidates[0])):
        previous = np.full(len(candidates[0]), np.inf, float)
        previous[start_index] = float(candidates[0][start_index][1])
        back: List[Optional[np.ndarray]] = [None]
        for i in range(1, len(candidates)):
            current = np.full(len(candidates[i]), np.inf, float)
            current_back = np.zeros(len(candidates[i]), np.int32)
            for k, current_item in enumerate(candidates[i]):
                current_radius = float(current_item[0])
                transition = []
                for j, previous_item in enumerate(candidates[i - 1]):
                    delta = abs(current_radius - float(previous_item[0]))
                    penalty = 0.35 * delta + 2.2 * max(0.0, delta - maximum_soft_jump) ** 2
                    transition.append(previous[j] + penalty)
                jbest = int(np.argmin(transition))
                current[k] = float(current_item[1]) + float(transition[jbest])
                current_back[k] = jbest
            previous = current
            back.append(current_back)
        closure = []
        first_radius = float(candidates[0][start_index][0])
        for j, item in enumerate(candidates[-1]):
            delta = abs(float(item[0]) - first_radius)
            closure.append(previous[j] + 0.35 * delta + 2.2 * max(0.0, delta - maximum_soft_jump) ** 2)
        end_index = int(np.argmin(closure))
        total = float(closure[end_index])
        if total < best_total:
            path = [0] * len(candidates)
            path[-1] = end_index
            for i in range(len(candidates) - 1, 0, -1):
                path[i - 1] = int(back[i][path[i]])
            if path[0] != start_index:
                continue
            best_total, best_path = total, path
    if best_path is None:
        raise RuntimeError("Circular threshold-candidate tracking failed")
    radii = np.array([float(candidates[i][best_path[i]][0]) for i in range(len(candidates))], float)
    real = np.array([bool(candidates[i][best_path[i]][2]) for i in range(len(candidates))], bool)
    return radii, real, best_total

# V32: V27 localization + Yakun signal-supported edges + V30 missing-site recovery.
import copy as _copy_v32

_V32_MEASUREMENT_GRAY = None
_V32_PROCESS_BASE = _process_one_image_v19
_V32_SLOT_BASE = _select_and_measure_bottom_slots_v19
_V32_VIA_BASE = _select_and_measure_vias_v19
_V32_TRENCH_BASE = _select_and_measure_trenches_v19


def _v32_raw_gray(gray):
    return _V32_MEASUREMENT_GRAY if _V32_MEASUREMENT_GRAY is not None else gray


def _v32_read_measurement_gray(path):
    """Decode measurement pixels without percentile clipping or quantization.

    Keep 8-bit gray values exact. Higher bit depths use an affine float scale
    only, so local half-gray crossings retain the original intensity geometry.
    """
    data = np.fromfile(str(path), dtype=np.uint8)
    decoded = cv2.imdecode(data, cv2.IMREAD_UNCHANGED)
    if decoded is None:
        raise RuntimeError(f'Cannot read image: {path}')
    if decoded.ndim == 3:
        code = cv2.COLOR_BGRA2GRAY if decoded.shape[2] == 4 else cv2.COLOR_BGR2GRAY
        decoded = cv2.cvtColor(decoded, code)
    source = decoded.astype(np.float32)
    if not np.all(np.isfinite(source)):
        raise ValueError(f'Non-finite measurement pixels: {path}')
    if decoded.dtype != np.uint8:
        low, high = np.percentile(source, [1., 99.])
        if high <= low:
            low, high = float(source.min()), float(source.max())
        if high > low:
            source = (source-low)*np.float32(255./(high-low))
    return source


def _process_one_image_v19(image_path, stage, condition, ann_dir):
    global _V32_MEASUREMENT_GRAY
    previous = _V32_MEASUREMENT_GRAY
    _V32_MEASUREMENT_GRAY = _v32_read_measurement_gray(image_path)
    try:
        return _V32_PROCESS_BASE(image_path, stage, condition, ann_dir)
    finally:
        _V32_MEASUREMENT_GRAY = previous


def _recover_local_slot_v18_deep(gray, expected_xy, x_pitch, row_pitch,
                                expected_w, expected_h, area_limits, aggressive_level):
    # Preserve V27's accepted shapes. V30 only supplies a missing-site proposal;
    # every final contour is independently checked against original image pixels.
    d = _recover_local_slot_v18_legacy(gray, expected_xy, x_pitch, row_pitch,
                                      expected_w, expected_h, area_limits, aggressive_level)
    if d is not None:
        return d
    d = _v32_recover_flatfield_slot(gray, expected_xy, x_pitch, row_pitch,
                                   expected_w, expected_h, area_limits)
    if d is not None:
        _mark_source(d, 'V30_missing_site_proposal')
    return d


def _v32_contour_samples(contour):
    points = np.asarray(contour, float).reshape(-1, 2)
    points = points[np.linalg.norm(points - np.roll(points, 1, axis=0), axis=1) > 1e-6]
    if len(points) < 3:
        return None
    lengths = np.linalg.norm(np.roll(points, -1, axis=0) - points, axis=1)
    perimeter = float(np.sum(lengths))
    if perimeter < 12:
        return None
    count = int(np.clip(math.ceil(perimeter / 1.4), 48, 240))
    cumulative = np.r_[0., np.cumsum(lengths)]
    distance = np.arange(count, dtype=float) * perimeter / count
    indices = np.minimum(np.searchsorted(cumulative, distance, side='right') - 1, len(points) - 1)
    fraction = (distance - cumulative[indices]) / lengths[indices]
    sampled = points[indices] + fraction[:, None] * (points[(indices + 1) % len(points)] - points[indices])
    tangent = np.roll(sampled, -1, axis=0) - np.roll(sampled, 1, axis=0)
    tangent /= np.maximum(np.linalg.norm(tangent, axis=1)[:, None], 1e-12)
    normal = np.column_stack([tangent[:, 1], -tangent[:, 0]])
    if cv2.contourArea(points.astype(np.float32).reshape(-1, 1, 2), oriented=True) < 0:
        normal *= -1
    return sampled, normal, tangent


def _v32_max_cyclic_gap(valid):
    run = longest = 0
    for missing in np.tile(~np.asarray(valid, bool), 2):
        run = run + 1 if missing else 0
        longest = max(longest, run)
    return min(len(valid), longest)


def _v32_edge_source(gray, contour, extent):
    """Subtract a robust exterior background plane, without flattening the feature.

    All fitting samples are original pixels outside the proposed contour and its
    blurred rim. The plane removes illumination tilt, not the dark/bright edge.
    """
    radius = max(8, int(math.ceil(extent + 4)))
    bx, by, bw, bh = cv2.boundingRect(np.asarray(contour, np.float32))
    x0, y0 = max(0, bx-radius), max(0, by-radius)
    x1, y1 = min(gray.shape[1], bx+bw+radius), min(gray.shape[0], by+bh+radius)
    roi = gray[y0:y1, x0:x1].astype(np.float32)
    shifted = np.rint(np.asarray(contour).reshape(-1, 1, 2) - [x0, y0]).astype(np.int32)
    mask = np.zeros(roi.shape, np.uint8)
    cv2.fillPoly(mask, [shifted], 1)
    excluded = cv2.dilate(mask, np.ones((7, 7), np.uint8))
    exterior = excluded == 0
    yy, xx = np.indices(roi.shape, dtype=float)
    design = np.column_stack([np.ones(np.count_nonzero(exterior)), xx[exterior], yy[exterior]])
    values = roi[exterior].astype(float)
    stride = max(1, len(values)//3000)
    design, values = design[::stride], values[::stride]
    keep = np.ones(len(values), bool)
    coefficient = np.zeros(3)
    for _ in range(4):
        if np.count_nonzero(keep) < 20:
            break
        coefficient = np.linalg.lstsq(design[keep], values[keep], rcond=None)[0]
        residual = values - design @ coefficient
        center = np.median(residual)
        scale = max(.15, 1.4826 * np.median(np.abs(residual-center)))
        keep = np.abs(residual-center) <= 2.5*scale
    roi = roi - (coefficient[1]*xx + coefficient[2]*yy).astype(np.float32)
    if FEATURE_POLARITY.lower() == 'bright':
        roi = 255. - roi
    return roi, np.array([x0, y0]), coefficient[1:]


def _v32_supported_contour(gray, contour):
    """Follow local gray half-level crossings without replacing the shape by a template.

    Adapt Yakun's side-band, signed topology and subpixel-crossing ideas to
    contour normals, so rectangular and noncircular boundaries remain intact.
    Yakun's original closed-loop dynamic tracker selects a continuous path.
    """
    samples = _v32_contour_samples(contour)
    if samples is None:
        return None, dict(edge_status='INVALID_COARSE_CONTOUR')
    points, normals, tangents = samples
    spans = np.ptp(np.asarray(contour).reshape(-1, 2), axis=0)
    width = float(min(spans))
    search = float(np.clip(.16 * width, 2., 6.))
    band = float(np.clip(.12 * width, 1.5, 4.))
    extent = max(search + 1., band + 2.5)
    offsets = np.arange(-extent, extent + .25, .5)
    h, w = gray.shape
    coordinates = points[:, None, :] + offsets[None, :, None] * normals[:, None, :]
    margin_valid = np.all((coordinates[:, :, 0] >= 1) & (coordinates[:, :, 0] < w - 2)
                          & (coordinates[:, :, 1] >= 1) & (coordinates[:, :, 1] < h - 2), axis=1)
    profiles = np.zeros(coordinates.shape[:2], dtype=float)
    source, roi_origin, background_slope = _v32_edge_source(gray, contour, extent)
    for lateral in (-.5, 0., .5):
        coords = coordinates + lateral * tangents[:, None, :] - roi_origin
        profiles += cv2.remap(source, coords[:, :, 0].astype(np.float32),
                              coords[:, :, 1].astype(np.float32), cv2.INTER_LINEAR,
                              borderMode=cv2.BORDER_REFLECT101) / 3.
    # Keep references away from the blurred transition. A coarse threshold
    # contour can sit outside the half-gray edge; sampling only one pixel inside
    # it would bias the dark reference bright and enlarge the final shape.
    inner_limit = max(1.5, .5*band)
    # Opposing robust quantiles estimate the two plateaus when a coarse contour
    # leaves part of one side band inside the blur transition. Symmetric 25/75
    # references avoid treating a mixed transition band as the dark plateau.
    inside = np.percentile(profiles[:, (offsets >= -band - 2.) & (offsets <= -inner_limit)], 25, axis=1)
    outside = np.percentile(profiles[:, (offsets >= inner_limit) & (offsets <= band + 2.)], 75, axis=1)
    contrast = outside - inside
    difference = np.diff(profiles, n=2, axis=1)
    noise = 1.4826 * np.median(np.abs(difference - np.median(difference, axis=1)[:, None]), axis=1) / math.sqrt(6.)
    levels = .5 * (inside + outside)
    candidates = []
    anchor_supported = np.zeros(len(points), bool)
    for i, profile in enumerate(profiles):
        options = []
        if margin_valid[i] and contrast[i] >= max(2.5, 4. * noise[i]):
            for j in range(2, len(offsets) - 3):
                if not profile[j] < levels[i] <= profile[j + 1]:
                    continue
                left = float(np.mean(profile[j - 2:j + 1]))
                right = float(np.mean(profile[j + 1:j + 4]))
                if not left <= levels[i] <= right or right - left < max(.15 * contrast[i], 2. * noise[i], .3):
                    continue
                shift = float(offsets[j] + .5 * (levels[i] - profile[j]) / (profile[j + 1] - profile[j]))
                if abs(shift) > search:
                    continue
                cost = .35 * abs(shift) / search + .2 * max(noise[i], .05) / max(right - left, .1)
                options.append((shift, float(cost), True, len(options)))
            anchor_supported[i] = any(abs(option[0]) <= .85 for option in options)
        candidates.append(sorted(options, key=lambda item: item[1])[:4])
    valid = np.array([bool(items) for items in candidates])
    support = float(np.mean(valid))
    gap = _v32_max_cyclic_gap(valid)
    audit = dict(edge_status='INSUFFICIENT_ORIGINAL_IMAGE_SUPPORT', edge_support_fraction=support,
                 edge_max_missing_fraction=gap / len(points), edge_threshold_percent=50.,
                 edge_reference_inside=float(np.median(inside)), edge_reference_outside=float(np.median(outside)),
                 edge_contrast_gray=float(np.median(contrast)), edge_noise_gray=float(np.median(noise)),
                 edge_candidate_count=len(points), edge_signal_source='original_image_gray',
                 edge_background_slope_x=float(background_slope[0]),
                 edge_background_slope_y=float(background_slope[1]),
                 edge_template_forced=False)
    if support < .85 or gap > max(3, int(.08 * len(points))):
        return None, audit
    shifts, measured, cost = _yakun_v4_dynamic_track(candidates, np.zeros(len(points)), search)
    final_points = points + shifts[:, None]*normals
    if not measured.all():
        keep = np.where(measured)[0]
        missing = np.where(~measured)[0]
        for dimension in (0, 1):
            values = final_points[keep, dimension]
            # Bridge actual supported endpoints. Interpolating displacement alone
            # would preserve an unsupported spike in the coarse contour and let
            # that spike incorrectly determine the measured X/Y extrema.
            final_points[missing, dimension] = np.interp(
                missing, np.r_[keep-len(points), keep, keep+len(points)], np.tile(values, 3))
    final = final_points.astype(np.float32).reshape(-1, 1, 2)
    displacement = np.linalg.norm(final_points-points, axis=1)
    original = np.asarray(contour, np.float32).reshape(-1, 1, 2)
    area_ratio = cv2.contourArea(final) / max(cv2.contourArea(original), 1e-6)
    span_ratio = np.ptp(final[:, 0], axis=0) / np.maximum(spans, 1e-6)
    center_shift = float(np.linalg.norm(final[:, 0].mean(axis=0) - points.mean(axis=0)))
    audit.update(edge_area_ratio=float(area_ratio), edge_center_shift_px=center_shift,
                 edge_median_shift_px=float(np.median(displacement)),
                 edge_max_shift_px=float(np.max(displacement)), edge_tracking_cost=float(cost),
                 edge_interpolation_method='linear_between_supported_edge_positions',
                 edge_interpolated_point_count=int(np.count_nonzero(~measured)))
    geometry_ok = (.70 <= area_ratio <= 1.40 and np.all((span_ratio >= .75) & (span_ratio <= 1.25))
                   and center_shift <= max(1.5, .10 * width))
    if not geometry_ok:
        # Keep a V27 boundary only when the original pixels already support it.
        # This is explicit pixel-level retention, not a fitted rectangle/ellipse.
        if float(np.mean(anchor_supported)) >= .92 and _v32_max_cyclic_gap(anchor_supported) <= 3:
            audit.update(edge_status='V27_SUPPORTED_CONTOUR_RETAINED', edge_method='V27_supported_contour',
                         edge_threshold_applied=False)
            return original, audit
        audit['edge_status'] = 'SHAPE_CHANGE_EXCEEDS_LIMIT'
        return None, audit
    audit.update(edge_status='OK', edge_method='Yakun_gray50_normal_tracking', edge_threshold_applied=True)
    return final, audit


def _v32_refine_desc(gray, desc, pattern):
    contour, audit = _v32_supported_contour(_v32_raw_gray(gray), desc.contour)
    if contour is None:
        return None, audit
    refined = describe_contour(contour)
    if refined is None:
        audit['edge_status'] = 'INVALID_REFINED_CONTOUR'
        return None, audit
    for name, value in vars(desc).items():
        if name.startswith('_'):
            setattr(refined, name, value)
    moments = cv2.moments(contour)
    if abs(moments['m00']) > 1e-6:
        refined.center = np.array([moments['m10'] / moments['m00'], moments['m01'] / moments['m00']])
    if pattern == 'slot' and not _slot_shape_ok_v17(refined, (SLOT_MIN_AREA_PX, np.inf), seed=False):
        audit['edge_status'] = 'REFINED_SLOT_SHAPE_QC_FAILED'
        return None, audit
    refined._v32_edge_audit = audit
    return refined, audit


def _v32_edge_reject(desc, pattern, audit):
    return _make_reject_from_contour(desc.contour, f'{pattern}_v32_edge_qc',
                                    'original_image_edge_not_supported', str(audit.get('edge_status')),
                                    desc=desc, **audit)


def _v32_valid_slot_columns(gray, result, locator):
    rows, selected, _, rejects, meta = result
    templates = {int(row['object_id']): row for row in rows}
    columns, failed = {}, set()
    for desc in selected:
        column = int(getattr(desc, '_v23_slot_col', 0))
        row_id = int(getattr(desc, '_v23_slot_row', 0))
        refined, audit = _v32_refine_desc(gray, desc, 'slot')
        if refined is None:
            xp = float(meta.get('slot_horizontal_pitch_px', np.nan))
            yp = float(meta.get('slot_row_pitch_px', np.nan))
            if np.isfinite(xp) and np.isfinite(yp) and min(xp, yp) > 5:
                # Only repair a contour that failed original-image evidence.
                # This avoids V30's unconditional replacement of good V27 shapes.
                bx, by, bw, bh = cv2.boundingRect(desc.contour.astype(np.int32))
                proposal = _v32_recover_flatfield_slot(_v32_raw_gray(gray), tuple(desc.center), xp, yp,
                                                      bw, bh, (SLOT_MIN_AREA_PX, np.inf))
                if (proposal is not None and abs(proposal.center[0]-desc.center[0]) <= max(3., .25*xp)
                        and abs(proposal.center[1]-desc.center[1]) <= max(3., .30*yp)):
                    alternative, alternative_audit = _v32_refine_desc(gray, proposal, 'slot')
                    if alternative is not None:
                        refined, audit = alternative, alternative_audit
                        audit['edge_coarse_source'] = 'V30_repair_of_unsupported_V27_candidate'
        if refined is None:
            failed.add(column)
            rejects.append(_v32_edge_reject(desc, 'slot', audit))
        else:
            columns.setdefault(column, {})[row_id] = (refined, {'slot_row': row_id})
    accepted = []
    for column, members in columns.items():
        if column not in failed and set(members) == {1, 2, 3, 4} and column in templates:
            accepted.append((members, templates[column], locator))
        elif column not in failed:
            failed.add(column)
    meta = dict(meta, v32_edge_rejected_columns=len(failed))
    return accepted, rejects, meta


def _select_and_measure_bottom_slots_v19(gray, pool, px_nm, image_name, stage, condition, pattern_key):
    print(f'      [{image_name}] V32: V27 slot localization + original-image edge validation', flush=True)
    primary = _V32_SLOT_BASE(gray, _copy_v32.deepcopy(pool), px_nm, image_name, stage, condition, pattern_key)
    columns, rejects, meta = _v32_valid_slot_columns(gray, primary, 'V27')
    expected = int(meta.get('array_expected_positions', 0)) // 4
    # V30 supplies missing localization only; a validated V27 column is kept.
    if len(columns) < max(1, expected):
        extra = _v32_slot_flatfield_candidates(gray, pool)
        supplemented, added = merge_candidate_descs_v12(_copy_v32.deepcopy(pool), extra, 'slot')
        meta['v32_v30_locator_added'] = added
        if extra:
            print(f'      [{image_name}] V32: supplemental V30 proposals={len(extra)}; recheck original edges', flush=True)
            secondary = _V32_SLOT_BASE(gray, supplemented, px_nm, image_name, stage, condition, pattern_key)
            candidates, extra_rejects, secondary_meta = _v32_valid_slot_columns(gray, secondary, 'V30_recovery')
            rejects.extend(extra_rejects)
            if not columns and candidates:
                meta.update(secondary_meta)
            pitch = float(meta.get('slot_horizontal_pitch_px', np.nan))
            if not np.isfinite(pitch):
                pitch = float(secondary_meta.get('slot_horizontal_pitch_px', 20.))
            tolerance = max(4., .28 * pitch) if np.isfinite(pitch) else 6.
            for candidate in candidates:
                center = candidate[0][4][0].center
                if not any(np.linalg.norm(center - item[0][4][0].center) <= tolerance for item in columns):
                    columns.append(candidate)
    if not columns:
        meta.update(v23_complete_slot_columns=0, slot_bottom_row_used_count=0)
        return [], [], [], rejects, meta
    axis = _slot_axis_from_columns_v23([item[0] for item in columns], gray.shape)
    xa, yu, origin = (np.asarray(axis[key], float) for key in ('x_axis', 'y_axis_up', 'origin'))
    columns.sort(key=lambda item: float((item[0][4][0].center - origin) @ xa))
    rows, all_selected, bottom = [], [], []
    for column_id, (members, template, locator) in enumerate(columns, 1):
        for row_id in range(1, 5):
            desc = members[row_id][0]
            desc._v23_slot_col, desc._v23_slot_row, desc._v22_slot_axis = column_id, row_id, axis
            all_selected.append(desc)
        desc = members[4][0]
        bottom.append(desc)
        row = dict(template, object_id=column_id, slot_column_id_v23=column_id,
                   slot_complete_column_count=len(columns), v32_locator=locator,
                   v27_x_width_nm=template['x_width_nm'], v27_y_height_nm=template['y_height_nm'])
        row.update(desc._v32_edge_audit)
        row.update(x_width_px=_projection_span_v19(desc.contour, xa), y_height_px=_projection_span_v19(desc.contour, yu),
                   center_image_x_px=float(desc.center[0]), center_image_y_px=float(desc.center[1]),
                   center_axis_x_px=float((desc.center-origin) @ xa), center_axis_y_up_px=float((desc.center-origin) @ yu),
                   slot_pixel_area_px=contour_pixel_count(desc.contour),
                   x_axis_angle_from_image_x_deg=axis['x_axis_angle_from_image_x_deg'],
                   y_axis_down_angle_from_image_x_deg=axis['y_axis_down_angle_from_image_x_deg'])
        row['x_width_nm'], row['y_height_nm'] = row['x_width_px'] * px_nm, row['y_height_px'] * px_nm
        regularity = slot_regularity_metrics(desc)
        row.update(slot_regularity_solidity=regularity['regularity_solidity'],
                   slot_regularity_convexity=regularity['regularity_convexity'],
                   slot_best_template_iou=regularity['regularity_template_iou'],
                   selection_source=_source_of(desc))
        rows.append(row)
    meta.update(v23_complete_slot_columns=len(rows), slot_bottom_row_used_count=len(rows),
                v23_axis_pair_count=axis['pair_count'], v23_axis_fallback=axis['fallback'],
                v23_x_axis_angle_from_image_x_deg=axis['x_axis_angle_from_image_x_deg'],
                v23_y_axis_down_angle_from_image_x_deg=axis['y_axis_down_angle_from_image_x_deg'],
                v32_edge_method='V27 localization / Yakun original-gray normal tracking / V30 missing-site proposals')
    return rows, all_selected, bottom, rejects, meta


def _v32_via_search_scale(pool, image_shape, px_nm, pattern_key):
    """Estimate a search scale from repeated, complete near-round candidates.

    The nominal 60 nm identifies the pattern; it is not a measured-CD ceiling.
    This only seeds localization. Lattice, shape and original-gray edge checks
    still decide which objects are measured. PixelSize never changes the contour.
    """
    nominal_nm, nominal_px = _via_nominal_diam_px_v21(pattern_key, px_nm)
    height, width = image_shape[:2]
    candidates = []
    for desc in pool:
        if desc is None or not (desc.solidity >= .70 and desc.circularity >= .35 and desc.aspect <= 1.8):
            continue
        diameter = 2.0 * math.sqrt(max(0., desc.area) / math.pi)
        if not np.isfinite(diameter) or diameter < 6.0 or not _bbox_complete_v19(desc, image_shape, .10):
            continue
        if desc.area > .10 * height * width:
            continue
        # Candidate pools can contain several thresholds of the same hole.
        if any(np.linalg.norm(desc.center - center) < .35 * min(diameter, old_d)
               for center, old_d in candidates):
            continue
        candidates.append((desc.center, diameter))
    diameters = np.asarray([d for _, d in candidates], dtype=float)
    if len(diameters) >= 4:
        median = float(np.median(diameters))
        support = diameters[(diameters >= .65 * median) & (diameters <= 1.50 * median)]
    else:
        support = np.empty(0)
    if len(support) >= 4:
        search_px, source = float(np.median(support)), 'candidate_pool_median'
    else:
        search_px, source = nominal_px, 'nominal_fallback'
    return search_px, dict(via_search_diameter_px=search_px,
                           via_search_diameter_nm=search_px * px_nm,
                           via_search_scale_source=source,
                           via_search_scale_candidate_count=int(len(support)),
                           via_nominal_diameter_nm=nominal_nm,
                           via_cd_policy='image_scale_relative; no fixed nominal-nm cutoff')


def _select_and_measure_vias_v19(gray, pool, px_nm, image_name, stage, condition, pattern_key):
    search_px, scale_meta = _v32_via_search_scale(pool, gray.shape, px_nm, pattern_key)
    rows, selected, rejects, bundle = _V32_VIA_BASE(
        gray, pool, px_nm, image_name, stage, condition, pattern_key, search_diam_px=search_px)
    bundle['meta'].update(scale_meta)
    for row in rows:
        row.update(scale_meta)
    accepted, output = [], []
    for desc, old in zip(selected, rows):
        refined, audit = _v32_refine_desc(gray, desc, 'via')
        if refined is not None:
            ok, shape, _ = _via_shape_qc_v24(refined.contour, old['via_search_diameter_px'], stage, pattern_key)
            if not ok:
                refined = None
                audit['edge_status'] = 'REFINED_VIA_SHAPE_QC_FAILED'
        if refined is None:
            rejects.append(_v32_edge_reject(desc, 'via', audit))
            continue
        row = dict(old, v27_x_width_nm=old['x_width_nm'], v27_y_height_nm=old['y_height_nm'], **audit)
        row.update(axis_ratio=shape.get('aspect', refined.aspect),
                   via_ellipse_iou=shape.get('ellipse_iou', np.nan), via_radial_cv=shape.get('radial_cv', np.nan))
        accepted.append(refined)
        output.append(row)
    # Refit after rejecting a changed boundary; compare with the image scale, not 60 nm.
    while accepted:
        axis = _fit_via_bottom_row_axes_v19(accepted, gray.shape)
        kept = []
        for desc, row in zip(accepted, output):
            spans = [_projection_span_v19(desc.contour, np.asarray(axis[key]))*px_nm
                     for key in ('x_axis', 'y_axis_up')]
            reference_nm = row['via_search_diameter_nm']
            if all(V23_VIA_MIN_SPAN_FACTOR*reference_nm <= span <= V23_VIA_MAX_SPAN_FACTOR*reference_nm for span in spans):
                kept.append((desc, row))
            else:
                audit = dict(desc._v32_edge_audit, edge_status='REFINED_VIA_PHYSICAL_SPAN_QC_FAILED')
                rejects.append(_v32_edge_reject(desc, 'via', audit))
        if len(kept) == len(accepted):
            break
        accepted, output = [item[0] for item in kept], [item[1] for item in kept]
    axis = _fit_via_bottom_row_axes_v19(accepted, gray.shape)
    xa, yu, origin = (np.asarray(axis[key], float) for key in ('x_axis', 'y_axis_up', 'origin'))
    for index, (desc, row) in enumerate(zip(accepted, output)):
        row.update(object_id=index+1, is_axis_bottom_row=index in axis.get('bottom_indices', []),
                   x_width_px=_projection_span_v19(desc.contour, xa), y_height_px=_projection_span_v19(desc.contour, yu),
                   center_image_x_px=float(desc.center[0]), center_image_y_px=float(desc.center[1]),
                   center_axis_x_px=float((desc.center-origin) @ xa), center_axis_y_up_px=float((desc.center-origin) @ yu),
                   circularity=desc.circularity, solidity=desc.solidity,
                   via_equivalent_diameter_nm=2*math.sqrt(desc.area/math.pi)*px_nm,
                   via_axis_angle_from_image_x_deg=axis.get('row_angle_deg', np.nan),
                   via_axis_fit_rms_px=axis.get('row_fit_rms_px', np.nan),
                   x_edge_method=row['edge_method'], y_edge_method=row['edge_method'])
        row['x_width_nm'], row['y_height_nm'] = row['x_width_px']*px_nm, row['y_height_px']*px_nm
    bundle = dict(bundle, axis=axis, meta=dict(bundle['meta'], v32_edge_validated_count=len(output),
                  via_v24_final_count=len(output), via_bottom_row_count=len(axis.get('bottom_indices', [])),
                  via_axis_angle_deg=axis.get('row_angle_deg', np.nan),
                  via_axis_fit_rms_px=axis.get('row_fit_rms_px', np.nan),
                  via_axis_fallback=bool(axis.get('fallback', False))))
    return output, accepted, rejects, bundle


def _select_and_measure_trenches_v19(gray, pool, px_nm, image_name, stage, condition, pattern_key):
    rows, pairs, rejects, meta = _V32_TRENCH_BASE(gray, pool, px_nm, image_name, stage, condition, pattern_key)
    raw = _v32_raw_gray(gray).astype(float)
    if FEATURE_POLARITY.lower() == 'bright':
        raw = 255. - raw
    accepted, output = [], []
    for original, old in zip(pairs, rows):
        pair = dict(original)
        widths = []
        for prefix, body_range in (('top', V20_TRENCH_TOP_BODY_RANGE), ('bottom', V20_TRENCH_BOTTOM_BODY_RANGE)):
            start, end = [int(round(value * raw.shape[0])) for value in body_range]
            profile = np.mean(raw[start:end], axis=0)
            left, right = pair[f'{prefix}_x_left_px'], pair[f'{prefix}_x_right_px']
            edge = _yakun_v111_edge_pair(profile, left, right)
            if not edge.get('valid') or abs(edge['width']-(right-left)) > max(3., .20*(right-left)):
                break
            pair[f'{prefix}_x_left_px'], pair[f'{prefix}_x_right_px'] = edge['left'], edge['right']
            pair[f'{prefix}_width_px_v22'] = edge['width']
            widths.append(edge['width'])
        if len(widths) == 2:
            x0 = max(0, int(math.ceil(max(pair['top_x_left_px'], pair['bottom_x_left_px']))))
            x1 = min(raw.shape[1], int(math.floor(min(pair['top_x_right_px'], pair['bottom_x_right_px']))))
            if x1 > x0:
                profile = np.mean(raw[:, x0:x1], axis=1)
                scale = max(8., pair['bottom_tip_y_px'] - pair['top_tip_y_px'])
                upper = _yakun_v111_one_sided_dark_edge(profile, pair['top_tip_y_px'], 'right', scale)
                lower = _yakun_v111_one_sided_dark_edge(profile, pair['bottom_tip_y_px'], 'left', scale)
                if upper.get('valid') and lower.get('valid') and upper['edge'] < lower['edge']:
                    shifts = [abs(upper['edge']-pair['top_tip_y_px']), abs(lower['edge']-pair['bottom_tip_y_px'])]
                    if max(shifts) <= max(3., .15 * scale):
                        pair['top_tip_y_px'], pair['bottom_tip_y_px'] = upper['edge'], lower['edge']
                        top = np.array([.5*(pair['top_x_left_px']+pair['top_x_right_px']), upper['edge']])
                        bottom = np.array([.5*(pair['bottom_x_left_px']+pair['bottom_x_right_px']), lower['edge']])
                        axis = pair['_v22_axis']; mid = .5*(top+bottom)
                        gap = abs(float((bottom-top) @ np.asarray(axis['y_axis_down'])))
                        width = .5 * sum(widths)
                        origin = np.asarray(axis['origin'])
                        axis_x = float((mid-origin) @ np.asarray(axis['x_axis']))
                        axis_y = float((mid-origin) @ np.asarray(axis['y_axis_up']))
                        pair.update(x_width_axis_px=width, tip_gap_axis_px=gap,
                                    axis_center_x_px=axis_x, axis_center_y_up_px=axis_y,
                                    top_width_px_v22=widths[0], bottom_width_px_v22=widths[1])
                        row = dict(old, object_id=len(output)+1, v27_x_width_nm=old['x_width_nm'],
                                   v27_tip_gap_y_nm=old['tip_gap_y_nm'], x_width_px=width, x_width_nm=width*px_nm,
                                   tip_gap_y_nm=gap*px_nm, trench_mean_rectangle_width_px=width,
                                   trench_mean_rectangle_width_nm=width*px_nm,
                                   top_trench_width_nm=widths[0]*px_nm, bottom_trench_width_nm=widths[1]*px_nm,
                                   center_image_x_px=float(mid[0]), center_image_y_px=float(mid[1]),
                                   center_axis_x_px=axis_x, center_axis_y_up_px=axis_y,
                                   width_mismatch_frac=abs(widths[0]-widths[1])/max(width, 1e-9),
                                   trench_detector='V27_geometry+Yakun_anchored_gradient_sideband50',
                                   edge_status='OK', edge_method='Yakun_anchored_gradient_sideband50',
                                   edge_threshold_percent=50., edge_threshold_applied=True,
                                   edge_signal_source='original_image_gray', edge_template_forced=False)
                        accepted.append(pair); output.append(row)
                        continue
        # Do not turn a synthetic geometric rectangle into a successful measurement.
        c = np.array([[original['top_x_left_px'], original['top_tip_y_px']-8],
                      [original['top_x_right_px'], original['top_tip_y_px']-8],
                      [original['top_x_right_px'], original['top_tip_y_px']],
                      [original['top_x_left_px'], original['top_tip_y_px']]], dtype=np.int32).reshape(-1, 1, 2)
        rejects.append(_make_reject_from_contour(c, 'trench_v32_edge_qc', 'original_image_edge_not_supported',
                                                 'Yakun width/tip signal checks failed or changed V27 geometry too far'))
    return output, accepted, rejects, dict(meta, v32_edge_validated_count=len(output))


# ============================================================================
# V32: MULTIFORMAT PIPELINE -- SLOT RECOVERY AND HARD PER-IMAGE TIMEOUT
# ============================================================================
# V27 geometry and output definitions, Yakun gray edges, V30 bounded recovery.
# No sem_cd_measure_200k_batch_V1_6.py code, import, or measurement rule is used.
#
# Dataset mapping per condition folder (natural filename order):
#   1-3 : Region 2 -> trench, slot, via
#   4-6 : Region 3 -> trench, slot, via
#   7-9 : Region 4 -> trench, slot, via
#   10+ : continue the same three-pattern order in Regions 5, 6, ...
# Any image count is supported, including an incomplete final group, provided
# the corresponding before/after condition folders contain the same count.
#
# Each supported image must have a same-stem TXT containing:
#   PixelSize=<number>
#
# Before/after pairing is by condition + region + pattern, NOT by equal filename.

import argparse as _argparse_v32
import json as _json_v32
import traceback as _traceback_v32
import pickle as _pickle_v32
import shutil as _shutil_v32
import subprocess as _subprocess_v32
import tempfile as _tempfile_v32
from typing import Any as _Any_v32, Iterable as _Iterable_v32, Sequence as _Sequence_v32

try:
    from scipy.optimize import linear_sum_assignment as _linear_sum_assignment_v32
except ImportError as _exc_v32:
    raise SystemExit(
        "V32 需要 scipy（逐对象 before/after 配对使用）。请运行：\n"
        "  pip install numpy pandas scipy opencv-python matplotlib openpyxl\n"
        f"原始错误：{_exc_v32}"
    ) from _exc_v32


V32_SCRIPT_VERSION = "2026-09-24-V32-VIA60-IMAGE-SCALE-RELAXED-CD"
V32_DEFAULT_BEFORE_ROOT = Path(r"C:\Users\z00027644\Documents\倾斜刻蚀\SEM\830SEM_before_treat")
V32_DEFAULT_AFTER_ROOT = Path(r"C:\Users\z00027644\Documents\倾斜刻蚀\SEM\0831SEM_10-14_topview_after")
V32_FIRST_REGION = 2
V32_PATTERN_SEQUENCE = ("trench", "slot", "via")
V32_STAGES = ("before", "after")
V32_STAGE_ZH = {"before": "处理前", "after": "处理后"}
V32_PATTERN_ZH = {"trench": "Trench", "slot": "Slot最下排", "via": "Via"}
V32_IMAGE_SUFFIXES = {".tif", ".tiff", ".png", ".jpg", ".jpeg"}

# Historical internal pattern keys. Trench/slot are fixed. Via defaults to via60;
# auto model selection is used only when explicitly requested with --via-pattern.
V32_FIXED_PATTERN_KEY = {
    "trench": "trench160",
    "slot": "slot210",
}
V32_VIA_KEYS = ("via40", "via60")
V32_VIA_NOMINAL_NM = {"via40": 40.0, "via60": 60.0}

# Requested output metrics.
V32_METRIC_SPECS = {
    "trench": [
        ("x_width_nm", "trench_width_nm", "Trench宽度"),
        ("tip_gap_y_nm", "tip_to_tip_y_nm", "Tip-to-tip Y距离"),
    ],
    "slot": [
        ("x_width_nm", "slot_x_width_nm", "Slot X宽度"),
        ("y_height_nm", "slot_y_length_nm", "Slot Y长度"),
    ],
    "via": [
        ("x_width_nm", "via_x_width_nm", "Via X宽度"),
        ("y_height_nm", "via_y_height_nm", "Via Y高度"),
    ],
}

# Per-object normalized-position matching thresholds. This matching is ONLY for
# comparison tables; it does not change any V25 detection or measurement result.
V32_MATCH_THRESHOLD_VIA = 0.30
V32_MATCH_THRESHOLD_1D = 0.28


# -----------------------------------------------------------------------------
# V32 generic helpers
# -----------------------------------------------------------------------------
def _v32_condition_key(condition: str) -> tuple:
    name = str(condition)
    return int(name), name


def _v32_discover_conditions(after_root: Path) -> List[str]:
    """Use exact digit-only folder names under AFTER as the condition list."""
    after_root = Path(after_root)
    if not after_root.exists():
        raise FileNotFoundError(f"after 根目录不存在：{after_root}")
    if not after_root.is_dir():
        raise NotADirectoryError(f"after 路径不是文件夹：{after_root}")
    conditions = sorted(
        [p.name for p in after_root.iterdir() if p.is_dir() and re.fullmatch(r"[0-9]+", p.name)],
        key=_v32_condition_key,
    )
    if not conditions:
        raise ValueError(
            f"after 根目录 {after_root} 下没有以纯数字命名的子文件夹"
            "（例如 0、1、02、100）；仅扫描直接子文件夹。"
        )
    return conditions


def _v32_sort_conditions(df: pd.DataFrame, columns: _Sequence_v32[str]) -> pd.DataFrame:
    """Sort condition names numerically while retaining exact names, e.g. 007."""
    if df.empty:
        return df
    out = df.copy()
    names = sorted(out["condition"].astype(str).unique(), key=_v32_condition_key)
    order = {name: index for index, name in enumerate(names)}
    out["_condition_order_v32"] = out["condition"].astype(str).map(order)
    sort_columns = ["_condition_order_v32"] + [c for c in columns if c != "condition"]
    return out.sort_values(sort_columns, kind="stable").drop(columns=["_condition_order_v32"]).reset_index(drop=True)


def _v32_finite_float(value: _Any_v32, default: float = math.nan) -> float:
    try:
        x = float(value)
    except Exception:
        return default
    return x if np.isfinite(x) else default


def _v32_safe_mean(values: _Iterable_v32) -> float:
    a = pd.to_numeric(pd.Series(list(values), dtype="object"), errors="coerce").dropna().to_numpy(float)
    return float(np.mean(a)) if len(a) else math.nan


def _v32_safe_median(values: _Iterable_v32) -> float:
    a = pd.to_numeric(pd.Series(list(values), dtype="object"), errors="coerce").dropna().to_numpy(float)
    return float(np.median(a)) if len(a) else math.nan


def _v32_natural_key(path: Path) -> tuple:
    parts = re.split(r"(\d+)", path.stem.lower())
    return (tuple(int(p) if p.isdigit() else p for p in parts), path.name.lower(), path.name)


def _v32_df(rows) -> pd.DataFrame:
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def _v32_write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8-sig")


def _v32_expected_positions(image_count: int) -> List[dict]:
    """Extend the original trench/slot/via order to every discovered image."""
    return [
        {
            "sequence_index": index + 1,
            "region": V32_FIRST_REGION + index // len(V32_PATTERN_SEQUENCE),
            "pattern": V32_PATTERN_SEQUENCE[index % len(V32_PATTERN_SEQUENCE)],
        }
        for index in range(image_count)
    ]


# -----------------------------------------------------------------------------
# Inventory / natural-order mapping with matching before/after counts
# -----------------------------------------------------------------------------
def _v32_discover_one_condition(stage: str, root: Path, condition: str):
    folder = Path(root) / str(condition)
    rows: List[dict] = []
    errors: List[dict] = []

    if not folder.exists() or not folder.is_dir():
        errors.append({
            "stage": stage,
            "condition": str(condition),
            "error_type": "MISSING_CONDITION_FOLDER",
            "message": f"找不到目标子文件夹：{folder}",
        })
        return rows, errors

    images = sorted(
        [p for p in folder.iterdir() if p.is_file() and p.suffix.lower() in V32_IMAGE_SUFFIXES],
        key=_v32_natural_key,
    )

    for pos, image_path in zip(_v32_expected_positions(len(images)), images):
        # Validate same-stem TXT + PixelSize now, but leave actual parsing to the exact V25
        # read_pixel_size_nm() again during measurement.
        txt_path = image_path.with_suffix(".txt")
        pixel_size_nm = math.nan
        pixel_note = ""
        try:
            info = read_pixel_size_nm(image_path)
            pixel_size_nm = float(info.value_nm)
        except Exception as exc:
            pixel_note = f"{type(exc).__name__}: {exc}"
            errors.append({
                "stage": stage,
                "condition": str(condition),
                "region": int(pos["region"]),
                "pattern": str(pos["pattern"]),
                "image": image_path.name,
                "path": str(image_path),
                "error_type": "PIXELSIZE_PREFLIGHT_ERROR",
                "message": pixel_note,
            })

        rows.append({
            "stage": stage,
            "stage_zh": V32_STAGE_ZH[stage],
            "condition": str(condition),
            "condition_label": f"Condition {condition}",
            "region": int(pos["region"]),
            "pattern": str(pos["pattern"]),
            "pattern_zh": V32_PATTERN_ZH[str(pos["pattern"])],
            "sequence_index": int(pos["sequence_index"]),
            "image": image_path.name,
            "path": str(image_path),
            "txt_path_expected": str(txt_path),
            "pixel_size_nm_preflight": pixel_size_nm,
            "pixel_size_preflight_note": pixel_note,
            "folder": str(folder),
            "pair_key": f"C{condition}_R{pos['region']}_{pos['pattern']}",
        })
    return rows, errors


def _v32_build_inventory(
    before_root: Path, after_root: Path, conditions: Optional[_Sequence_v32[str]] = None,
):
    if conditions is None:
        conditions = _v32_discover_conditions(after_root)
    rows: List[dict] = []
    errors: List[dict] = []
    for condition in conditions:
        condition = str(condition)
        before_rows, before_errors = _v32_discover_one_condition("before", before_root, condition)
        after_rows, after_errors = _v32_discover_one_condition("after", after_root, condition)
        condition_errors = before_errors + after_errors
        errors.extend(condition_errors)
        if any(e["error_type"] == "MISSING_CONDITION_FOLDER" for e in condition_errors):
            continue
        if len(before_rows) != len(after_rows):
            errors.append({
                "stage": "before/after",
                "condition": condition,
                "error_type": "IMAGE_COUNT_MISMATCH",
                "before_count": len(before_rows),
                "after_count": len(after_rows),
                "before_folder": str(Path(before_root) / condition),
                "after_folder": str(Path(after_root) / condition),
                "message": (
                    f"Condition {condition} 的图片数量不一致："
                    f"before={len(before_rows)}，after={len(after_rows)}。"
                    "跳过该 condition 的前后两组照片，避免错位配对。"
                ),
            })
            continue
        rows.extend(before_rows)
        rows.extend(after_rows)
    inv = pd.DataFrame(rows)
    err = pd.DataFrame(errors)
    if not inv.empty:
        inv = _v32_sort_conditions(inv, ["condition", "region", "sequence_index", "stage"])
    return inv, err


# -----------------------------------------------------------------------------
# Exact V25 processing with an externally assigned pattern key
# -----------------------------------------------------------------------------
def _v32_process_v25_assigned(
    image_path: Path,
    stage: str,
    condition: str,
    region: int,
    public_pattern: str,
    v25_pattern_key: str,
    ann_dir: Path,
):
    """Assign pattern keys to the embedded pipeline, including V32 slot recovery."""
    if v25_pattern_key not in PATTERN_KEYS_V19:
        raise ValueError(f"非法 V25 pattern key: {v25_pattern_key}")

    global _classify_pattern_key_v19
    original_classifier = _classify_pattern_key_v19
    target = Path(image_path).resolve()

    def assigned_classifier(p: Path):
        try:
            if Path(p).resolve() == target:
                return v25_pattern_key
        except Exception:
            pass
        return original_classifier(p)

    _classify_pattern_key_v19 = assigned_classifier
    try:
        rows, status, rejects = _process_one_image_v19(
            Path(image_path), stage, str(condition), Path(ann_dir)
        )
    finally:
        _classify_pattern_key_v19 = original_classifier

    pair_key = f"C{condition}_R{region}_{public_pattern}"
    normalized_rows: List[dict] = []
    for r0 in rows:
        r = dict(r0)
        r["source_pattern_key"] = r.get("pattern", v25_pattern_key)
        r["pattern"] = public_pattern
        r["pattern_zh"] = V32_PATTERN_ZH[public_pattern]
        r["region"] = int(region)
        r["pair_key"] = pair_key
        r["source_algorithm"] = "V32_V27_Yakun_fusion"
        r["x_width_nm"] = _v32_finite_float(r.get("x_width_nm"))
        r["y_height_nm"] = _v32_finite_float(r.get("y_height_nm"))
        r["tip_gap_y_nm"] = _v32_finite_float(r.get("tip_gap_y_nm"))
        normalized_rows.append(r)

    status = dict(status)
    status["source_pattern_key"] = status.get("pattern", v25_pattern_key)
    status["pattern"] = public_pattern
    status["pattern_zh"] = V32_PATTERN_ZH[public_pattern]
    status["region"] = int(region)
    status["pair_key"] = pair_key
    status["source_algorithm"] = "V32_V27_Yakun_fusion"

    normalized_rejects: List[dict] = []
    for q0 in rejects:
        q = dict(q0)
        q["source_pattern_key"] = q.get("pattern", v25_pattern_key)
        q["pattern"] = public_pattern
        q["region"] = int(region)
        q["pair_key"] = pair_key
        q["source_algorithm"] = "V32_V27_Yakun_fusion"
        normalized_rejects.append(q)

    return normalized_rows, status, normalized_rejects


# -----------------------------------------------------------------------------
# Optional Via40 / Via60 model choice using V32
# -----------------------------------------------------------------------------
def _v32_via_trial_quality(rows: List[dict], status: dict, key: str) -> tuple:
    """Lexicographic quality used only for explicit auto mode on BEFORE.

    Rows have passed V32 original-gray edge validation. Higher is better.
    """
    if not rows:
        return (-1, -1, -1e9, -1e9, -1e9, -1e9, -1e9)

    df = pd.DataFrame(rows)
    n = int(len(df))
    fallback_level = int(_v32_finite_float(status.get("auto_fallback_level_used", 0), 0.0))
    strict = 1 if fallback_level == 0 else 0
    occupancy = _v32_finite_float(status.get("measure_array_occupancy_seed"), -1.0)
    good_ray = _v32_safe_mean(df.get("via_good_ray_fraction", pd.Series(dtype=float)).tolist())
    scale_score = _v32_safe_mean(df.get("via_scale_match_score", pd.Series(dtype=float)).tolist())
    leak = _v32_safe_mean(df.get("via_outer_dark_leak_fraction", pd.Series(dtype=float)).tolist())

    x = pd.to_numeric(df.get("x_width_nm", pd.Series(dtype=float)), errors="coerce")
    y = pd.to_numeric(df.get("y_height_nm", pd.Series(dtype=float)), errors="coerce")
    d = (0.5 * (x + y)).dropna().to_numpy(float)
    med = float(np.median(d)) if len(d) else math.nan
    nominal = V32_VIA_NOMINAL_NM[key]
    nominal_rel_error = abs(med - nominal) / nominal if np.isfinite(med) else 99.0

    # First ensure a repeated valid array, then favor strict detection, lattice
    # consistency, number of accepted vias, edge quality, and physical-size consistency.
    return (
        1 if n >= 3 else 0,
        strict,
        occupancy,
        float(n),
        good_ray if np.isfinite(good_ray) else -1.0,
        scale_score if np.isfinite(scale_score) else -1.0,
        -(leak if np.isfinite(leak) else 1.0) - nominal_rel_error,
    )


def _v32_choose_via_key_from_before(
    before_path: Path,
    condition: str,
    region: int,
    trial_root: Path,
    forced: str,
):
    if forced in V32_VIA_KEYS:
        return forced, [{
            "condition": str(condition),
            "region": int(region),
            "before_image": before_path.name,
            "trial_pattern_key": forced,
            "mode": "FORCED",
            "selected": True,
            "n_measurements": math.nan,
            "quality_key": "forced",
        }]

    trial_rows: List[dict] = []
    best_key = "via40"
    best_quality = (-1, -1, -1e9, -1e9, -1e9, -1e9, -1e9)
    for key in V32_VIA_KEYS:
        ann_dir = trial_root / str(condition) / f"Region_{region}" / key
        try:
            rows, status, _ = _v32_process_v25_assigned(
                before_path, "before", str(condition), int(region), "via", key, ann_dir
            )
            quality = _v32_via_trial_quality(rows, status, key)
            err = ""
            n = len(rows)
            fallback = int(_v32_finite_float(status.get("auto_fallback_level_used", 0), 0.0))
            medx = _v32_safe_median(pd.DataFrame(rows).get("x_width_nm", pd.Series(dtype=float)).tolist()) if rows else math.nan
            medy = _v32_safe_median(pd.DataFrame(rows).get("y_height_nm", pd.Series(dtype=float)).tolist()) if rows else math.nan
        except Exception as exc:
            quality = (-1, -1, -1e9, -1e9, -1e9, -1e9, -1e9)
            err = f"{type(exc).__name__}: {exc}"
            n = 0
            fallback = math.nan
            medx = math.nan
            medy = math.nan

        trial_rows.append({
            "condition": str(condition),
            "region": int(region),
            "before_image": before_path.name,
            "trial_pattern_key": key,
            "trial_nominal_nm": V32_VIA_NOMINAL_NM[key],
            "mode": "AUTO_V32",
            "selected": False,
            "n_measurements": n,
            "fallback_level": fallback,
            "median_x_width_nm": medx,
            "median_y_height_nm": medy,
            "quality_key": repr(quality),
            "error": err,
        })
        if quality > best_quality:
            best_quality = quality
            best_key = key

    for r in trial_rows:
        r["selected"] = bool(r["trial_pattern_key"] == best_key)
    return best_key, trial_rows


# -----------------------------------------------------------------------------
# Image-level statistics and before/after comparison
# -----------------------------------------------------------------------------
def _v32_build_image_metrics(objects_df: pd.DataFrame) -> pd.DataFrame:
    if objects_df.empty:
        return pd.DataFrame()
    records: List[dict] = []
    group_cols = [
        "stage", "stage_zh", "condition", "condition_label", "region",
        "pair_key", "pattern", "pattern_zh", "image",
    ]
    for keys, g in objects_df.groupby(group_cols, dropna=False, sort=True):
        base = dict(zip(group_cols, keys if isinstance(keys, tuple) else (keys,)))
        pattern = str(base["pattern"])
        for source_col, metric, metric_zh in V32_METRIC_SPECS[pattern]:
            vals = pd.to_numeric(g.get(source_col, pd.Series(dtype=float)), errors="coerce").dropna().to_numpy(float)
            if not len(vals):
                continue
            records.append({
                **base,
                "metric": metric,
                "metric_zh": metric_zh,
                "source_column": source_col,
                "n": int(len(vals)),
                "mean_nm": float(np.mean(vals)),
                "std_nm": float(np.std(vals, ddof=1)) if len(vals) >= 2 else math.nan,
                "sem_nm": float(np.std(vals, ddof=1) / math.sqrt(len(vals))) if len(vals) >= 2 else math.nan,
                "median_nm": float(np.median(vals)),
                "min_nm": float(np.min(vals)),
                "max_nm": float(np.max(vals)),
            })
    return _v32_sort_conditions(pd.DataFrame(records), ["condition", "region", "pattern", "stage", "image", "metric"])


def _v32_build_image_pair_comparison(image_metrics: pd.DataFrame) -> pd.DataFrame:
    if image_metrics.empty:
        return pd.DataFrame()
    idx = [
        "condition", "condition_label", "region", "pair_key", "pattern",
        "pattern_zh", "metric", "metric_zh", "source_column",
    ]
    out = image_metrics[idx].drop_duplicates().copy()
    for stage in V32_STAGES:
        s = image_metrics[image_metrics["stage"] == stage].copy()
        keep = idx + ["image", "n", "mean_nm", "std_nm", "sem_nm", "median_nm", "min_nm", "max_nm"]
        s = s[keep].rename(columns={
            "image": f"{stage}_image",
            "n": f"{stage}_n",
            "mean_nm": f"{stage}_mean_nm",
            "std_nm": f"{stage}_std_nm",
            "sem_nm": f"{stage}_sem_nm",
            "median_nm": f"{stage}_median_nm",
            "min_nm": f"{stage}_min_nm",
            "max_nm": f"{stage}_max_nm",
        })
        out = out.merge(s, on=idx, how="left")

    before = pd.to_numeric(out.get("before_mean_nm"), errors="coerce")
    after = pd.to_numeric(out.get("after_mean_nm"), errors="coerce")
    out["delta_after_minus_before_nm"] = after - before
    out["abs_delta_nm"] = np.abs(out["delta_after_minus_before_nm"])
    out["change_percent"] = np.where(
        np.isfinite(before) & (np.abs(before) > 1e-12),
        100.0 * out["delta_after_minus_before_nm"] / before,
        np.nan,
    )

    porder = {p: i for i, p in enumerate(V32_PATTERN_SEQUENCE)}
    out["_po"] = out["pattern"].map(porder).fillna(999)
    out = _v32_sort_conditions(out, ["condition", "region", "_po", "metric"]).drop(columns=["_po"])
    return out


def _v32_build_overall_comparison(objects_df: pd.DataFrame) -> pd.DataFrame:
    """Pool all discovered regions per condition, separate from image-by-image comparison."""
    if objects_df.empty:
        return pd.DataFrame()
    long_rows: List[dict] = []
    for (stage, condition, pattern), g in objects_df.groupby(["stage", "condition", "pattern"], sort=True):
        for source_col, metric, metric_zh in V32_METRIC_SPECS[str(pattern)]:
            vals = pd.to_numeric(g.get(source_col, pd.Series(dtype=float)), errors="coerce").dropna().to_numpy(float)
            if not len(vals):
                continue
            long_rows.append({
                "stage": stage,
                "condition": str(condition),
                "condition_label": f"Condition {condition}",
                "pattern": pattern,
                "pattern_zh": V32_PATTERN_ZH[str(pattern)],
                "metric": metric,
                "metric_zh": metric_zh,
                "source_column": source_col,
                "n": int(len(vals)),
                "mean_nm": float(np.mean(vals)),
                "std_nm": float(np.std(vals, ddof=1)) if len(vals) >= 2 else math.nan,
                "median_nm": float(np.median(vals)),
            })
    long = pd.DataFrame(long_rows)
    if long.empty:
        return long
    idx = ["condition", "condition_label", "pattern", "pattern_zh", "metric", "metric_zh", "source_column"]
    out = long[idx].drop_duplicates().copy()
    for stage in V32_STAGES:
        s = long[long["stage"] == stage][idx + ["n", "mean_nm", "std_nm", "median_nm"]].rename(columns={
            "n": f"{stage}_n",
            "mean_nm": f"{stage}_mean_nm",
            "std_nm": f"{stage}_std_nm",
            "median_nm": f"{stage}_median_nm",
        })
        out = out.merge(s, on=idx, how="left")
    b = pd.to_numeric(out.get("before_mean_nm"), errors="coerce")
    a = pd.to_numeric(out.get("after_mean_nm"), errors="coerce")
    out["delta_after_minus_before_nm"] = a - b
    out["change_percent"] = np.where(np.isfinite(b) & (np.abs(b) > 1e-12), 100.0 * (a - b) / b, np.nan)
    return _v32_sort_conditions(out, ["condition", "pattern", "metric"])


# -----------------------------------------------------------------------------
# Per-object spatial before/after matching (comparison only)
# -----------------------------------------------------------------------------
def _v32_normalized_1d(g: pd.DataFrame) -> np.ndarray:
    if "center_axis_x_px" in g:
        s = g["center_axis_x_px"]
    elif "center_image_x_px" in g:
        s = g["center_image_x_px"]
    else:
        return np.linspace(0.0, 1.0, len(g)) if len(g) else np.empty(0)
    x = pd.to_numeric(s, errors="coerce").to_numpy(float)
    if len(x) <= 1:
        return np.full(len(x), 0.5, float)
    if not np.all(np.isfinite(x)) or np.nanmax(x) - np.nanmin(x) < 1e-9:
        return np.linspace(0.0, 1.0, len(x))
    return (x - np.nanmin(x)) / max(np.nanmax(x) - np.nanmin(x), 1e-9)


def _v32_normalized_2d(g: pd.DataFrame) -> np.ndarray:
    out = np.zeros((len(g), 2), float)
    for j, name in enumerate(("center_image_x_px", "center_image_y_px")):
        if name not in g:
            out[:, j] = 0.5
            continue
        arr = pd.to_numeric(g[name], errors="coerce").to_numpy(float)
        if len(arr) <= 1 or not np.all(np.isfinite(arr)) or np.nanmax(arr) - np.nanmin(arr) < 1e-9:
            out[:, j] = 0.5
        else:
            out[:, j] = (arr - np.nanmin(arr)) / max(np.nanmax(arr) - np.nanmin(arr), 1e-9)
    return out


def _v32_match_objects(before: pd.DataFrame, after: pd.DataFrame, pattern: str):
    b = before.copy().reset_index(drop=True)
    a = after.copy().reset_index(drop=True)
    if b.empty:
        return [(None, a.iloc[j], math.nan) for j in range(len(a))]
    if a.empty:
        return [(b.iloc[i], None, math.nan) for i in range(len(b))]

    if pattern == "via":
        pb, pa = _v32_normalized_2d(b), _v32_normalized_2d(a)
        threshold = V32_MATCH_THRESHOLD_VIA
    else:
        pb = _v32_normalized_1d(b)[:, None]
        pa = _v32_normalized_1d(a)[:, None]
        threshold = V32_MATCH_THRESHOLD_1D

    cost = np.linalg.norm(pb[:, None, :] - pa[None, :, :], axis=2)
    bi, aj = _linear_sum_assignment_v32(cost)
    accepted = {}
    used_a = set()
    dist = {}
    for i, j in zip(bi.tolist(), aj.tolist()):
        d = float(cost[i, j])
        if d <= threshold:
            accepted[i] = j
            used_a.add(j)
            dist[(i, j)] = d

    pairs = []
    for i in range(len(b)):
        if i in accepted:
            j = accepted[i]
            pairs.append((b.iloc[i], a.iloc[j], dist[(i, j)]))
        else:
            pairs.append((b.iloc[i], None, math.nan))
    for j in range(len(a)):
        if j not in used_a:
            pairs.append((None, a.iloc[j], math.nan))
    return pairs


def _v32_build_object_comparison(objects_df: pd.DataFrame) -> pd.DataFrame:
    if objects_df.empty:
        return pd.DataFrame()
    records: List[dict] = []
    keys = objects_df[["condition", "region", "pair_key", "pattern"]].drop_duplicates()
    for _, k in keys.iterrows():
        condition = str(k["condition"])
        region = int(k["region"])
        pair_key = str(k["pair_key"])
        pattern = str(k["pattern"])
        g = objects_df[
            (objects_df["condition"].astype(str) == condition)
            & (objects_df["region"] == region)
            & (objects_df["pair_key"] == pair_key)
            & (objects_df["pattern"] == pattern)
        ]
        before = g[g["stage"] == "before"]
        after = g[g["stage"] == "after"]
        pairs = _v32_match_objects(before, after, pattern)
        for pair_id, (b, a, match_distance) in enumerate(pairs, 1):
            base = {
                "condition": condition,
                "condition_label": f"Condition {condition}",
                "region": region,
                "pair_key": pair_key,
                "pattern": pattern,
                "pattern_zh": V32_PATTERN_ZH[pattern],
                "matched_pair_id": pair_id,
                "match_distance_normalized": match_distance,
                "before_object_id": b.get("object_id") if b is not None else math.nan,
                "after_object_id": a.get("object_id") if a is not None else math.nan,
                "before_image": b.get("image") if b is not None else "",
                "after_image": a.get("image") if a is not None else "",
                "match_status": "MATCHED" if b is not None and a is not None else ("BEFORE_ONLY" if b is not None else "AFTER_ONLY"),
            }
            for source_col, metric, metric_zh in V32_METRIC_SPECS[pattern]:
                bv = _v32_finite_float(b.get(source_col)) if b is not None else math.nan
                av = _v32_finite_float(a.get(source_col)) if a is not None else math.nan
                delta = av - bv if np.isfinite(av) and np.isfinite(bv) else math.nan
                rec = dict(base)
                rec.update({
                    "metric": metric,
                    "metric_zh": metric_zh,
                    "before_nm": bv,
                    "after_nm": av,
                    "delta_after_minus_before_nm": delta,
                    "change_percent": 100.0 * delta / bv if np.isfinite(delta) and np.isfinite(bv) and abs(bv) > 1e-12 else math.nan,
                })
                records.append(rec)
    return _v32_sort_conditions(pd.DataFrame(records), ["condition", "region", "pattern", "matched_pair_id", "metric"])


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def _v32_safe_plot_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s))


def _v32_save_pair_plots(pair_df: pd.DataFrame, plot_dir: Path) -> None:
    if pair_df.empty:
        return
    plot_dir.mkdir(parents=True, exist_ok=True)
    for (condition, pattern, metric), g in pair_df.groupby(["condition", "pattern", "metric"], sort=True):
        g = g.sort_values("region")
        x = np.arange(len(g), dtype=float)
        width = 0.36
        b = pd.to_numeric(g["before_mean_nm"], errors="coerce").to_numpy(float)
        a = pd.to_numeric(g["after_mean_nm"], errors="coerce").to_numpy(float)
        be = pd.to_numeric(g.get("before_std_nm"), errors="coerce").to_numpy(float)
        ae = pd.to_numeric(g.get("after_std_nm"), errors="coerce").to_numpy(float)
        be = np.where(np.isfinite(be), be, 0.0)
        ae = np.where(np.isfinite(ae), ae, 0.0)

        fig, ax = plt.subplots(figsize=(8.2, 5.2))
        ax.bar(x - width / 2, b, width, yerr=be, capsize=3, label="Before")
        ax.bar(x + width / 2, a, width, yerr=ae, capsize=3, label="After")
        ax.set_xticks(x, [f"Region {int(r)}" for r in g["region"]])
        ax.set_ylabel("nm")
        ax.set_title(f"Condition {condition} | {str(pattern).upper()} | {metric}")
        ax.legend()
        ax.grid(axis="y", alpha=0.25)
        fig.tight_layout()
        out = plot_dir / f"C{condition}__{pattern}__{_v32_safe_plot_name(metric)}__before_after.png"
        fig.savefig(out, dpi=220)
        plt.close(fig)

    # One delta plot per metric across all conditions/regions.
    for (pattern, metric), g in pair_df.groupby(["pattern", "metric"], sort=True):
        g = g.copy()
        g["label"] = g.apply(lambda r: f"C{r['condition']}-R{int(r['region'])}", axis=1)
        g = _v32_sort_conditions(g, ["condition", "region"])
        d = pd.to_numeric(g["delta_after_minus_before_nm"], errors="coerce").to_numpy(float)
        x = np.arange(len(g), dtype=float)
        fig, ax = plt.subplots(figsize=(max(9.0, 0.55 * len(g)), 5.2))
        ax.bar(x, d)
        ax.axhline(0.0, linewidth=1.0)
        ax.set_xticks(x, g["label"].tolist(), rotation=55, ha="right")
        ax.set_ylabel("After - Before (nm)")
        ax.set_title(f"{str(pattern).upper()} | {metric} | Before/After delta")
        ax.grid(axis="y", alpha=0.25)
        fig.tight_layout()
        out = plot_dir / f"ALL__{pattern}__{_v32_safe_plot_name(metric)}__delta.png"
        fig.savefig(out, dpi=220)
        plt.close(fig)


# -----------------------------------------------------------------------------
# Main V32 pipeline
# -----------------------------------------------------------------------------
class V32ImageTimeout(TimeoutError):
    """The image's complete processing budget, including via trials, expired."""


def _v32_measure_image_job(job: dict):
    """Run inside one disposable worker; return only exportable data."""
    image_path = Path(job['image_path'])
    ann_dir = Path(job['ann_dir'])
    stage, pattern = job['stage'], job['public_pattern']
    key = job.get('v25_pattern_key')
    trials = []
    if pattern == 'via' and stage == 'before' and not key:
        best = None
        for candidate_key in V32_VIA_KEYS:
            try:
                payload = _v32_process_v25_assigned(
                    image_path, stage, job['condition'], job['region'], pattern,
                    candidate_key, ann_dir / '_via_trials' / candidate_key)
                quality = _v32_via_trial_quality(payload[0], payload[1], candidate_key)
                trials.append(dict(trial_pattern_key=candidate_key, n_measurements=len(payload[0]),
                                   quality_key=repr(quality), error='',
                                   fallback_level=int(_v32_finite_float(payload[1].get('auto_fallback_level_used', 0), 0)),
                                   median_x_width_nm=_v32_safe_median([r.get('x_width_nm', np.nan) for r in payload[0]]),
                                   median_y_height_nm=_v32_safe_median([r.get('y_height_nm', np.nan) for r in payload[0]])))
                if best is None or quality > best[0]:
                    best = (quality, candidate_key, payload)
            except Exception as exc:
                trials.append(dict(trial_pattern_key=candidate_key, n_measurements=0,
                                   quality_key='', error=f'{type(exc).__name__}: {exc}'))
        if best is None:
            raise RuntimeError(f'Both V32 via trials failed: {trials}')
        _, key, payload = best
        rows, status, rejects = payload
        # Reuse the winning trial. A second run would spend the same image budget
        # again and could yield a different stochastic segmentation.
        status = dict(status)
        annotation = Path(status['annotated_path']) if status.get('annotated_path') else None
        if annotation is not None and annotation.is_file():
            selected_path = ann_dir / annotation.name
            _shutil_v32.copy2(annotation, selected_path)
            status['annotated_path'] = str(selected_path)
    else:
        if key not in PATTERN_KEYS_V19:
            raise ValueError(f'No valid pattern model for {image_path.name}: {key}')
        rows, status, rejects = _v32_process_v25_assigned(
            image_path, stage, job['condition'], job['region'], pattern, key, ann_dir)
        if pattern == 'via' and stage == 'before':
            trials.append(dict(trial_pattern_key=key, n_measurements=len(rows), quality_key='forced', error=''))
    for trial in trials:
        trial.update(condition=job['condition'], region=job['region'], before_image=image_path.name,
                     selected=trial['trial_pattern_key'] == key,
                     trial_nominal_nm=V32_VIA_NOMINAL_NM[trial['trial_pattern_key']],
                     mode='AUTO_V32' if not job.get('v25_pattern_key') else 'FORCED')
        for field in ('fallback_level', 'median_x_width_nm', 'median_y_height_nm'):
            trial.setdefault(field, math.nan)
    return rows, status, rejects, key, trials


def _v32_worker_main(request_path: str) -> int:
    request = _json_v32.loads(Path(request_path).read_text(encoding='utf-8'))
    try:
        result = {'ok': True, 'payload': _v32_measure_image_job(request)}
    except Exception as exc:
        result = {'ok': False, 'error': f'{type(exc).__name__}: {exc}',
                  'traceback': _traceback_v32.format_exc()}
    with Path(request['result_path']).open('wb') as stream:
        _pickle_v32.dump(result, stream, protocol=_pickle_v32.HIGHEST_PROTOCOL)
    return 0


def _v32_process_image_with_timeout(
    image_path: Path, stage: str, condition: str, region: int,
    public_pattern: str, v25_pattern_key: Optional[str], ann_dir: Path,
    timeout_seconds: float = 300.0,
):
    """Hard per-image wall clock limit, portable to Windows and native OpenCV."""
    if not math.isfinite(timeout_seconds) or not 0 < timeout_seconds <= 300:
        raise ValueError('image timeout must be > 0 and <= 300 seconds')
    started = time.monotonic()
    deadline = started + timeout_seconds
    with _tempfile_v32.TemporaryDirectory(prefix='sem_v32_image_') as temporary:
        root = Path(temporary)
        scratch_annotations = root / 'annotated'
        request = dict(image_path=str(Path(image_path).resolve()), stage=stage,
                       condition=str(condition), region=int(region), public_pattern=public_pattern,
                       v25_pattern_key=v25_pattern_key, ann_dir=str(scratch_annotations),
                       result_path=str(root / 'result.pkl'))
        request_path = root / 'request.json'
        request_path.write_text(_json_v32.dumps(request, ensure_ascii=False), encoding='utf-8')
        process = _subprocess_v32.Popen([sys.executable, '-u', str(Path(__file__).resolve()),
                                        '--_v32-worker', str(request_path)])
        last_progress = started
        try:
            while True:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise V32ImageTimeout(f'{Path(image_path).name}: exceeded {timeout_seconds:g}s; image skipped')
                try:
                    process.wait(timeout=min(1.0, remaining))
                    break
                except _subprocess_v32.TimeoutExpired:
                    now = time.monotonic()
                    if now - last_progress >= 15:
                        print(f'    V32 {Path(image_path).name}: elapsed={now-started:.0f}s / {timeout_seconds:g}s', flush=True)
                        last_progress = now
            if time.monotonic() > deadline:
                raise V32ImageTimeout(f'{Path(image_path).name}: exceeded {timeout_seconds:g}s; image skipped')
            result_path = Path(request['result_path'])
            if process.returncode != 0 or not result_path.is_file():
                raise RuntimeError(f'Image worker exited without a result (exit={process.returncode})')
            # This file is created only by our child in a private temporary directory.
            with result_path.open('rb') as stream:
                result = _pickle_v32.load(stream)
            if not result['ok']:
                raise RuntimeError(result['error'] + '\n' + result['traceback'])
            rows, status, rejects, key, trials = result['payload']
            # Publish annotations only for completed jobs. Timed-out work cannot
            # leave a partial annotated image that appears to be a valid result.
            for source in scratch_annotations.rglob('*'):
                if source.is_file():
                    destination = Path(ann_dir) / source.relative_to(scratch_annotations)
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    _shutil_v32.copy2(source, destination)
            if status.get('annotated_path'):
                status['annotated_path'] = str(Path(ann_dir) / Path(status['annotated_path']).relative_to(scratch_annotations))
            status['processing_elapsed_seconds'] = time.monotonic() - started
            status['image_timeout_seconds'] = timeout_seconds
            return rows, status, rejects, key, trials
        finally:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=2)
                except _subprocess_v32.TimeoutExpired:
                    process.kill()
                    process.wait(timeout=2)


def main_v32(
    before_root: Path,
    after_root: Path,
    output_root: Optional[Path] = None,
    via_pattern_mode: str = "via60",
    image_timeout_seconds: float = 300.0,
):
    if not math.isfinite(image_timeout_seconds) or not 0 < image_timeout_seconds <= 300:
        raise ValueError('image_timeout_seconds must be > 0 and <= 300')
    if via_pattern_mode not in ("auto", "via40", "via60"):
        raise ValueError("via_pattern_mode must be via60, via40 or auto")
    before_root = Path(before_root)
    after_root = Path(after_root)
    # AFTER determines the condition list, including folders missing from BEFORE.
    conditions = _v32_discover_conditions(after_root)
    if output_root is None:
        output_root = after_root.parent / "SEM_0830_0831_compare_V32"
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    # V25 uses these globals for labels; keep the discovered names intact.
    global CONDITION_IDS_V19, CONDITION_LABELS_V19
    CONDITION_IDS_V19 = [str(c) for c in conditions]
    CONDITION_LABELS_V19 = {str(c): f"Condition {c}" for c in conditions}

    ann_root = output_root / "annotated"
    plot_dir = output_root / "plots"
    trial_root = output_root / "_via_pattern_trials"
    ann_root.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 100)
    print("V32 | SEM comparison | V27 geometry + Yakun edge evidence + per-image timeout")
    print(f"Image time limit: {image_timeout_seconds:g}s (including all fallback / via trials)")
    print(f"Before : {before_root}")
    print(f"After  : {after_root}")
    print(f"Output : {output_root}")
    print(f"Conditions (from after): {', '.join(conditions)}")
    print(f"Via mode: {via_pattern_mode}")
    print("=" * 100)

    inventory_df, preflight_errors_df = _v32_build_inventory(before_root, after_root, conditions)
    _v32_write_csv(inventory_df, output_root / "inventory_mapping.csv")
    _v32_write_csv(preflight_errors_df, output_root / "processing_errors.csv")
    if not preflight_errors_df.empty:
        for error in preflight_errors_df.to_dict("records"):
            print(f"  PREFLIGHT [{error['error_type']}]: {error['message']}")
    if inventory_df.empty:
        raise RuntimeError(
            "没有建立任何有效 inventory。请检查对应子文件夹是否存在、是否包含 TIF/TIFF/PNG/JPG/JPEG 图片，"
            "以及 before/after 照片数量是否一致；详情见 processing_errors.csv。"
        )

    all_rows: List[dict] = []
    status_rows: List[dict] = []
    rejected_rows: List[dict] = []
    processing_errors: List[dict] = preflight_errors_df.to_dict("records") if not preflight_errors_df.empty else []
    via_choice_rows: List[dict] = []

    # Process pair-by-pair so the same auto-selected via40/via60 model is used for
    # BEFORE and AFTER of the corresponding condition + region.
    pair_keys = inventory_df[["condition", "region", "pattern", "pair_key"]].drop_duplicates().copy()
    _pattern_order_v32 = {p: i for i, p in enumerate(V32_PATTERN_SEQUENCE)}
    pair_keys["_pattern_order_v32"] = pair_keys["pattern"].map(_pattern_order_v32).fillna(999)
    pair_keys = _v32_sort_conditions(pair_keys, ["condition", "region", "_pattern_order_v32"]).drop(columns=["_pattern_order_v32"])

    for _, pk in pair_keys.iterrows():
        condition = str(pk["condition"])
        region = int(pk["region"])
        pattern = str(pk["pattern"])
        pair_key = str(pk["pair_key"])
        pair_inv = inventory_df[inventory_df["pair_key"] == pair_key]

        print(f"\n[{pair_key}] {V32_PATTERN_ZH[pattern]}")
        stage_paths = {}
        valid = True
        for stage in V32_STAGES:
            q = pair_inv[pair_inv["stage"] == stage]
            if q.empty:
                valid = False
                processing_errors.append({
                    "stage": stage,
                    "condition": condition,
                    "region": region,
                    "pattern": pattern,
                    "pair_key": pair_key,
                    "error_type": "MISSING_STAGE_INVENTORY",
                    "message": "该 pair 缺少 before 或 after inventory。",
                })
                continue
            row0 = q.iloc[0]
            if str(row0.get("pixel_size_preflight_note", "")):
                valid = False
                processing_errors.append({
                    "stage": stage,
                    "condition": condition,
                    "region": region,
                    "pattern": pattern,
                    "pair_key": pair_key,
                    "image": row0.get("image", ""),
                    "path": row0.get("path", ""),
                    "error_type": "SKIP_DUE_TO_PIXELSIZE_PREFLIGHT_ERROR",
                    "message": row0.get("pixel_size_preflight_note", ""),
                })
                continue
            stage_paths[stage] = Path(row0["path"])
        if not valid or len(stage_paths) != 2:
            print("  SKIP: before/after inventory 或 PixelSize 预检不完整")
            continue

        # Via auto-choice and its winning result share ONE before-image worker
        # and deadline. AFTER uses the same key only if BEFORE completed.
        if pattern in V32_FIXED_PATTERN_KEY:
            selected_key = V32_FIXED_PATTERN_KEY[pattern]
        else:
            selected_key = via_pattern_mode if via_pattern_mode in V32_VIA_KEYS else None

        for stage in V32_STAGES:
            image_path = stage_paths[stage]
            ann_dir = ann_root / stage / condition / f"Region_{region}"
            seq_index = int(pair_inv[pair_inv["stage"] == stage].iloc[0]["sequence_index"])
            if pattern == 'via' and stage == 'after' and selected_key is None:
                note = 'BEFORE did not complete via model selection; matching AFTER skipped.'
                status_rows.append(dict(stage=stage, condition=condition, region=region,
                                        pattern=pattern, pair_key=pair_key, image=image_path.name,
                                        path=str(image_path), status='SKIPPED_VIA_MODEL_UNAVAILABLE', note=note))
                processing_errors.append(dict(stage=stage, condition=condition, region=region,
                                               pattern=pattern, pair_key=pair_key, image=image_path.name,
                                               error_type='VIA_MODEL_UNAVAILABLE', message=note))
                print(f'  SKIP {image_path.name}: {note}', flush=True)
                continue
            print(f"  {V32_STAGE_ZH[stage]} {image_path.name} -> {selected_key}")
            try:
                rows, status, rejects, used_key, trials = _v32_process_image_with_timeout(
                    image_path=image_path,
                    stage=stage,
                    condition=condition,
                    region=region,
                    public_pattern=pattern,
                    v25_pattern_key=selected_key,
                    ann_dir=ann_dir,
                    timeout_seconds=image_timeout_seconds,
                )
                selected_key = used_key
                for trial in trials:
                    trial.update(pair_key=pair_key, after_image=stage_paths['after'].name)
                via_choice_rows.extend(trials)
                for r in rows:
                    r["sequence_index"] = seq_index
                    r["selected_pattern_key"] = selected_key
                    r["selected_v25_pattern_key"] = selected_key
                status["sequence_index"] = seq_index
                status["selected_pattern_key"] = selected_key
                status["selected_v25_pattern_key"] = selected_key
                all_rows.extend(rows)
                status_rows.append(status)
                rejected_rows.extend(rejects)
                print(f"    accepted={len(rows)} status={status.get('status', '')}")
            except Exception as exc:
                processing_errors.append({
                    "stage": stage,
                    "condition": condition,
                    "region": region,
                    "pattern": pattern,
                    "pair_key": pair_key,
                    "image": image_path.name,
                    "path": str(image_path),
                    "selected_pattern_key": selected_key,
                    "selected_v25_pattern_key": selected_key,
                    "error_type": 'IMAGE_TIMEOUT' if isinstance(exc, V32ImageTimeout) else type(exc).__name__,
                    "message": str(exc),
                    "traceback": _traceback_v32.format_exc(),
                })
                status_rows.append({
                    "stage": stage,
                    "stage_zh": V32_STAGE_ZH[stage],
                    "condition": condition,
                    "condition_label": f"Condition {condition}",
                    "region": region,
                    "pair_key": pair_key,
                    "pattern": pattern,
                    "pattern_zh": V32_PATTERN_ZH[pattern],
                    "image": image_path.name,
                    "path": str(image_path),
                    "selected_pattern_key": selected_key,
                    "selected_v25_pattern_key": selected_key,
                    "source_algorithm": "V32_V27_Yakun_fusion",
                    "status": "SKIPPED_TIMEOUT" if isinstance(exc, V32ImageTimeout) else "ERROR",
                    "image_timeout_seconds": image_timeout_seconds,
                    "note": f"{type(exc).__name__}: {exc}",
                })
                label = 'TIMEOUT -> SKIP' if isinstance(exc, V32ImageTimeout) else 'ERROR'
                print(f"    {label} {type(exc).__name__}: {exc}", flush=True)
            # Keep completed statuses and skip reasons available during long runs.
            _v32_write_csv(_v32_df(status_rows), output_root / 'image_status.csv')
            _v32_write_csv(_v32_df(processing_errors), output_root / 'processing_errors.csv')

    objects_df = _v32_df(all_rows)
    status_df = _v32_df(status_rows)
    rejected_df = _v32_df(rejected_rows)
    errors_df = _v32_df(processing_errors)
    via_choice_df = _v32_df(via_choice_rows)

    if not objects_df.empty:
        sort_cols = [c for c in ["condition", "region", "pattern", "stage", "image", "object_id"] if c in objects_df.columns]
        objects_df = _v32_sort_conditions(objects_df, sort_cols)

    image_metrics_df = _v32_build_image_metrics(objects_df)
    image_pair_df = _v32_build_image_pair_comparison(image_metrics_df)
    overall_df = _v32_build_overall_comparison(objects_df)
    object_pair_df = _v32_build_object_comparison(objects_df)

    outputs = {
        "excel": output_root / "SEM_0830_0831_before_after_V32_results.xlsx",
        "inventory": output_root / "inventory_mapping.csv",
        "image_status": output_root / "image_status.csv",
        "all_objects": output_root / "all_object_measurements.csv",
        "image_metrics": output_root / "image_metric_summary.csv",
        "image_pair_comparison": output_root / "image_pair_comparison.csv",
        "overall_comparison": output_root / "overall_comparison.csv",
        "object_comparison": output_root / "object_before_after_comparison.csv",
        "rejected_objects": output_root / "rejected_objects.csv",
        "via_model_choice": output_root / "via_model_choice.csv",
        "errors": output_root / "processing_errors.csv",
        "settings": output_root / "settings.json",
    }

    _v32_write_csv(inventory_df, outputs["inventory"])
    _v32_write_csv(status_df, outputs["image_status"])
    _v32_write_csv(objects_df, outputs["all_objects"])
    _v32_write_csv(image_metrics_df, outputs["image_metrics"])
    _v32_write_csv(image_pair_df, outputs["image_pair_comparison"])
    _v32_write_csv(overall_df, outputs["overall_comparison"])
    _v32_write_csv(object_pair_df, outputs["object_comparison"])
    _v32_write_csv(rejected_df, outputs["rejected_objects"])
    _v32_write_csv(via_choice_df, outputs["via_model_choice"])
    _v32_write_csv(errors_df, outputs["errors"])

    settings = {
        "script_version": V32_SCRIPT_VERSION,
        "algorithm_base": "V27 geometry/output + Yakun gray edge primitives + V30 bounded recovery/timeout",
        "v32_edge_rule": "original-image local gray 50%; contour normals with circular tracking; no forced ellipse/rectangle",
        "v32_edge_min_support_fraction": 0.85,
        "v32_edge_max_missing_fraction": 0.08,
        "v32_illumination_rule": "robust exterior plane subtraction for contour profiles; enhanced images are localization only",
        "v32_slot_recovery_rule": "V27 first; V30 only missing columns or contours failing original-gray evidence; four supported rows required",
        "image_timeout_seconds": image_timeout_seconds,
        "timeout_scope": "one image, including fallback and both BEFORE via model trials; hard subprocess deadline",
        "v1_6_used": False,
        "before_root": str(before_root),
        "after_root": str(after_root),
        "output_root": str(output_root),
        "conditions": list(conditions),
        "condition_discovery_rule": "digit-only immediate subfolders of after_root, sorted numerically; match exact names in before_root",
        "condition_name_rule": "preserve leading zeros; ignore before-only folders and nonnumeric folders",
        "regions": sorted(int(r) for r in inventory_df["region"].unique()),
        "pattern_sequence_per_region": list(V32_PATTERN_SEQUENCE),
        "supported_image_extensions": sorted(V32_IMAGE_SUFFIXES),
        "required_image_count_per_condition_stage": None,
        "image_count_rule": "equal supported-image counts in corresponding before/after condition folders; formats may differ; no fixed count",
        "pairing_rule": "condition + natural filename order; region/pattern extend in groups of three from Region 2; filenames need not match",
        "pixel_size_rule": "same-stem TXT, PixelSize=<number>, parsed by V25 read_pixel_size_nm",
        "trench_v25_key": "trench160",
        "slot_v25_key": "slot210",
        "via_pattern_mode": via_pattern_mode,
        "cd_filter_policy": {
            "via": "nominal pattern retained; search scale from image candidates; no fixed 40/60-nm CD cutoff",
            "via_span_factors_of_search_scale": [V23_VIA_MIN_SPAN_FACTOR, V23_VIA_MAX_SPAN_FACTOR],
            "via_scale_requires_candidates": 4,
            "slot_relative_area_factors": [SLOT_AREA_MIN_MEAN_FACTOR, SLOT_AREA_MAX_MEAN_FACTOR],
            "slot_lowgray_max_image_area_fraction": V18_SLOT_LOWGRAY_MODEL_AREA_MAX_FRAC,
            "trench_width_fraction_tolerance": V20_TRENCH_WIDTH_FRAC_TOL,
            "trench_gap_fraction_tolerance": V20_TRENCH_GAP_FRAC_TOL,
            "trench_width_gap_mad_factors": [V20_TRENCH_WIDTH_MAD_FACTOR, V20_TRENCH_GAP_MAD_FACTOR],
            "slot_trench_fixed_nm_limits": None,
            "preserved_checks": "original-image edges, shape, complete four-row slots, lattice and image/neighbor boundaries",
        },
        "via_auto_rule": "test via40 and via60 using V32 on BEFORE only; select by validated result quality; freeze same key for AFTER",
        "metric_specs": V32_METRIC_SPECS,
        "accepted_object_count": int(len(objects_df)),
        "image_metric_rows": int(len(image_metrics_df)),
        "image_pair_rows": int(len(image_pair_df)),
        "object_pair_rows": int(len(object_pair_df)),
        "error_or_warning_rows": int(len(errors_df)),
    }
    outputs["settings"].write_text(_json_v32.dumps(settings, ensure_ascii=False, indent=2), encoding="utf-8")

    with pd.ExcelWriter(outputs["excel"], engine="openpyxl") as writer:
        inventory_df.to_excel(writer, sheet_name="inventory", index=False)
        status_df.to_excel(writer, sheet_name="image_status", index=False)
        objects_df.to_excel(writer, sheet_name="all_objects", index=False)
        image_metrics_df.to_excel(writer, sheet_name="image_metrics", index=False)
        image_pair_df.to_excel(writer, sheet_name="image_pair_compare", index=False)
        overall_df.to_excel(writer, sheet_name="overall_compare", index=False)
        object_pair_df.to_excel(writer, sheet_name="object_compare", index=False)
        rejected_df.to_excel(writer, sheet_name="rejected_objects", index=False)
        via_choice_df.to_excel(writer, sheet_name="via_model_choice", index=False)
        errors_df.to_excel(writer, sheet_name="errors", index=False)
        pd.DataFrame([settings]).to_excel(writer, sheet_name="settings", index=False)
        _autofit_excel_v19(writer)

    _v32_save_pair_plots(image_pair_df, plot_dir)

    print("\n" + "=" * 100)
    print("V32 处理完成")
    print(f"Accepted objects        : {len(objects_df)}")
    print(f"Image metric rows       : {len(image_metrics_df)}")
    print(f"Image pair compare rows : {len(image_pair_df)}")
    print(f"Object compare rows     : {len(object_pair_df)}")
    print(f"Errors/warnings         : {len(errors_df)}")
    print(f"Excel                   : {outputs['excel']}")
    print(f"Annotated               : {ann_root}")
    print(f"Plots                   : {plot_dir}")
    print("=" * 100)
    return outputs


# -----------------------------------------------------------------------------
# V32 CLI
# -----------------------------------------------------------------------------
def _v32_build_parser():
    parser = _argparse_v32.ArgumentParser(
        description=(
            "V32: SEM comparison using numeric-named folders discovered under AFTER, "
            "matched by exact folder name in BEFORE, with equal image counts (TIF/TIFF/PNG/JPG/JPEG). "
            "V32 combines V27 geometry, Yakun edges and V30 bounded recovery; 300s image limit."
        )
    )
    parser.add_argument(
        "--before", type=Path, default=V32_DEFAULT_BEFORE_ROOT,
        help=f"处理前根目录；默认：{V32_DEFAULT_BEFORE_ROOT}",
    )
    parser.add_argument(
        "--after", type=Path, default=V32_DEFAULT_AFTER_ROOT,
        help=f"处理后根目录；自动识别其纯数字名称子文件夹；默认：{V32_DEFAULT_AFTER_ROOT}",
    )
    parser.add_argument(
        "--output", type=Path, default=None,
        help="输出目录；默认在两个输入目录同级建立 SEM_0830_0831_compare_V32",
    )
    parser.add_argument(
        "--image-timeout-seconds", type=float, default=300.0,
        help="单张图处理上限，含自动补检及 via 模型尝试；默认 300 秒，可设为更短的正数。",
    )
    parser.add_argument(
        "--via-pattern", choices=("auto", "via40", "via60"), default="via60",
        help=(
            "Via 模型，默认 via60；60 是图案标称值，实际 CD 按原图测量。auto=分别试 via40/via60，"
            "选更可靠者并固定给对应 AFTER；只有明确需要时才使用 auto 或 via40。"
        ),
    )
    return parser


def _v32_cli(argv=None) -> int:
    args = _v32_build_parser().parse_args(argv)
    try:
        before_root = Path(args.before)
        after_root = Path(args.after)
        output_root = Path(args.output) if args.output is not None else after_root.parent / "SEM_0830_0831_compare_V32"
        main_v32(
            before_root=before_root,
            after_root=after_root,
            output_root=output_root,
            via_pattern_mode=str(args.via_pattern),
            image_timeout_seconds=float(args.image_timeout_seconds),
        )
        return 0
    except KeyboardInterrupt:
        print("用户中断。", file=sys.stderr)
        return 130
    except Exception as exc:
        print(f"V32 程序失败：{type(exc).__name__}: {exc}", file=sys.stderr)
        _traceback_v32.print_exc()
        return 1


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == '--_v32-worker':
        raise SystemExit(_v32_worker_main(sys.argv[2]))
    raise SystemExit(_v32_cli())
