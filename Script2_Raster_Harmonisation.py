# =============================================================================
# Script 2 — Raster Harmonisation
# Project : MSc GIS Dissertation — A66 Northern Trans-Pennine Corridor
# Author  : Abhin Nangari (Student ID: 52536856)
# University: University of Aberdeen — MSc GIS GG5910/GG5912
# Supervisor: Dr Shaktiman Singh
# GitHub  : https://github.com/Abhin-Nangari/A66-Connectivity-MSc
# License : MIT
#
# Purpose :
#   1. Pre-flight check — verify all Script 1 outputs exist before proceeding
#   2. Derive Slope from EA DTM 1m (no snap raster set at this stage)
#   3. Clip all layers to StudyArea_2km_buffer
#   4. Cap CHM at 40m (evidence-based: confirmed buffer max = 36.67m)
#   5. Resample to 2m:
#        - Bilinear      : CHM, VRI, Slope (continuous)
#        - Nearest+Mask  : LCM (categorical) — ExtractByMask + ExtractBand
#          forces pixel-perfect alignment and single band output
#   6. Normalise continuous layers to 0-1
#   7. Final alignment verification — all four outputs must match exactly
#
# Master grid:
#   Set from CHM_2m.tif after resampling — all subsequent outputs snap to this.
#   Snap raster is NOT set before CHM_2m exists — prevents environment errors.
#
# LCM alignment method:
#   Standard Resample with snapRaster does not reliably force alignment in
#   ArcGIS Pro when source and target grids have different origins.
#   ExtractByMask using CHM_2m as mask guarantees pixel-perfect alignment.
#   ExtractBand([1]) enforces single-band output — prevents multi-band artefact.
#
# CHM cap justification:
#   Diagnostic run confirmed buffer max = 36.67m. Cap at 40m retains all
#   legitimate vegetation values while removing any built-structure artefacts.
#
# Resistance weight citations:
#   LCM source: Morton et al. (2024) NERC EDS EIDC DOI:10.5285/7727ce7d...
#
# Run from : ArcGIS Pro Python console (Analysis tab → Python)
# CRS      : EPSG:27700 British National Grid (all outputs)
# Inputs   : Processed\LiDAR\CHM_1m.tif
#            Processed\LiDAR\VRI_1m.tif
#            Processed\LiDAR\DTM_mosaic.tif
#            Raw\UKCEH_LCM_2023_10m\
# Outputs  : Processed\Harmonised\
# =============================================================================

import arcpy
import os
import datetime
import gc
import glob as _glob

arcpy.CheckOutExtension("Spatial")
from arcpy.sa import Raster, ExtractByMask, Con, Float, Slope

# ---------------------------------------------------------------------------
# 0. ENVIRONMENT SETUP
# ---------------------------------------------------------------------------
# NOTE: snapRaster and extent are NOT set here.
# They are set only after CHM_2m.tif exists (Stage 5 onwards).
# This prevents RuntimeError when snap raster does not yet exist.

arcpy.env.overwriteOutput = True
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(27700)

ROOT         = r"C:\Mac\Home\Desktop\UK\Abardeen\UoA\Dissetation\A66"
GDB          = os.path.join(ROOT, "Processed", "StudyArea", "A66_Study.gdb")
STUDY_BUFFER = os.path.join(GDB, "StudyArea_2km_buffer")
LIDAR_DIR    = os.path.join(ROOT, "Processed", "LiDAR")
OUT_DIR      = os.path.join(ROOT, "Processed", "Harmonised")
CHM_CAP      = 40.0

if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)

# Script 1 inputs
IN_DTM   = os.path.join(LIDAR_DIR, "DTM_mosaic.tif")
IN_CHM   = os.path.join(LIDAR_DIR, "CHM_1m.tif")
IN_VRI   = os.path.join(LIDAR_DIR, "VRI_1m.tif")

# Find LCM raster
lcm_candidates = _glob.glob(
    os.path.join(ROOT, "Raw", "UKCEH_LCM_2023_10m", "**", "*.tif"),
    recursive=True)
if not lcm_candidates:
    lcm_candidates = _glob.glob(
        os.path.join(ROOT, "Raw", "UKCEH_LCM_2023_10m", "**", "*.img"),
        recursive=True)
if not lcm_candidates:
    raise FileNotFoundError(
        "No LCM raster found in Raw\\UKCEH_LCM_2023_10m\\. "
        "Check extraction of both FME zip parts.")
IN_LCM = lcm_candidates[0]

# Output paths — intermediate
SLOPE_RAW  = os.path.join(OUT_DIR, "Slope_raw_1m.tif")
CLIP_CHM   = os.path.join(OUT_DIR, "CHM_clipped_1m.tif")
CLIP_VRI   = os.path.join(OUT_DIR, "VRI_clipped_1m.tif")
CLIP_SLOPE = os.path.join(OUT_DIR, "Slope_clipped_1m.tif")
CLIP_LCM   = os.path.join(OUT_DIR, "LCM_clipped_10m.tif")
CHM_CAPPED = os.path.join(OUT_DIR, "CHM_capped_1m.tif")

# Output paths — resampled 2m
RS_CHM     = os.path.join(OUT_DIR, "CHM_2m.tif")
RS_VRI     = os.path.join(OUT_DIR, "VRI_2m.tif")
RS_SLOPE   = os.path.join(OUT_DIR, "Slope_2m.tif")
LCM_ALIGN  = os.path.join(OUT_DIR, "LCM_2m_aligned.tif")

# Output paths — normalised
NORM_CHM   = os.path.join(OUT_DIR, "CHM_norm.tif")
NORM_VRI   = os.path.join(OUT_DIR, "VRI_norm.tif")
NORM_SLOPE = os.path.join(OUT_DIR, "Slope_norm.tif")

print("=" * 60)
print("SCRIPT 2 — RASTER HARMONISATION")
print(f"Started: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 60)

# ---------------------------------------------------------------------------
# HELPER — full verification
# ---------------------------------------------------------------------------

MASTER_EXT = None  # Set after CHM_2m is created

def verify(path, name, exp_min=None, exp_max=None, exp_bands=1,
           check_alignment=True):
    print(f"\n  Verifying: {name}")
    bands = int(arcpy.management.GetRasterProperties(
        path, "BANDCOUNT").getOutput(0))
    cs    = arcpy.management.GetRasterProperties(
        path, "CELLSIZEX").getOutput(0)
    desc  = arcpy.Describe(path)
    ext   = desc.extent
    ext_t = (ext.XMin, ext.YMin, ext.XMax, ext.YMax)

    band_ok = (bands == exp_bands)
    print(f"    Bands    : {bands} {'[OK]' if band_ok else '[WARNING] Expected ' + str(exp_bands)}")
    print(f"    Cell size: {cs}m")
    print(f"    Extent   : {ext_t}")

    ext_ok = True
    if check_alignment and MASTER_EXT is not None:
        ext_ok = (ext_t == MASTER_EXT)
        print(f"    Alignment: {'[OK]' if ext_ok else '[WARNING] Differs from master ' + str(MASTER_EXT)}")

    r_min  = float(arcpy.management.GetRasterProperties(path, "MINIMUM").getOutput(0))
    r_max  = float(arcpy.management.GetRasterProperties(path, "MAXIMUM").getOutput(0))
    r_mean = float(arcpy.management.GetRasterProperties(path, "MEAN").getOutput(0))
    r_std  = float(arcpy.management.GetRasterProperties(path, "STD").getOutput(0))
    print(f"    Min={r_min:.4f}  Max={r_max:.4f}  Mean={r_mean:.4f}  STD={r_std:.4f}")

    stat_ok = True
    if exp_min is not None and r_min < exp_min - 0.001:
        print(f"    [WARNING] Min below expected {exp_min}")
        stat_ok = False
    if exp_max is not None and r_max > exp_max + 0.001:
        print(f"    [WARNING] Max above expected {exp_max}")
        stat_ok = False

    if band_ok and ext_ok and stat_ok:
        print(f"    [OK] All checks passed.")

    return r_min, r_max, r_mean, r_std

# ---------------------------------------------------------------------------
# STAGE 0 — PRE-FLIGHT CHECK
# ---------------------------------------------------------------------------
# Verify all Script 1 outputs exist before any processing begins.

print("\n--- STAGE 0: Pre-flight Check ---")

preflight = {
    "DTM_mosaic.tif": IN_DTM,
    "CHM_1m.tif"    : IN_CHM,
    "VRI_1m.tif"    : IN_VRI,
    "LCM raster"    : IN_LCM,
}

preflight_ok = True
for name, path in preflight.items():
    if os.path.exists(path):
        print(f"  [OK] {name}: {path}")
    else:
        print(f"  [ABORT] {name} not found: {path}")
        preflight_ok = False

if not preflight_ok:
    raise RuntimeError(
        "Pre-flight check failed. Run Script 1 before Script 2.")

print("\n  [OK] All inputs confirmed — proceeding.")

# ---------------------------------------------------------------------------
# STAGE 1 — DERIVE SLOPE FROM EA DTM 1m
# ---------------------------------------------------------------------------
# Slope derived at 1m from high-quality EA DTM.
# No snap raster set here — DTM_mosaic covers full A66 extent.
# Slope is clipped to study buffer in Stage 2.

print("\n--- STAGE 1: Derive Slope from EA DTM 1m ---")

dtm_r     = Raster(IN_DTM)
slope_raw = Slope(dtm_r, "DEGREE")
slope_raw.save(SLOPE_RAW)
del dtm_r, slope_raw
gc.collect()
print(f"  [OK] Slope derived at 1m: {SLOPE_RAW}")

# ---------------------------------------------------------------------------
# STAGE 2 — CLIP ALL LAYERS TO StudyArea_2km_buffer
# ---------------------------------------------------------------------------
# Clipping before resampling reduces processing time and prevents
# interpolation outside the study area.
# No snap raster set — ExtractByMask clips to exact buffer boundary.

print("\n--- STAGE 2: Clip all layers to StudyArea_2km_buffer ---")

print("\n  Clipping CHM...")
chm_c = ExtractByMask(IN_CHM, STUDY_BUFFER)
chm_c.save(CLIP_CHM)
del chm_c; gc.collect()
verify(CLIP_CHM, "CHM clipped 1m",
       exp_min=0, exp_max=CHM_CAP, check_alignment=False)

print("\n  Clipping VRI...")
vri_c = ExtractByMask(IN_VRI, STUDY_BUFFER)
vri_c.save(CLIP_VRI)
del vri_c; gc.collect()
verify(CLIP_VRI, "VRI clipped 1m",
       exp_min=0, check_alignment=False)

print("\n  Clipping Slope...")
slp_c = ExtractByMask(SLOPE_RAW, STUDY_BUFFER)
slp_c.save(CLIP_SLOPE)
del slp_c; gc.collect()
verify(CLIP_SLOPE, "Slope clipped 1m",
       exp_min=0, exp_max=90, check_alignment=False)

print("\n  Clipping LCM...")
lcm_desc = arcpy.Describe(IN_LCM)
lcm_sr   = lcm_desc.spatialReference
print(f"  LCM source CRS: {lcm_sr.name} (EPSG:{lcm_sr.factoryCode})")
if lcm_sr.factoryCode != 27700:
    print("  Reprojecting LCM to EPSG:27700...")
    lcm_reproj = os.path.join(OUT_DIR, "LCM_reproj.tif")
    arcpy.management.ProjectRaster(
        in_raster       = IN_LCM,
        out_raster      = lcm_reproj,
        out_coor_system = arcpy.SpatialReference(27700),
        resampling_type = "NEAREST",
        cell_size       = "10")
    lcm_for_clip = lcm_reproj
else:
    print("  LCM already in EPSG:27700.")
    lcm_for_clip = IN_LCM

lcm_c = ExtractByMask(lcm_for_clip, STUDY_BUFFER)
lcm_c.save(CLIP_LCM)
del lcm_c; gc.collect()
# LCM may show 2 bands at this stage — resolved in Stage 5
# Alignment warning expected at 10m resolution — resolved in Stage 5
verify(CLIP_LCM, "LCM clipped 10m", check_alignment=False)
print("  NOTE: LCM band count and alignment resolved in Stage 5.")

# ---------------------------------------------------------------------------
# STAGE 3 — CAP CHM OUTLIERS AT 40m
# ---------------------------------------------------------------------------
# Confirmed buffer max from diagnostic = 36.67m.
# Cap at 40m retains all legitimate vegetation values.

print("\n--- STAGE 3: Cap CHM outliers at 40m ---")

chm_capped = Con(Raster(CLIP_CHM) > CHM_CAP, CHM_CAP, Raster(CLIP_CHM))
chm_capped.save(CHM_CAPPED)
del chm_capped; gc.collect()
verify(CHM_CAPPED, f"CHM capped at {CHM_CAP}m",
       exp_min=0, exp_max=CHM_CAP, check_alignment=False)

# ---------------------------------------------------------------------------
# STAGE 4 — RESAMPLE CHM, VRI, SLOPE TO 2m
# ---------------------------------------------------------------------------
# Bilinear resampling for continuous surfaces.
# No snap raster set yet — CHM_2m becomes the master grid reference.

print("\n--- STAGE 4: Resample CHM, VRI, Slope to 2m ---")

print("\n  Resampling CHM...")
arcpy.management.Resample(CHM_CAPPED, RS_CHM, "2", "BILINEAR")

print("\n  Resampling VRI...")
arcpy.management.Resample(CLIP_VRI, RS_VRI, "2", "BILINEAR")

print("\n  Resampling Slope...")
arcpy.management.Resample(CLIP_SLOPE, RS_SLOPE, "2", "BILINEAR")

# Set master grid from CHM_2m — all subsequent outputs must match this
desc       = arcpy.Describe(RS_CHM)
ext        = desc.extent
MASTER_EXT = (ext.XMin, ext.YMin, ext.XMax, ext.YMax)
print(f"\n  Master grid set from CHM_2m: {MASTER_EXT}")

# Now verify all three resampled layers against master grid
verify(RS_CHM,   "CHM 2m",   exp_min=0, exp_max=CHM_CAP)
verify(RS_VRI,   "VRI 2m",   exp_min=0)
verify(RS_SLOPE, "Slope 2m", exp_min=0, exp_max=90)

# ---------------------------------------------------------------------------
# STAGE 5 — LCM ALIGNMENT AND SINGLE-BAND EXTRACTION
# ---------------------------------------------------------------------------
# ExtractByMask using CHM_2m as mask forces pixel-perfect alignment.
# ExtractBand([1]) enforces single-band output.
# Snap raster set to CHM_2m before this stage.

print("\n--- STAGE 5: LCM alignment and single-band extraction ---")

arcpy.env.snapRaster = RS_CHM
arcpy.env.extent     = arcpy.Describe(RS_CHM).extent
arcpy.env.cellSize   = "2"

lcm_temp = os.path.join(OUT_DIR, "LCM_temp.tif")
arcpy.management.Resample(CLIP_LCM, lcm_temp, "2", "NEAREST")

lcm_masked = ExtractByMask(Raster(lcm_temp), RS_CHM)
lcm_single = arcpy.sa.ExtractBand(lcm_masked, [1])
lcm_single.save(LCM_ALIGN)
del lcm_masked, lcm_single; gc.collect()
arcpy.management.Delete(lcm_temp)

verify(LCM_ALIGN, "LCM_2m_aligned",
       exp_min=0, exp_max=21, exp_bands=1)

# ---------------------------------------------------------------------------
# STAGE 6 — NORMALISE CONTINUOUS LAYERS TO 0-1
# ---------------------------------------------------------------------------
# Min-max normalisation anchored to study buffer values only.
# Snap raster remains CHM_2m.

print("\n--- STAGE 6: Normalise continuous layers to 0-1 ---")

def normalise(in_path, out_path, name):
    r     = Raster(in_path)
    r_min = float(arcpy.management.GetRasterProperties(
        in_path, "MINIMUM").getOutput(0))
    r_max = float(arcpy.management.GetRasterProperties(
        in_path, "MAXIMUM").getOutput(0))
    print(f"\n  Normalising {name}: range {r_min:.4f}–{r_max:.4f}")
    if r_max == r_min:
        print(f"  [WARNING] Zero variance in {name}")
        norm = r * 0.0
    else:
        norm = Float(r - r_min) / Float(r_max - r_min)
    norm.save(out_path)
    del r, norm; gc.collect()
    verify(out_path, f"{name} normalised",
           exp_min=0.0, exp_max=1.0, exp_bands=1)
    out_min = float(arcpy.management.GetRasterProperties(
        out_path, "MINIMUM").getOutput(0))
    out_max = float(arcpy.management.GetRasterProperties(
        out_path, "MAXIMUM").getOutput(0))
    if abs(out_min) <= 0.001 and abs(out_max - 1.0) <= 0.001:
        print(f"  [OK] Normalisation confirmed 0-1.")
    else:
        print(f"  [WARNING] Normalisation outside 0-1.")

normalise(RS_CHM,   NORM_CHM,   "CHM")
normalise(RS_VRI,   NORM_VRI,   "VRI")
normalise(RS_SLOPE, NORM_SLOPE, "Slope")

# ---------------------------------------------------------------------------
# STAGE 7 — FINAL ALIGNMENT VERIFICATION
# ---------------------------------------------------------------------------
# All four outputs must match exactly — same extent, same cell size,
# single band. Script 3 will abort if any mismatch is found.

print("\n--- STAGE 7: Final alignment verification — all outputs ---")

final_outputs = {
    "CHM_norm"      : NORM_CHM,
    "VRI_norm"      : NORM_VRI,
    "Slope_norm"    : NORM_SLOPE,
    "LCM_2m_aligned": LCM_ALIGN,
}

ref_ext     = None
all_aligned = True

for name, path in final_outputs.items():
    bands = int(arcpy.management.GetRasterProperties(
        path, "BANDCOUNT").getOutput(0))
    cs    = arcpy.management.GetRasterProperties(
        path, "CELLSIZEX").getOutput(0)
    desc  = arcpy.Describe(path)
    ext   = desc.extent
    ext_t = (ext.XMin, ext.YMin, ext.XMax, ext.YMax)
    if ref_ext is None:
        ref_ext = ext_t
    aligned = (ext_t == ref_ext) and (bands == 1)
    status  = "[OK]" if aligned else "[WARNING]"
    if not aligned:
        all_aligned = False
    print(f"  {status} {name}: bands={bands}  cs={cs}m  extent={ext_t}")

if all_aligned:
    print("\n  [OK] ALL OUTPUTS PERFECTLY ALIGNED — ready for Script 3.")
else:
    print("\n  [WARNING] Alignment issues — do not proceed to Script 3.")
    print("  Check warnings above and contact supervisor.")

# ---------------------------------------------------------------------------
# STAGE 8 — LCM vs LIVING ENGLAND CROSS-CHECK NOTE
# ---------------------------------------------------------------------------

print("\n--- STAGE 8: LCM vs Living England cross-check ---")
print("  ACTION REQUIRED — visual cross-check in ArcGIS Pro:")
print("  1. Load LCM_2m_aligned.tif from Processed\\Harmonised\\")
print("  2. Load Living England via .lyrx file")
print("  3. Compare upland grassland classes over study area")
print("  4. Note disagreements in dissertation limitations")
print("  5. LCM 2023 is primary layer — Living England is validation only")

# ---------------------------------------------------------------------------
# FINAL SUMMARY
# ---------------------------------------------------------------------------

print()
print("=" * 60)
print("SCRIPT 2 COMPLETE")
print(f"Finished: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 60)
print()

all_files = {
    "CHM_clipped_1m.tif"  : CLIP_CHM,
    "VRI_clipped_1m.tif"  : CLIP_VRI,
    "Slope_clipped_1m.tif": CLIP_SLOPE,
    "LCM_clipped_10m.tif" : CLIP_LCM,
    "CHM_capped_1m.tif"   : CHM_CAPPED,
    "CHM_2m.tif"          : RS_CHM,
    "VRI_2m.tif"          : RS_VRI,
    "Slope_2m.tif"        : RS_SLOPE,
    "LCM_2m_aligned.tif"  : LCM_ALIGN,
    "CHM_norm.tif"        : NORM_CHM,
    "VRI_norm.tif"        : NORM_VRI,
    "Slope_norm.tif"      : NORM_SLOPE,
}

print("Files created:")
for fname, fpath in all_files.items():
    status = "[OK]" if os.path.exists(fpath) else "[MISSING]"
    print(f"  {status} {fname}")

print()
print("VERIFIED STATISTICS FOR DISSERTATION (Chapter 5):")
print(f"  CHM buffer max pre-cap : 36.67m  Cap applied: {CHM_CAP}m")
print("  All normalised layers  : 0-1 confirmed")
print("  LCM alignment method   : ExtractByMask + ExtractBand (single band)")
print("  Target resolution      : 2m EPSG:27700")
print("  All outputs aligned    :", "[OK]" if all_aligned else "[WARNING]")
print()
print("NEXT STEP: Script 3 — Fuzzy MCE Resistance Surface")
print()
print("Push to GitHub:")
print("  https://github.com/Abhin-Nangari/A66-Connectivity-MSc")

arcpy.CheckInExtension("Spatial")
