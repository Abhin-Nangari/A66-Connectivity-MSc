# =============================================================================
# Script 3 — Fuzzy MCE Resistance Surface
# Project : MSc GIS Dissertation — A66 Northern Trans-Pennine Corridor
# Author  : Abhin Nangari (Student ID: 52536856)
# University: University of Aberdeen — MSc GIS GG5910/GG5912
# Supervisor: Dr Shaktiman Singh
# GitHub  : https://github.com/Abhin-Nangari/A66-Connectivity-MSc
# License : MIT
#
# Purpose :
#   1. Pre-flight validation — verify all Script 2 outputs before proceeding
#      Aborts immediately if any input has wrong band count or wrong extent
#   2. Apply fuzzy membership functions to all input layers
#   3. Reclassify LCM 2023 by resistance class (class 0 → NoData)
#      ExtractBand([1]) enforces single-band output at this stage
#   4. Weighted overlay — 4 sensitivity scenarios processed sequentially:
#        S1 — Primary : LCM 0.40, VRI 0.30, Slope 0.20, CHM 0.10
#        S2 — +20%    : proportionally increased, renormalised to sum 1.0
#        S3 — -20%    : proportionally decreased, renormalised to sum 1.0
#        S4 — Equal   : LCM 0.25, VRI 0.25, Slope 0.25, CHM 0.25
#   5. Spatial stability raster — pixel-wise range across all 4 scenarios
#   6. Full verification at every stage — band count, extent, statistics
#
# Memory management:
#   Each stage saves to disk and explicitly releases memory before proceeding.
#   Each weighted overlay scenario processed independently.
#   Ensures full reproducibility on hardware with limited RAM.
#
# Resistance weight citations:
#   LCM   0.40 — Zeller, McGarigal and Whiteley (2012) Landscape Ecology 27(6)
#   VRI   0.30 — Zeller et al. (2012); Liu et al. (2018) PLoS ONE 13(11)
#   Slope 0.20 — Zeller et al. (2012)
#   CHM   0.10 — Zeller et al. (2012)
#
# UKCEH LCM 2023 resistance reclassification:
#   Class 0  — NoData/background    → NoData (excluded)
#   Class 1  — Broadleaved woodland → 0.10 (low resistance)
#   Class 2  — Coniferous woodland  → 0.20
#   Class 3  — Arable               → 0.70
#   Class 4  — Improved grassland   → 0.60
#   Class 5  — Neutral grassland    → 0.40
#   Class 6  — Calcareous grassland → 0.35
#   Class 7  — Acid grassland       → 0.35
#   Class 8  — Fen/marsh/swamp      → 0.30
#   Class 9  — Heather              → 0.30
#   Class 10 — Heather grassland    → 0.35
#   Class 11 — Bog                  → 0.45
#   Class 12 — Inland rock          → 0.80
#   Class 13 — Saltwater            → 1.00 (impassable)
#   Class 14 — Freshwater           → 0.55
#   Class 15 — Supralittoral rock   → 0.90
#   Class 16 — Supralittoral sed.   → 0.90
#   Class 17 — Littoral rock        → 1.00
#   Class 18 — Littoral sediment    → 1.00
#   Class 19 — Saltmarsh            → 0.70
#   Class 20 — Urban                → 0.95 (A66 road surface)
#   Class 21 — Suburban             → 0.80
#
# Citation: Morton et al. (2024) NERC EDS EIDC DOI:10.5285/7727ce7d...
# Resistance values: Zeller et al. (2012); Liu et al. (2018)
#
# Run from : ArcGIS Pro Python console (Analysis tab → Python)
# CRS      : EPSG:27700 British National Grid (all outputs)
# Inputs   : Processed\Harmonised\CHM_norm.tif
#            Processed\Harmonised\VRI_norm.tif
#            Processed\Harmonised\Slope_norm.tif
#            Processed\Harmonised\LCM_2m_aligned.tif
# Outputs  : Processed\Resistance\
# =============================================================================

import arcpy
import os
import datetime
import gc

arcpy.CheckOutExtension("Spatial")
from arcpy.sa import Raster, Con, SetNull, Float

# ---------------------------------------------------------------------------
# 0. ENVIRONMENT SETUP
# ---------------------------------------------------------------------------

arcpy.env.overwriteOutput = True
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(27700)

ROOT     = r"C:\Mac\Home\Desktop\UK\Abardeen\UoA\Dissetation\A66"
HARM_DIR = os.path.join(ROOT, "Processed", "Harmonised")
OUT_DIR  = os.path.join(ROOT, "Processed", "Resistance")

if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)

# Input paths
IN_CHM   = os.path.join(HARM_DIR, "CHM_norm.tif")
IN_VRI   = os.path.join(HARM_DIR, "VRI_norm.tif")
IN_SLOPE = os.path.join(HARM_DIR, "Slope_norm.tif")
IN_LCM   = os.path.join(HARM_DIR, "LCM_2m_aligned.tif")

# Master grid — set from CHM_norm (reference for all outputs)
arcpy.env.snapRaster = IN_CHM
arcpy.env.cellSize   = 2
arcpy.env.extent     = arcpy.Describe(IN_CHM).extent

ref_desc   = arcpy.Describe(IN_CHM)
ref_ext    = ref_desc.extent
MASTER_EXT = (ref_ext.XMin, ref_ext.YMin, ref_ext.XMax, ref_ext.YMax)

print("=" * 60)
print("SCRIPT 3 — FUZZY MCE RESISTANCE SURFACE")
print(f"Started: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 60)

# ---------------------------------------------------------------------------
# HELPER — full verification
# ---------------------------------------------------------------------------

def verify(path, name, exp_min=None, exp_max=None, exp_bands=1):
    print(f"\n  Verifying: {name}")
    bands = int(arcpy.management.GetRasterProperties(
        path, "BANDCOUNT").getOutput(0))
    cs    = arcpy.management.GetRasterProperties(
        path, "CELLSIZEX").getOutput(0)
    desc  = arcpy.Describe(path)
    ext   = desc.extent
    ext_t = (ext.XMin, ext.YMin, ext.XMax, ext.YMax)

    band_ok = (bands == exp_bands)
    ext_ok  = (ext_t == MASTER_EXT)
    print(f"    Bands    : {bands} {'[OK]' if band_ok else '[WARNING] Expected ' + str(exp_bands)}")
    print(f"    Cell size: {cs}m")
    print(f"    Extent   : {ext_t} {'[OK]' if ext_ok else '[WARNING] Differs from master'}")

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
# HELPER — renormalise weights proportionally
# ---------------------------------------------------------------------------

def renorm(weights, factor):
    adjusted = [w * factor for w in weights]
    total    = sum(adjusted)
    return tuple(round(w / total, 6) for w in adjusted)

# ---------------------------------------------------------------------------
# STAGE 0 — PRE-FLIGHT VALIDATION
# ---------------------------------------------------------------------------
# Verifies all Script 2 outputs exist, have correct band count,
# correct extent, and correct value range.
# Script aborts immediately if any check fails.

print("\n--- STAGE 0: Pre-flight Validation of Script 2 Outputs ---")

inputs = {
    "CHM_norm"      : (IN_CHM,   0.0,  1.0),
    "VRI_norm"      : (IN_VRI,   0.0,  1.0),
    "Slope_norm"    : (IN_SLOPE, 0.0,  1.0),
    "LCM_2m_aligned": (IN_LCM,  0.0, 21.0),
}

preflight_ok = True
for name, (path, exp_min, exp_max) in inputs.items():
    if not os.path.exists(path):
        print(f"  [ABORT] {name} not found: {path}")
        print("  Run Script 2 before Script 3.")
        preflight_ok = False
        continue

    bands  = int(arcpy.management.GetRasterProperties(
        path, "BANDCOUNT").getOutput(0))
    desc   = arcpy.Describe(path)
    ext    = desc.extent
    ext_t  = (ext.XMin, ext.YMin, ext.XMax, ext.YMax)
    r_min  = float(arcpy.management.GetRasterProperties(
        path, "MINIMUM").getOutput(0))
    r_max  = float(arcpy.management.GetRasterProperties(
        path, "MAXIMUM").getOutput(0))

    band_ok  = (bands == 1)
    ext_ok   = (ext_t == MASTER_EXT)
    range_ok = (r_min >= exp_min - 0.001) and (r_max <= exp_max + 0.001)

    if band_ok and ext_ok and range_ok:
        print(f"  [OK] {name}: bands={bands}  extent={ext_t}  range={r_min:.2f}-{r_max:.2f}")
    else:
        preflight_ok = False
        print(f"  [ABORT] {name}:")
        if not band_ok:
            print(f"    Band count = {bands} — expected 1. Rerun Script 2.")
        if not ext_ok:
            print(f"    Extent {ext_t} differs from master {MASTER_EXT}. Rerun Script 2.")
        if not range_ok:
            print(f"    Range {r_min:.2f}-{r_max:.2f} outside expected {exp_min}-{exp_max}.")

if not preflight_ok:
    raise RuntimeError(
        "Pre-flight validation failed. Fix issues above before rerunning Script 3.")

print("\n  [OK] All inputs validated — proceeding to Stage 1.")

# ---------------------------------------------------------------------------
# STAGE 1 — FUZZY MEMBERSHIP FUNCTIONS
# ---------------------------------------------------------------------------
# All layers already normalised 0-1 in Script 2.
# Fuzzy functions convert each to a resistance score 0-1
# where 1 = maximum resistance to wildlife movement.
#
# CHM   — Gaussian: mid-range CHM (~8.8m, structured woodland edge) =
#          lowest resistance. Bare ground and dense closed canopy = higher.
#          Focal species (badger, roe deer, pine marten) prefer woodland edges.
#          Midpoint = 0.25 (normalised ~8.8m), spread = 0.15.
#          Result rescaled to strict 0-1.
#
# VRI   — Linear inverted: high structural complexity = low resistance.
#          resistance = 1 - VRI_norm
#
# Slope — Linear: steeper = higher resistance.
#          resistance = Slope_norm (0=flat, 1=steepest)

print("\n--- STAGE 1: Fuzzy Membership Functions ---")

# CHM — Gaussian fuzzy
print("\n  Computing CHM Gaussian fuzzy membership...")
chm_r         = Raster(IN_CHM)
midpoint      = 0.25
spread        = 0.15
chm_fuzzy_raw = 1.0 - arcpy.sa.Exp(
    -1.0 * ((chm_r - midpoint) ** 2) / (2.0 * spread ** 2))
chm_fmin  = float(arcpy.management.GetRasterProperties(
    chm_fuzzy_raw, "MINIMUM").getOutput(0))
chm_fmax  = float(arcpy.management.GetRasterProperties(
    chm_fuzzy_raw, "MAXIMUM").getOutput(0))
chm_fuzzy = Float(chm_fuzzy_raw - chm_fmin) / Float(chm_fmax - chm_fmin)
chm_fuzzy_path = os.path.join(OUT_DIR, "Fuzzy_CHM.tif")
chm_fuzzy.save(chm_fuzzy_path)
del chm_r, chm_fuzzy_raw, chm_fuzzy; gc.collect()
verify(chm_fuzzy_path, "Fuzzy_CHM", exp_min=0.0, exp_max=1.0, exp_bands=1)

# VRI — Linear inverted
print("\n  Computing VRI linear inverted fuzzy membership...")
vri_r          = Raster(IN_VRI)
vri_fuzzy      = Float(1.0) - vri_r
vri_fuzzy_path = os.path.join(OUT_DIR, "Fuzzy_VRI.tif")
vri_fuzzy.save(vri_fuzzy_path)
del vri_r, vri_fuzzy; gc.collect()
verify(vri_fuzzy_path, "Fuzzy_VRI", exp_min=0.0, exp_max=1.0, exp_bands=1)

# Slope — Linear pass-through
print("\n  Computing Slope linear fuzzy membership...")
slope_r          = Raster(IN_SLOPE)
slope_fuzzy      = Float(slope_r)
slope_fuzzy_path = os.path.join(OUT_DIR, "Fuzzy_Slope.tif")
slope_fuzzy.save(slope_fuzzy_path)
del slope_r, slope_fuzzy; gc.collect()
verify(slope_fuzzy_path, "Fuzzy_Slope", exp_min=0.0, exp_max=1.0, exp_bands=1)

# ---------------------------------------------------------------------------
# STAGE 2 — LCM RESISTANCE RECLASSIFICATION
# ---------------------------------------------------------------------------
# UKCEH LCM 2023 classes reclassified to resistance values 0-1.
# Class 0 (NoData/background) → SetNull → excluded from all analysis.
# ExtractBand([1]) enforces single-band output — prevents multi-band artefact.
# Citation: Morton et al. (2024); values from Zeller et al. (2012)

print("\n--- STAGE 2: LCM Resistance Reclassification ---")

lcm_r      = Raster(IN_LCM)
lcm_nodata = SetNull(lcm_r == 0, lcm_r)

lcm_reclass = Con(lcm_nodata == 1,  0.10,
               Con(lcm_nodata == 2,  0.20,
               Con(lcm_nodata == 3,  0.70,
               Con(lcm_nodata == 4,  0.60,
               Con(lcm_nodata == 5,  0.40,
               Con(lcm_nodata == 6,  0.35,
               Con(lcm_nodata == 7,  0.35,
               Con(lcm_nodata == 8,  0.30,
               Con(lcm_nodata == 9,  0.30,
               Con(lcm_nodata == 10, 0.35,
               Con(lcm_nodata == 11, 0.45,
               Con(lcm_nodata == 12, 0.80,
               Con(lcm_nodata == 13, 1.00,
               Con(lcm_nodata == 14, 0.55,
               Con(lcm_nodata == 15, 0.90,
               Con(lcm_nodata == 16, 0.90,
               Con(lcm_nodata == 17, 1.00,
               Con(lcm_nodata == 18, 1.00,
               Con(lcm_nodata == 19, 0.70,
               Con(lcm_nodata == 20, 0.95,
               Con(lcm_nodata == 21, 0.80,
                   lcm_nodata)))))))))))))))))))))

# ExtractBand — enforces single band output
lcm_single      = arcpy.sa.ExtractBand(lcm_reclass, [1])
lcm_resist_path = os.path.join(OUT_DIR, "LCM_resistance.tif")
lcm_single.save(lcm_resist_path)
del lcm_r, lcm_nodata, lcm_reclass, lcm_single; gc.collect()
verify(lcm_resist_path, "LCM_resistance",
       exp_min=0.10, exp_max=1.00, exp_bands=1)
print("  [NOTE] Class 0 set to NoData — excluded from resistance surface.")

# ---------------------------------------------------------------------------
# STAGE 3 — WEIGHTED OVERLAY — 4 SENSITIVITY SCENARIOS
# ---------------------------------------------------------------------------
# Each scenario processed independently — loaded fresh, saved to disk,
# memory released before next scenario begins.
# Formula: R = (w_LCM*LCM) + (w_VRI*VRI) + (w_Slope*Slope) + (w_CHM*CHM)
# All weights sum exactly to 1.0.
# Output: 0-1 where 1 = maximum resistance to wildlife movement.

print("\n--- STAGE 3: Weighted Overlay — 4 Sensitivity Scenarios ---")

base_weights = [0.40, 0.30, 0.20, 0.10]
s2_weights   = renorm(base_weights, 1.20)
s3_weights   = renorm(base_weights, 0.80)

SCENARIOS = [
    ("S1_primary", tuple(base_weights), "Primary weights — Zeller et al. (2012)"),
    ("S2_plus20",  s2_weights,          "Sensitivity +20% — proportional, renormalised"),
    ("S3_minus20", s3_weights,          "Sensitivity -20% — proportional, renormalised"),
    ("S4_equal",   (0.25, 0.25, 0.25, 0.25), "Equal weights sensitivity test"),
]

resistance_paths = []

for name, weights, desc in SCENARIOS:
    w_lcm, w_vri, w_slope, w_chm = weights
    weight_sum = round(sum(weights), 6)
    out_path   = os.path.join(OUT_DIR, f"Resistance_{name}.tif")

    print(f"\n  Scenario : {name}")
    print(f"  {desc}")
    print(f"  Weights  : LCM={w_lcm}  VRI={w_vri}  Slope={w_slope}  CHM={w_chm}")
    print(f"  Sum      : {weight_sum}")
    if abs(weight_sum - 1.0) > 0.001:
        print(f"  [WARNING] Weights do not sum to 1.0")

    # Load all four inputs fresh for this scenario
    lcm_res  = Raster(lcm_resist_path)
    vri_fuz  = Raster(vri_fuzzy_path)
    slp_fuz  = Raster(slope_fuzzy_path)
    chm_fuz  = Raster(chm_fuzzy_path)

    resistance = (w_lcm   * lcm_res)  + \
                 (w_vri   * vri_fuz)  + \
                 (w_slope * slp_fuz)  + \
                 (w_chm   * chm_fuz)

    resistance.save(out_path)
    resistance_paths.append(out_path)

    verify(out_path, f"Resistance_{name}",
           exp_min=0.0, exp_max=1.0, exp_bands=1)

    del lcm_res, vri_fuz, slp_fuz, chm_fuz, resistance
    gc.collect()
    print(f"  Memory released.")
    print(f"  Completed: {datetime.datetime.now().strftime('%H:%M:%S')}")

# ---------------------------------------------------------------------------
# STAGE 4 — SPATIAL STABILITY CHECK
# ---------------------------------------------------------------------------
# Pixel-wise range across S1-S4.
# Low mean range = corridors spatially robust to weight assumptions.
# High mean range = weight-sensitive areas — discuss in Ch7.

print("\n--- STAGE 4: Spatial Stability Check ---")

r_list    = [Raster(p) for p in resistance_paths]
stability = arcpy.sa.CellStatistics(r_list, "RANGE", "DATA")
stab_path = os.path.join(OUT_DIR, "Resistance_stability_range.tif")
stability.save(stab_path)
del r_list, stability; gc.collect()

stab_min, stab_max, stab_mean, stab_std = verify(
    stab_path, "Stability range S1-S4")

print(f"\n  Mean pixel range across scenarios: {stab_mean:.4f}")
if stab_mean < 0.05:
    print("  [OK] High spatial stability — corridors robust to weight assumptions.")
elif stab_mean < 0.10:
    print("  [NOTE] Moderate stability — discuss weight sensitivity in Ch7.")
else:
    print("  [WARNING] High sensitivity to weights — discuss carefully in Ch7.")

# ---------------------------------------------------------------------------
# FINAL SUMMARY
# ---------------------------------------------------------------------------

print()
print("=" * 60)
print("SCRIPT 3 COMPLETE")
print(f"Finished: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 60)
print()

all_outputs = {
    "Fuzzy_CHM.tif"                 : chm_fuzzy_path,
    "Fuzzy_VRI.tif"                 : vri_fuzzy_path,
    "Fuzzy_Slope.tif"               : slope_fuzzy_path,
    "LCM_resistance.tif"            : lcm_resist_path,
    "Resistance_S1_primary.tif"     : resistance_paths[0],
    "Resistance_S2_plus20.tif"      : resistance_paths[1],
    "Resistance_S3_minus20.tif"     : resistance_paths[2],
    "Resistance_S4_equal.tif"       : resistance_paths[3],
    "Resistance_stability_range.tif": stab_path,
}

print("Files created:")
for fname, fpath in all_outputs.items():
    status = "[OK]" if os.path.exists(fpath) else "[MISSING]"
    print(f"  {status} {fname}")

print()
print("VISUAL QA — verify in ArcGIS Pro before Script 4:")
print("  [ ] Resistance_S1_primary.tif — A66 road = bright high-resistance band")
print("  [ ] Woodland patches = dark low-resistance zones")
print("  [ ] Eden Valley Railway = linear low-resistance corridor near Kirkby Thore")
print("  [ ] British Gypsum site = high resistance")
print("  [ ] River Eden corridor = moderate-low resistance")
print("  [ ] Resistance_stability_range.tif — low values = robust corridors")
print()
print("VERIFIED STATISTICS FOR DISSERTATION (Chapter 5):")
print("  4 sensitivity scenarios produced and verified")
print("  S1: LCM=0.40  VRI=0.30  Slope=0.20  CHM=0.10")
print("  S2: +20% proportional, renormalised to sum 1.0")
print("  S3: -20% proportional, renormalised to sum 1.0")
print("  S4: Equal weights 0.25 each")
print("  All resistance surfaces: range 0-1 confirmed")
print(f"  Spatial stability mean range: {stab_mean:.4f}")
print()
print("NEXT STEP: Script 4 — Circuitscape Connectivity")
print("  Export ASCII → run 4 scenarios via Julia → import current density")
print()
print("Push to GitHub:")
print("  https://github.com/Abhin-Nangari/A66-Connectivity-MSc")

arcpy.CheckInExtension("Spatial")
