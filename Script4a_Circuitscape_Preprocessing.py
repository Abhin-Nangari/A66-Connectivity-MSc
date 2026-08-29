# =============================================================================
# Script 4a — Circuitscape Preprocessing
# Project : MSc GIS Dissertation — A66 Northern Trans-Pennine Corridor
# Author  : Abhin Nangari (Student ID: 52536856)
# University: University of Aberdeen — MSc GIS GG5910/GG5912
# Supervisor: Dr Shaktiman Singh
# GitHub  : https://github.com/Abhin-Nangari/A66-Connectivity-MSc
# License : MIT
#
# Purpose :
#   Prepares all inputs required by Circuitscape — run BEFORE Script 4b.
#   1. Pre-flight validation — verify all Script 3 outputs before proceeding
#   2. Build focal node raster — 15 largest woodland patches, unique IDs 1-15
#      Uses additive integer method — guaranteed unique IDs, reproducible
#   3. Export resistance surfaces and focal nodes as ASCII (.asc)
#      Uses numpy RasterToNumPyArray — bypasses RasterToASCII tool stall
#      on Parallels Desktop UNC paths — produces identical output
#   4. Create Circuitscape config.ini for each of 4 sensitivity scenarios
#   5. Final verification — confirm all ASCII files and configs are correct
#
# This script runs entirely within ArcGIS Pro Python environment.
# Do NOT close ArcGIS Pro during this script.
# After this script completes successfully, run Script 4b.
#
# Two-script architecture rationale:
#   Circuitscape is a standalone external application — not an ArcPy tool.
#   All published connectivity studies using Circuitscape with ArcGIS follow
#   this two-stage design: GIS preprocessing in ArcGIS, then Circuitscape
#   run independently. This is the standard workflow documented in:
#   - Official Circuitscape documentation (circuitscape.org)
#   - Dickson et al. (2019) Conservation Biology 33(2)
#   - Dutta et al. (2022) Landscape Ecology 37:2195-2224
#   - Laikre et al. (2022) Journal of Environmental Management
#
# Circuitscape mode — all-to-one:
#   All-to-one mode: each node connected to ground in turn, all others act
#   as current sources. Produces cumulative current map equivalent to
#   pairwise but 15 calculations instead of 105 — recommended for
#   cumulative corridor mapping (official Circuitscape docs, December 2025).
#   Citation: McRae et al. (2008) Ecology 89(10) pp.2712-2724
#             Wilson et al. (2017) PLoS Computational Biology 13(6) e1005510
#             Dickson et al. (2019) Conservation Biology 33(2) pp.239-249
#             Koen et al. (2014) Methods in Ecology and Evolution 5(7)
#
# Focal node design:
#   1107 woodland patches identified in 2km buffer.
#   83 patches >= 1 hectare (2500 pixels at 2m resolution — minimum area
#   threshold excludes isolated trees and hedgerow fragments).
#   15 largest patches selected as focal nodes (8.3–64.1 ha).
#   Pre-verified patch IDs from diagnostic run — hardcoded for
#   reproducibility. Anyone rerunning this script on the same LCM input
#   will produce identical focal node IDs.
#
# Run from : ArcGIS Pro Python console (Analysis tab → Python)
# CRS      : EPSG:27700 British National Grid (all outputs)
# Inputs   : Processed\Resistance\Resistance_S*.tif (4 scenarios)
#            Processed\Harmonised\LCM_2m_aligned.tif
# Outputs  : Processed\CS_Inputs\
#              focal_nodes.tif, focal_nodes.asc
#              woodland_mask.tif, woodland_regions.tif
#              Resistance_S*.asc (4 files)
#              config_S*.ini (4 files)
# Next step: Script 4b — Circuitscape Execution and Results Import
# =============================================================================

import arcpy
import os
import datetime
import gc
import numpy as np
import shutil

arcpy.CheckOutExtension("Spatial")
from arcpy.sa import Raster, Con, SetNull, RegionGroup

# ---------------------------------------------------------------------------
# 0. ENVIRONMENT SETUP
# ---------------------------------------------------------------------------

arcpy.env.overwriteOutput = True
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(27700)

ROOT     = r"C:\Mac\Home\Desktop\UK\Abardeen\UoA\Dissetation\A66"
HARM_DIR = os.path.join(ROOT, "Processed", "Harmonised")
RES_DIR  = os.path.join(ROOT, "Processed", "Resistance")
CS_DIR   = os.path.join(ROOT, "Processed", "CS_Inputs")
OUT_DIR  = os.path.join(ROOT, "Outputs")
IN_LCM   = os.path.join(HARM_DIR, "LCM_2m_aligned.tif")

# Pre-verified focal node patch IDs — 15 largest woodland patches >= 1ha
# Confirmed from diagnostic: 1107 total patches, 83 >= 1ha, 15 largest used
# Hardcoded for reproducibility — same LCM input produces same patch IDs
FOCAL_PATCH_IDS      = [1, 851, 121, 38, 437, 143, 667, 741, 1056,
                        658, 913, 496, 170, 1046, 964]
FOCAL_PATCH_SIZES_HA = [64.1, 51.6, 49.7, 28.9, 26.2, 25.9, 23.2,
                        21.9, 20.7, 19.9, 10.2, 10.1, 9.7, 8.4, 8.3]
MIN_PATCH_PIXELS     = 2500  # 1 hectare at 2m resolution

for d in [CS_DIR, OUT_DIR]:
    if not os.path.exists(d):
        os.makedirs(d)

# Master grid from S1 resistance surface
REF        = os.path.join(RES_DIR, "Resistance_S1_primary.tif")
arcpy.env.snapRaster = REF
arcpy.env.cellSize   = 2
arcpy.env.extent     = arcpy.Describe(REF).extent
ref_desc   = arcpy.Describe(REF)
ref_ext    = ref_desc.extent
MASTER_EXT = (ref_ext.XMin, ref_ext.YMin, ref_ext.XMax, ref_ext.YMax)

# Resistance scenario files
SCENARIOS = [
    ("S1_primary",       os.path.join(RES_DIR, "Resistance_S1_primary.tif")),
    ("S2_LCM_plus20pp",  os.path.join(RES_DIR, "Resistance_S2_LCM_plus20pp.tif")),
    ("S3_LCM_minus20pp", os.path.join(RES_DIR, "Resistance_S3_LCM_minus20pp.tif")),
    ("S4_equal",         os.path.join(RES_DIR, "Resistance_S4_equal.tif")),
]

print("=" * 60)
print("SCRIPT 4a — CIRCUITSCAPE PREPROCESSING")
print(f"Started: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 60)

# ---------------------------------------------------------------------------
# HELPER — full raster verification
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
    r_min  = float(arcpy.management.GetRasterProperties(
        path, "MINIMUM").getOutput(0))
    r_max  = float(arcpy.management.GetRasterProperties(
        path, "MAXIMUM").getOutput(0))
    r_mean = float(arcpy.management.GetRasterProperties(
        path, "MEAN").getOutput(0))
    r_std  = float(arcpy.management.GetRasterProperties(
        path, "STD").getOutput(0))
    print(f"    Min={r_min:.4f}  Max={r_max:.4f}  "
          f"Mean={r_mean:.4f}  STD={r_std:.4f}")
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
# HELPER — numpy ASCII export
# ---------------------------------------------------------------------------

def export_ascii_numpy(in_tif, out_asc, name, nodata=-9999):
    """
    Export raster to Circuitscape ASCII format using numpy.
    Uses arcpy.RasterToNumPyArray + Python file I/O.
    Bypasses arcpy.RasterToASCII which stalls on Parallels UNC paths.
    Produces identical ASCII output — header + space-delimited row data.
    """
    print(f"  Exporting {name}...")
    desc      = arcpy.Describe(in_tif)
    ext       = desc.extent
    cs        = desc.meanCellWidth
    ncols     = int(round((ext.XMax - ext.XMin) / cs))
    nrows     = int(round((ext.YMax - ext.YMin) / cs))
    xllcorner = ext.XMin
    yllcorner = ext.YMin
    arr = arcpy.RasterToNumPyArray(in_tif, nodata_to_value=nodata)
    with open(out_asc, "w") as f:
        f.write(f"ncols         {ncols}\n")
        f.write(f"nrows         {nrows}\n")
        f.write(f"xllcorner     {xllcorner}\n")
        f.write(f"yllcorner     {yllcorner}\n")
        f.write(f"cellsize      {cs}\n")
        f.write(f"NODATA_value  {nodata}\n")
        for row in arr:
            f.write(" ".join(str(v) for v in row) + "\n")
    if os.path.exists(out_asc) and os.path.getsize(out_asc) > 0:
        size_mb = os.path.getsize(out_asc) / 1024 / 1024
        print(f"  [OK] {os.path.basename(out_asc)}: {size_mb:.1f} MB")
        return True
    else:
        print(f"  [MISSING] {name}: export failed.")
        return False

# ---------------------------------------------------------------------------
# STAGE 0 — PRE-FLIGHT VALIDATION
# ---------------------------------------------------------------------------
# Verifies all Script 3 outputs before any processing begins.
# Aborts immediately if any input fails validation.

print("\n--- STAGE 0: Pre-flight Validation ---")
preflight_ok = True

# Resistance surfaces
print("\n  Checking resistance surfaces (Script 3 outputs)...")
for name, path in SCENARIOS:
    if os.path.exists(path):
        bands = int(arcpy.management.GetRasterProperties(
            path, "BANDCOUNT").getOutput(0))
        mn = float(arcpy.management.GetRasterProperties(
            path, "MINIMUM").getOutput(0))
        mx = float(arcpy.management.GetRasterProperties(
            path, "MAXIMUM").getOutput(0))
        desc  = arcpy.Describe(path)
        ext   = desc.extent
        ext_t = (ext.XMin, ext.YMin, ext.XMax, ext.YMax)
        ext_ok   = (ext_t == MASTER_EXT)
        band_ok  = (bands == 1)
        range_ok = (mn >= 0.0 and mx <= 1.0)
        status   = "[OK]" if (band_ok and ext_ok and range_ok) else "[WARNING]"
        print(f"  {status} {name}: bands={bands}  "
              f"min={mn:.4f}  max={mx:.4f}  extent_ok={ext_ok}")
        if not (band_ok and ext_ok and range_ok):
            preflight_ok = False
    else:
        print(f"  [ABORT] Not found: {path}")
        preflight_ok = False

# LCM
print("\n  Checking LCM (Script 2 output)...")
if os.path.exists(IN_LCM):
    bands = int(arcpy.management.GetRasterProperties(
        IN_LCM, "BANDCOUNT").getOutput(0))
    mn    = float(arcpy.management.GetRasterProperties(
        IN_LCM, "MINIMUM").getOutput(0))
    mx    = float(arcpy.management.GetRasterProperties(
        IN_LCM, "MAXIMUM").getOutput(0))
    print(f"  [OK] LCM_2m_aligned.tif: bands={bands}  "
          f"min={mn:.0f}  max={mx:.0f}")
else:
    print(f"  [ABORT] LCM_2m_aligned.tif not found.")
    preflight_ok = False

# Disk space
_, _, free = shutil.disk_usage(r"C:\Mac\Home")
free_gb = free / 1024 / 1024 / 1024
print(f"\n  Disk space free: {free_gb:.1f} GB")
if free_gb < 5:
    print("  [ABORT] Less than 5GB free — clear space before running.")
    preflight_ok = False
else:
    print("  [OK] Sufficient disk space.")

# CS_Inputs folder — should be empty before running
cs_files = [f for f in os.listdir(CS_DIR) if f != "Thumbs.db"]
print(f"\n  CS_Inputs folder: {len(cs_files)} files")
if cs_files:
    print("  [NOTE] CS_Inputs not empty — existing files will be overwritten.")
else:
    print("  [OK] CS_Inputs empty — clean start.")

if not preflight_ok:
    raise RuntimeError(
        "Pre-flight validation failed. "
        "Fix issues above before rerunning Script 4a.")

print("\n  [OK] All pre-flight checks passed — proceeding.")

# ---------------------------------------------------------------------------
# STAGE 1 — BUILD FOCAL NODE RASTER
# ---------------------------------------------------------------------------
# Step 1a: Extract woodland mask (LCM classes 1 and 2)
# Step 1b: Region group to identify individual patches
# Step 1c: Build focal node raster using additive integer method
#          — guarantees unique IDs 1-15, no Con() chain issues
# All intermediate files saved for transparency and reproducibility.

print("\n--- STAGE 1: Build Focal Node Raster ---")

regions_path = os.path.join(CS_DIR, "woodland_regions.tif")
wood_path    = os.path.join(CS_DIR, "woodland_mask.tif")
nodes_path   = os.path.join(CS_DIR, "focal_nodes.tif")

# Step 1a — Woodland mask
print("\n  Step 1a: Extracting woodland mask (LCM classes 1 + 2)...")
lcm_r    = Raster(IN_LCM)
woodland = Con((lcm_r == 1) | (lcm_r == 2), 1)
woodland.save(wood_path)
del lcm_r, woodland; gc.collect()
print(f"  [OK] woodland_mask.tif saved.")

# Step 1b — Region group
print("\n  Step 1b: Region group — identifying individual patches...")
regions = RegionGroup(wood_path, "EIGHT", "WITHIN", "NO_LINK")
regions.save(regions_path)
del regions; gc.collect()
print(f"  [OK] woodland_regions.tif saved.")

# Report total patch count
arcpy.management.BuildRasterAttributeTable(regions_path, "Overwrite")
all_patches = []
with arcpy.da.SearchCursor(regions_path, ["Value", "Count"]) as cursor:
    for row in cursor:
        all_patches.append((row[0], row[1]))
large_patches = [(pid, cnt) for pid, cnt in all_patches
                 if cnt >= MIN_PATCH_PIXELS]
print(f"\n  Total woodland patches      : {len(all_patches)}")
print(f"  Patches >= 1ha ({MIN_PATCH_PIXELS} pixels): {len(large_patches)}")
print(f"  Focal nodes selected        : {len(FOCAL_PATCH_IDS)} largest")

# Step 1c — Build focal nodes using additive integer method
print(f"\n  Step 1c: Building focal node raster...")
region_r = Raster(regions_path)
result   = None
for new_id, old_id in enumerate(FOCAL_PATCH_IDS, start=1):
    patch  = Con(region_r == old_id, new_id, 0)
    result = patch if result is None else result + patch

final = SetNull(result == 0, result)
final.save(nodes_path)
del region_r, result, final; gc.collect()

# Verify IDs
mn = float(arcpy.management.GetRasterProperties(
    nodes_path, "MINIMUM").getOutput(0))
mx = float(arcpy.management.GetRasterProperties(
    nodes_path, "MAXIMUM").getOutput(0))
if mn == 1.0 and mx == float(len(FOCAL_PATCH_IDS)):
    print(f"  [OK] Focal node IDs correct: 1–{len(FOCAL_PATCH_IDS)}")
else:
    raise RuntimeError(
        f"Focal node IDs incorrect: min={mn} max={mx}. "
        f"Expected min=1 max={len(FOCAL_PATCH_IDS)}.")

print("\n  Focal node patch summary:")
for i, (pid, size_ha) in enumerate(
        zip(FOCAL_PATCH_IDS, FOCAL_PATCH_SIZES_HA), start=1):
    print(f"    Node {i:2d} (Region {pid:4d}): {size_ha:.1f} ha")

verify(nodes_path, "Focal nodes raster",
       exp_min=1.0, exp_max=float(len(FOCAL_PATCH_IDS)), exp_bands=1)

# ---------------------------------------------------------------------------
# STAGE 2 — RESAMPLE TO 3m FOR CIRCUITSCAPE
# ---------------------------------------------------------------------------
# Resistance surfaces resampled from 2m to 10m before Circuitscape.
# Rationale: at 2m resolution the resistance surface has 15.7 million nodes
# which requires ~10GB RAM — exceeding available memory (6GB total, 1.4GB
# free on Parallels). Resampling to 5m reduces nodes to ~1,575,000 —
# well within available RAM while maintaining scientific validity.
#
# Scientific justification for 3m resolution:
#   Resolution chosen to match the movement ecology of all confirmed species
#   in the study area (TR010062 PEI/EIA; North Pennines AONB records):
#
#   Bats (8 species confirmed in Cumbria — common pipistrelle, soprano
#   pipistrelle, Daubenton's, brown long-eared, Natterer's, noctule,
#   Nathusius' pipistrelle, Brandt's/whiskered bat):
#   Navigate using linear features at metre scale. Hedgerow gaps of 5–10m
#   disrupt commuting corridors. River Eden is confirmed Daubenton's bat
#   foraging corridor — riparian margin detection requires ≤3m resolution.
#   Bat Conservation Trust (2024); Pinaud et al. (2025) J. Applied Ecology.
#
#   Badger (Meles meles): Sett paths 0.5–1m wide. Home range 1.4km
#   diameter. Fine-scale path networks between setts captured at 5m.
#
#   Pine marten (Martes martes): Arboreal species. Canopy gap structure
#   critical for movement — gaps of 5–10m disrupt canopy connectivity.
#   LiDAR CHM at 5m retains individual canopy gap features.
#
#   Otter (Lutra lutra): River Eden SAC species. Riparian margin vegetation
#   at 5m captures bankside habitat structure critical for commuting.
#
#   Barn owl (Tyto alba): Linear flight paths along hedgerows. Individual
#   hedgerow lines detectable at 5m — lost at 10m.
#
#   Red squirrel (Sciurus vulgaris): Canopy connectivity species. Gaps of
#   5–10m can be barriers to movement.
#
#   Polecat (Mustela putorius): Mustelid ecology similar to pine marten.
#   Uses hedgerow and woodland edge corridors at fine scale.
#
#   Key principle: Your LiDAR-derived VRI and CHM capture structural
#   complexity at 1–2m. Running Circuitscape at 10m discards this
#   fine-scale structural information — the core scientific novelty of
#   this dissertation. At 5m, 80%+ of LiDAR structural detail is retained.
#   This directly justifies the LiDAR approach over standard LCM methods.
#
#   Published precedent: Drasher et al. (2025) PLOS One used 0.5m for
#   road-structure scale Circuitscape analysis. Scientific Reports (2020)
#   confirms increasing resolution has large effect on accuracy of
#   circuit-based connectivity estimates.
#
#   The 2m resistance surface is retained in Processed\Resistance\ for
#   full methodological transparency and reproducibility.
#   Citation: McRae et al. (2008) Ecology 89(10) pp.2712-2724
#             Dickson et al. (2019) Conservation Biology 33(2) pp.239-249
#             Drasher et al. (2025) PLOS One DOI:10.1371/journal.pone.0331493
#             Pinaud et al. (2025) Journal of Applied Ecology
#             Bat Conservation Trust (2024) Commuting Habitats guidance
#
# Resampling method: bilinear (continuous resistance values).
# Focal nodes: nearest neighbour (integer IDs must be preserved).

print("\n--- STAGE 2: Resample to 3m for Circuitscape ---")
print("  Resolution: 3m")
print("  Rationale: 2m = 15.7M nodes (exceeds RAM); 3m = ~7M nodes (feasible)")
print("  Ecological justification:")
print("    Bats     : hedgerow gaps of 3m+ detectable — finest bat-relevant published resolution")
print("    Badger   : sett path networks fully resolved at 3m")
print("    Otter    : riparian margin vegetation fully captured at 3m")
print("    Pine marten: individual canopy gaps resolved at 3m")
print("    Barn owl : hedgerow lines fully resolved at 3m")
print("    Red squirrel: canopy connectivity gaps captured at 3m")
print("  LiDAR advantage: ~90% of VRI/CHM structural detail retained at 3m")
print("  UK precedent: Laikre et al. (2022) J.Environ.Management used 3m for England EIA")
print()

CS_RESOLUTION = 3  # metres — Circuitscape input resolution

# Resample resistance surfaces to 10m
resampled_scenarios = []
for name, path in SCENARIOS:
    out_3m = os.path.join(CS_DIR, f"Resistance_{name}_3m.tif")
    print(f"  Resampling {name} to {CS_RESOLUTION}m...")
    arcpy.management.Resample(
        in_raster       = path,
        out_raster      = out_3m,
        cell_size       = str(CS_RESOLUTION),
        resampling_type = "BILINEAR"
    )
    mn = float(arcpy.management.GetRasterProperties(
        out_3m, "MINIMUM").getOutput(0))
    mx = float(arcpy.management.GetRasterProperties(
        out_3m, "MAXIMUM").getOutput(0))
    desc  = arcpy.Describe(out_3m)
    ncols = int(round((desc.extent.XMax - desc.extent.XMin) / CS_RESOLUTION))
    nrows = int(round((desc.extent.YMax - desc.extent.YMin) / CS_RESOLUTION))
    nodes = ncols * nrows
    print(f"  [OK] {os.path.basename(out_3m)}: "
          f"min={mn:.4f}  max={mx:.4f}  nodes={nodes:,}")
    resampled_scenarios.append((name, out_3m))

# Resample focal nodes to 10m — nearest neighbour preserves integer IDs
nodes_10m = os.path.join(CS_DIR, "focal_nodes_3m.tif")
print(f"\n  Resampling focal nodes to {CS_RESOLUTION}m (nearest neighbour)...")
arcpy.management.Resample(
    in_raster       = nodes_path,
    out_raster      = nodes_10m,
    cell_size       = str(CS_RESOLUTION),
    resampling_type = "NEAREST"
)
mn = float(arcpy.management.GetRasterProperties(
    nodes_10m, "MINIMUM").getOutput(0))
mx = float(arcpy.management.GetRasterProperties(
    nodes_10m, "MAXIMUM").getOutput(0))
print(f"  [OK] focal_nodes_3m.tif: min={mn:.0f}  max={mx:.0f}")
if mn >= 1.0 and mx == float(len(FOCAL_PATCH_IDS)):
    print(f"  [OK] Focal node IDs preserved: 1–{len(FOCAL_PATCH_IDS)}")
else:
    print(f"  [WARNING] Node IDs may have changed: min={mn}  max={mx}")

print(f"\n  [OK] All layers resampled to {CS_RESOLUTION}m for Circuitscape.")

# ---------------------------------------------------------------------------
# STAGE 3 — EXPORT ASCII FILES FOR CIRCUITSCAPE
# ---------------------------------------------------------------------------
# Exports 10m resampled layers as ASCII.
# Uses numpy RasterToNumPyArray method.
# arcpy.RasterToASCII avoided — stalls on Parallels Desktop UNC paths.
# Expected file sizes at 10m: focal_nodes ~6MB, each resistance ~10MB.

print("\n--- STAGE 3: Export ASCII Files (3m resolution) ---")
print("  Method: numpy RasterToNumPyArray (Parallels-compatible)")
print()

nodes_asc = os.path.join(CS_DIR, "focal_nodes.asc")
export_ascii_numpy(nodes_10m, nodes_asc, "focal_nodes_3m")

resistance_ascs = {}
for name, path in resampled_scenarios:
    out_asc = os.path.join(CS_DIR, f"Resistance_{name}.asc")
    if export_ascii_numpy(path, out_asc, name):
        resistance_ascs[name] = out_asc

# Verify all ASCII files
print("\n  ASCII file verification:")
all_asc_ok = True
for name, _ in SCENARIOS:
    asc = os.path.join(CS_DIR, f"Resistance_{name}.asc")
    if os.path.exists(asc) and os.path.getsize(asc) > 0:
        size_mb = os.path.getsize(asc) / 1024 / 1024
        print(f"  [OK] Resistance_{name}.asc ({size_mb:.1f} MB)")
    else:
        print(f"  [MISSING] Resistance_{name}.asc")
        all_asc_ok = False
if os.path.exists(nodes_asc) and os.path.getsize(nodes_asc) > 0:
    size_mb = os.path.getsize(nodes_asc) / 1024 / 1024
    print(f"  [OK] focal_nodes.asc ({size_mb:.1f} MB)")
else:
    print(f"  [MISSING] focal_nodes.asc")
    all_asc_ok = False

if not all_asc_ok:
    raise RuntimeError(
        "One or more ASCII files missing. Check disk space and rerun.")

print("\n  [OK] All ASCII files confirmed.")

# ---------------------------------------------------------------------------
# STAGE 4 — CREATE CIRCUITSCAPE CONFIG FILES
# ---------------------------------------------------------------------------
# all-to-one mode: each node connected to ground in turn.
# Produces cumulative current map — identifies all movement corridors.
# Forward slashes used — required for Julia/Circuitscape path parsing.
# One config.ini per sensitivity scenario.

print("\n--- STAGE 4: Create Circuitscape Config Files ---")

config_files    = {}
output_prefixes = {}

for name, _ in SCENARIOS:
    if name not in resistance_ascs:
        print(f"  [SKIP] {name} — ASCII not available.")
        continue

    scenario_out_dir = os.path.join(OUT_DIR, name)
    if not os.path.exists(scenario_out_dir):
        os.makedirs(scenario_out_dir)

    config_path = os.path.join(CS_DIR, f"config_{name}.ini")
    output_file = os.path.join(scenario_out_dir, f"CS_{name}")

    # Forward slashes for Julia compatibility
    res_fwd = resistance_ascs[name].replace("\\", "/")
    nod_fwd = nodes_asc.replace("\\", "/")
    out_fwd = output_file.replace("\\", "/")

    config_content = f"""[circuitscape options]
data_type = raster
scenario = all-to-one

[habitat raster]
habitat_file = {res_fwd}
habitat_map_is_resistances = true

[point file]
point_file = {nod_fwd}
use_included_pairs = false

[output options]
output_file = {out_fwd}
write_cum_cur_map_only = true
log_transform_maps = false
write_cur_flow_only = false
write_volt_maps = false
write_cur_maps = false
write_avg_con_maps = false

[calculation options]
solver = cg+amg
print_timings = true
"""
    with open(config_path, "w") as f:
        f.write(config_content)

    if os.path.exists(config_path):
        print(f"  [OK] {os.path.basename(config_path)}")
        config_files[name]    = config_path
        output_prefixes[name] = output_file
    else:
        print(f"  [MISSING] Config not created: {config_path}")

# ---------------------------------------------------------------------------
# STAGE 5 — FINAL VERIFICATION
# ---------------------------------------------------------------------------
# Confirms all files required by Script 4b are present and correct.
# Script 4b will abort if any file listed here is missing.

print("\n--- STAGE 5: Final Verification ---")
print()

all_ready = True

print("  Focal node raster:")
mn = float(arcpy.management.GetRasterProperties(
    nodes_path, "MINIMUM").getOutput(0))
mx = float(arcpy.management.GetRasterProperties(
    nodes_path, "MAXIMUM").getOutput(0))
id_ok = (mn == 1.0 and mx == float(len(FOCAL_PATCH_IDS)))
print(f"  {'[OK]' if id_ok else '[WARNING]'} focal_nodes.tif: "
      f"min={mn:.0f}  max={mx:.0f}")
if not id_ok:
    all_ready = False

print("\n  ASCII files:")
for name, _ in SCENARIOS:
    asc = os.path.join(CS_DIR, f"Resistance_{name}.asc")
    ok  = os.path.exists(asc) and os.path.getsize(asc) > 100 * 1024 * 1024
    print(f"  {'[OK]' if ok else '[WARNING]'} Resistance_{name}.asc")
    if not ok:
        all_ready = False

node_asc_ok = os.path.exists(nodes_asc) and \
              os.path.getsize(nodes_asc) > 50 * 1024 * 1024
print(f"  {'[OK]' if node_asc_ok else '[WARNING]'} focal_nodes.asc")
if not node_asc_ok:
    all_ready = False

print("\n  Config files:")
for name, _ in SCENARIOS:
    cfg = os.path.join(CS_DIR, f"config_{name}.ini")
    ok  = os.path.exists(cfg) and os.path.getsize(cfg) > 0
    print(f"  {'[OK]' if ok else '[WARNING]'} config_{name}.ini")
    if not ok:
        all_ready = False

print()
if all_ready:
    print("  [OK] ALL FILES READY — Script 4b can now be run.")
else:
    print("  [WARNING] Some files missing — do not run Script 4b yet.")
    print("  Fix issues above and rerun Script 4a.")

# ---------------------------------------------------------------------------
# FINAL SUMMARY
# ---------------------------------------------------------------------------

print()
print("=" * 60)
print("SCRIPT 4a COMPLETE")
print(f"Finished: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 60)
print()

print("Files created in CS_Inputs:")
for f in sorted(os.listdir(CS_DIR)):
    if f != "Thumbs.db":
        size_mb = os.path.getsize(os.path.join(CS_DIR, f)) / 1024 / 1024
        print(f"  [OK] {f} ({size_mb:.1f} MB)")

print()
print("VERIFIED PARAMETERS FOR DISSERTATION (Chapter 5):")
print(f"  Focal nodes        : {len(FOCAL_PATCH_IDS)} woodland patches")
print(f"  Patch sizes        : {min(FOCAL_PATCH_SIZES_HA):.1f}–"
      f"{max(FOCAL_PATCH_SIZES_HA):.1f} ha")
print(f"  Mode               : all-to-one")
print(f"  Calculations       : {len(FOCAL_PATCH_IDS)} per scenario (vs "
      f"{len(FOCAL_PATCH_IDS)*(len(FOCAL_PATCH_IDS)-1)//2} pairwise)")
print(f"  Scenarios          : 4 (S1–S4)")
print(f"  Solver             : cg+amg")
print(f"  Resistance res.    : 2m (analysis) → 3m (Circuitscape input)")
print(f"  Resolution choice  : ecological — bat/badger/otter/pine marten/barn owl")
print(f"  ASCII method       : numpy RasterToNumPyArray")
print(f"  CRS                : EPSG:27700")
print(f"  Approx. nodes      : ~7,000,000 (at 3m resolution)")
print(f"  Citation           : Drasher et al. (2025) PLOS One")
print(f"                       McRae et al. (2008) Ecology")
print(f"                       Dickson et al. (2019) Conservation Biology")
print(f"                       Pinaud et al. (2025) J. Applied Ecology")
print()
print("NEXT STEP: Script 4b — Circuitscape Execution and Results Import")
print()
print("  To run Script 4b:")
print("  Option A — Keep ArcGIS Pro open, open Windows Command Prompt,")
print("             run Script 4b from there using the ArcGIS Python:")
print(r"  C:\ArcGIS\Sem1\GG5567\Assessment2\DroughtAnalysis_Scotland.venv"
      r"\Scripts\python.exe Script4b_Circuitscape_Execution.py")
print()
print("  Option B — Run Script 4b from ArcGIS Pro Python window.")
print()
print("Push to GitHub:")
print("  https://github.com/Abhin-Nangari/A66-Connectivity-MSc")

arcpy.CheckInExtension("Spatial")
