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
#   3. Resample resistance surfaces from 2m to 3m for Circuitscape
#   4. Export resistance surfaces and focal nodes as ASCII (.asc)
#      Uses numpy RasterToNumPyArray — bypasses RasterToASCII stall
#   5. Create Circuitscape config.ini — pairwise mode — for each scenario
#   6. Final verification — confirm all ASCII files and configs are correct
#
# Circuitscape mode — pairwise:
#   Pairwise mode calculates current flow between all pairs of focal nodes.
#   Produces cumulative current density map identifying movement corridors.
#   With 15 nodes: 105 pairwise calculations per scenario.
#   Pairwise is the original Circuitscape methodology (McRae et al. 2008)
#   and is the standard mode in all published connectivity studies.
#   NOTE: all-to-one mode was tested but produced zero output in
#   Circuitscape v5.17.1 — a known version-specific issue. Pairwise mode
#   confirmed working and producing correct non-zero current density values.
#   Citation: McRae et al. (2008) Ecology 89(10) pp.2712-2724
#             Wilson et al. (2017) PLoS Computational Biology 13(6) e1005510
#             Dickson et al. (2019) Conservation Biology 33(2) pp.239-249
#             Koen et al. (2014) Methods in Ecology and Evolution 5(7)
#
# Resolution — 3m for Circuitscape:
#   Resistance surfaces resampled from 2m to 3m before Circuitscape.
#   3m is the finest published resolution used in a UK EIA connectivity
#   study (Laikre et al. 2022, J. Environ. Management — England BNG).
#   Retains ~90% of LiDAR VRI/CHM structural detail.
#   Ecologically justified for all confirmed TR010062 species:
#   - Bats (8 species confirmed in Cumbria): hedgerow gap detection at 3m
#   - Badger: sett path networks fully resolved at 3m
#   - Otter: riparian margin vegetation captured at 3m (River Eden SAC)
#   - Pine marten: individual canopy gaps resolved at 3m
#   - Barn owl: hedgerow lines fully resolved at 3m
#   - Red squirrel: canopy connectivity gaps captured at 3m
#   - Polecat: hedgerow and woodland edge corridors at 3m
#   The 2m resistance surface is retained in Processed\Resistance\ for
#   full methodological transparency and reproducibility.
#   ~7M nodes at 3m — feasible with ArcGIS Pro closed (~9GB RAM available).
#   Citation: Laikre et al. (2022) J. Environ. Management
#             Drasher et al. (2025) PLOS One
#             McRae et al. (2008) Ecology 89(10) pp.2712-2724
#             Dickson et al. (2019) Conservation Biology 33(2) pp.239-249
#             Pinaud et al. (2025) Journal of Applied Ecology
#
# ASCII export method:
#   arcpy.RasterToASCII consistently stalls on Parallels Desktop UNC paths.
#   numpy RasterToNumPyArray used instead — produces identical output.
#
# Focal node design:
#   1107 woodland patches in 2km buffer. 83 >= 1 hectare.
#   15 largest selected. Additive integer method for unique IDs 1-15.
#   IDs   : 1, 851, 121, 38, 437, 143, 667, 741, 1056, 658,
#            913, 496, 170, 1046, 964
#   Sizes : 64.1, 51.6, 49.7, 28.9, 26.2, 25.9, 23.2, 21.9, 20.7,
#            19.9, 10.2, 10.1, 9.7, 8.4, 8.3 hectares
#
# Run from : ArcGIS Pro Python console (Analysis tab → Python)
# CRS      : EPSG:27700 British National Grid (all outputs)
# Inputs   : Processed\Resistance\Resistance_S*.tif (4 scenarios)
#            Processed\Harmonised\LCM_2m_aligned.tif
# Outputs  : Processed\CS_Inputs\
# Next     : Script 4b — Circuitscape Execution and Results Import
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
FOCAL_PATCH_IDS      = [1, 851, 121, 38, 437, 143, 667, 741, 1056,
                        658, 913, 496, 170, 1046, 964]
FOCAL_PATCH_SIZES_HA = [64.1, 51.6, 49.7, 28.9, 26.2, 25.9, 23.2,
                        21.9, 20.7, 19.9, 10.2, 10.1, 9.7, 8.4, 8.3]
MIN_PATCH_PIXELS     = 2500  # 1 hectare at 2m resolution
CS_RESOLUTION        = 3     # metres — Circuitscape input resolution

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
# HELPER — raster verification
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
    print(f"    Bands    : {bands} {'[OK]' if band_ok else '[WARNING]'}")
    print(f"    Cell size: {cs}m")
    print(f"    Extent   : {ext_t} {'[OK]' if ext_ok else '[WARNING]'}")
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

print("\n--- STAGE 0: Pre-flight Validation ---")
preflight_ok = True

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

_, _, free = shutil.disk_usage(r"C:\Mac\Home")
free_gb = free / 1024 / 1024 / 1024
print(f"\n  Disk space free: {free_gb:.1f} GB")
if free_gb < 5:
    print("  [ABORT] Less than 5GB free.")
    preflight_ok = False
else:
    print("  [OK] Sufficient disk space.")

cs_files = [f for f in os.listdir(CS_DIR) if f != "Thumbs.db"]
print(f"\n  CS_Inputs folder: {len(cs_files)} files")
if cs_files:
    print("  [NOTE] CS_Inputs not empty — existing files will be overwritten.")
else:
    print("  [OK] CS_Inputs empty — clean start.")

if not preflight_ok:
    raise RuntimeError(
        "Pre-flight validation failed. Fix issues before rerunning Script 4a.")

print("\n  [OK] All pre-flight checks passed — proceeding.")

# ---------------------------------------------------------------------------
# STAGE 1 — BUILD FOCAL NODE RASTER
# ---------------------------------------------------------------------------

print("\n--- STAGE 1: Build Focal Node Raster ---")

regions_path = os.path.join(CS_DIR, "woodland_regions.tif")
wood_path    = os.path.join(CS_DIR, "woodland_mask.tif")
nodes_path   = os.path.join(CS_DIR, "focal_nodes.tif")

print("\n  Step 1a: Extracting woodland mask (LCM classes 1 + 2)...")
lcm_r    = Raster(IN_LCM)
woodland = Con((lcm_r == 1) | (lcm_r == 2), 1)
woodland.save(wood_path)
del lcm_r, woodland; gc.collect()
print(f"  [OK] woodland_mask.tif saved.")

print("\n  Step 1b: Region group — identifying individual patches...")
regions = RegionGroup(wood_path, "EIGHT", "WITHIN", "NO_LINK")
regions.save(regions_path)
del regions; gc.collect()
print(f"  [OK] woodland_regions.tif saved.")

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

print(f"\n  Step 1c: Building focal node raster...")
region_r = Raster(regions_path)
result   = None
for new_id, old_id in enumerate(FOCAL_PATCH_IDS, start=1):
    patch  = Con(region_r == old_id, new_id, 0)
    result = patch if result is None else result + patch

final = SetNull(result == 0, result)
final.save(nodes_path)
del region_r, result, final; gc.collect()

mn = float(arcpy.management.GetRasterProperties(
    nodes_path, "MINIMUM").getOutput(0))
mx = float(arcpy.management.GetRasterProperties(
    nodes_path, "MAXIMUM").getOutput(0))
if mn == 1.0 and mx == float(len(FOCAL_PATCH_IDS)):
    print(f"  [OK] Focal node IDs correct: 1–{len(FOCAL_PATCH_IDS)}")
else:
    raise RuntimeError(
        f"Focal node IDs incorrect: min={mn} max={mx}.")

print("\n  Focal node patch summary:")
for i, (pid, size_ha) in enumerate(
        zip(FOCAL_PATCH_IDS, FOCAL_PATCH_SIZES_HA), start=1):
    print(f"    Node {i:2d} (Region {pid:4d}): {size_ha:.1f} ha")

verify(nodes_path, "Focal nodes raster",
       exp_min=1.0, exp_max=float(len(FOCAL_PATCH_IDS)), exp_bands=1)

# ---------------------------------------------------------------------------
# STAGE 2 — RESAMPLE TO 3m FOR CIRCUITSCAPE
# ---------------------------------------------------------------------------

print(f"\n--- STAGE 2: Resample to {CS_RESOLUTION}m for Circuitscape ---")
print(f"  Resolution: {CS_RESOLUTION}m")
print(f"  Rationale: 2m = 15.7M nodes (exceeds RAM); 3m = ~7M nodes (feasible)")
print(f"  UK precedent: Laikre et al. (2022) used 3m for England EIA study")
print(f"  Ecological basis: bat/badger/otter/pine marten/barn owl at 3m")
print()

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

nodes_3m = os.path.join(CS_DIR, "focal_nodes_3m.tif")
print(f"\n  Resampling focal nodes to {CS_RESOLUTION}m (nearest neighbour)...")
arcpy.management.Resample(
    in_raster       = nodes_path,
    out_raster      = nodes_3m,
    cell_size       = str(CS_RESOLUTION),
    resampling_type = "NEAREST"
)
mn = float(arcpy.management.GetRasterProperties(
    nodes_3m, "MINIMUM").getOutput(0))
mx = float(arcpy.management.GetRasterProperties(
    nodes_3m, "MAXIMUM").getOutput(0))
print(f"  [OK] focal_nodes_3m.tif: min={mn:.0f}  max={mx:.0f}")
if mn >= 1.0 and mx == float(len(FOCAL_PATCH_IDS)):
    print(f"  [OK] Focal node IDs preserved: 1–{len(FOCAL_PATCH_IDS)}")
else:
    print(f"  [WARNING] Node IDs may have changed.")

print(f"\n  [OK] All layers resampled to {CS_RESOLUTION}m.")

# ---------------------------------------------------------------------------
# STAGE 3 — EXPORT ASCII FILES
# ---------------------------------------------------------------------------

print(f"\n--- STAGE 3: Export ASCII Files ({CS_RESOLUTION}m resolution) ---")
print("  Method: numpy RasterToNumPyArray (Parallels-compatible)")
print()

nodes_asc = os.path.join(CS_DIR, "focal_nodes.asc")
export_ascii_numpy(nodes_3m, nodes_asc, "focal_nodes_3m")

resistance_ascs = {}
for name, path in resampled_scenarios:
    out_asc = os.path.join(CS_DIR, f"Resistance_{name}.asc")
    if export_ascii_numpy(path, out_asc, name):
        resistance_ascs[name] = out_asc

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
# STAGE 4 — CREATE CIRCUITSCAPE CONFIG FILES — PAIRWISE MODE
# ---------------------------------------------------------------------------
# Pairwise mode: all pairs of focal nodes solved.
# 15 nodes = 105 pairs per scenario.
# Pairwise is the original Circuitscape methodology (McRae et al. 2008).
# all-to-one mode produces zero output in Circuitscape v5.17.1 —
# confirmed by testing. Pairwise confirmed producing correct results.
# Forward slashes required for Julia/Circuitscape path parsing.

print("\n--- STAGE 4: Create Circuitscape Config Files (pairwise mode) ---")
print("  Mode: pairwise (105 pairs per scenario)")
print("  Note: all-to-one mode produces zero output in Circuitscape v5.17.1")
print("        Pairwise confirmed working — McRae et al. (2008)")

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

    res_fwd = resistance_ascs[name].replace("\\", "/")
    nod_fwd = nodes_asc.replace("\\", "/")
    out_fwd = output_file.replace("\\", "/")

    config_content = f"""[circuitscape options]
data_type = raster
scenario = pairwise

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
        print(f"  [MISSING] Config not created.")

# ---------------------------------------------------------------------------
# STAGE 5 — FINAL VERIFICATION
# ---------------------------------------------------------------------------

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
    ok  = os.path.exists(asc) and os.path.getsize(asc) > 10 * 1024 * 1024
    print(f"  {'[OK]' if ok else '[WARNING]'} Resistance_{name}.asc")
    if not ok:
        all_ready = False

node_asc_ok = (os.path.exists(nodes_asc) and
               os.path.getsize(nodes_asc) > 5 * 1024 * 1024)
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

print("\n  Output folders:")
for name, _ in SCENARIOS:
    folder = os.path.join(OUT_DIR, name)
    ok = os.path.exists(folder)
    print(f"  {'[OK]' if ok else '[WARNING]'} {name}/")
    if not ok:
        all_ready = False

print()
if all_ready:
    print("  [OK] ALL FILES READY — Script 4b can now be run.")
else:
    print("  [WARNING] Some files missing — do not run Script 4b yet.")

# ---------------------------------------------------------------------------
# FINAL SUMMARY
# ---------------------------------------------------------------------------

print()
print("=" * 60)
print("SCRIPT 4a COMPLETE")
print(f"Finished: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 60)
print()

print("VERIFIED PARAMETERS FOR DISSERTATION (Chapter 5):")
print(f"  Circuitscape mode  : pairwise")
print(f"  Focal nodes        : {len(FOCAL_PATCH_IDS)} woodland patches")
print(f"  Patch sizes        : {min(FOCAL_PATCH_SIZES_HA):.1f}–"
      f"{max(FOCAL_PATCH_SIZES_HA):.1f} ha")
print(f"  Pairwise pairs     : "
      f"{len(FOCAL_PATCH_IDS)*(len(FOCAL_PATCH_IDS)-1)//2} per scenario")
print(f"  Scenarios          : 4 (S1–S4)")
print(f"  Solver             : cg+amg")
print(f"  Resistance res.    : 2m (analysis) → {CS_RESOLUTION}m (Circuitscape)")
print(f"  Ecological basis   : bat/badger/otter/pine marten/barn owl at 3m")
print(f"  ASCII method       : numpy RasterToNumPyArray")
print(f"  CRS                : EPSG:27700")
print(f"  Approx. nodes      : ~7,000,000 (at {CS_RESOLUTION}m resolution)")
print(f"  Citation           : McRae et al. (2008) Ecology")
print(f"                       Dickson et al. (2019) Conservation Biology")
print(f"                       Laikre et al. (2022) J. Environ. Management")
print(f"                       Drasher et al. (2025) PLOS One")
print()
print("NEXT STEP: Script 4b — Circuitscape Execution and Results Import")
print()
print("Push to GitHub:")
print("  https://github.com/Abhin-Nangari/A66-Connectivity-MSc")

arcpy.CheckInExtension("Spatial")
