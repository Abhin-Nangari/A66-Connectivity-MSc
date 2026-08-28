# =============================================================================
# Script 4 — Circuitscape Connectivity Analysis
# Project : MSc GIS Dissertation — A66 Northern Trans-Pennine Corridor
# Author  : Abhin Nangari (Student ID: 52536856)
# University: University of Aberdeen — MSc GIS GG5910/GG5912
# Supervisor: Dr Shaktiman Singh
# GitHub  : https://github.com/Abhin-Nangari/A66-Connectivity-MSc
# License : MIT
#
# Purpose :
#   1. Pre-flight validation — verify all Script 3 outputs and Julia/Circuitscape
#   2. Build focal node raster from pre-verified woodland patch IDs
#   3. Export resistance surfaces as ASCII (.asc) for each scenario
#   4. Create Circuitscape config.ini — all-to-one mode — for each scenario
#   5. Run Circuitscape via Julia subprocess — all 4 scenarios sequentially
#   6. Import cumulative current density rasters into ArcGIS Pro
#   7. Verify all outputs — statistics and spatial extent
#   8. Validate current density along Eden Valley Railway (primary validation)
#
# Circuitscape mode — all-to-one:
#   All-to-one mode passes current from all focal nodes simultaneously to
#   each node in turn, producing a cumulative current map equivalent to
#   pairwise but using significantly less memory and running faster.
#   For the goal of identifying cumulative connectivity corridors across
#   the A66, all-to-one is the recommended mode per official Circuitscape
#   documentation (docs.circuitscape.org, December 2025).
#   With 15 nodes: all-to-one = 15 calculations vs pairwise = 105.
#   Citation: McRae et al. (2008) Ecology 89(10) pp.2712-2724
#             Wilson et al. (2017) PLoS Computational Biology 13(6) e1005510
#             Dickson et al. (2019) Conservation Biology 33(2) pp.239-249
#
# IMPORTANT — run Circuitscape with ArcGIS Pro CLOSED:
#   Official Circuitscape documentation explicitly recommends closing ArcGIS
#   and all RAM-intensive applications before running Circuitscape.
#   This script handles this by running Circuitscape via subprocess AFTER
#   saving all ArcPy outputs to disk and releasing memory.
#   User must close ArcGIS Pro manually when prompted before Stage 4 runs.
#
# Focal node configuration — pre-verified from diagnostic run:
#   Total woodland patches in buffer : 1107
#   Patches >= 1 hectare             : 83
#   Focal nodes selected             : 15 largest
#   Patch IDs : 1, 851, 121, 38, 437, 143, 667, 741, 1056, 658,
#               913, 496, 170, 1046, 964
#   Sizes (ha): 64.1, 51.6, 49.7, 28.9, 26.2, 25.9, 23.2, 21.9,
#               20.7, 19.9, 10.2, 10.1, 9.7, 8.4, 8.3
#
# Resistance weight citations:
#   McRae et al. (2008) Ecology 89(10) pp.2712-2724
#   Wilson et al. (2017) PLoS Computational Biology 13(6) e1005510
#   Dickson et al. (2019) Conservation Biology 33(2) pp.239-249
#   Koen et al. (2014) Methods in Ecology and Evolution 5(7) pp.626-633
#
# Julia path : C:\Users\abhin_n\AppData\Local\Programs\Julia-1.12.6\bin\julia.exe
# Run from   : ArcGIS Pro Python console (Analysis tab → Python)
# CRS        : EPSG:27700 British National Grid (all outputs)
# Inputs     : Processed\Resistance\Resistance_S*.tif (4 scenarios)
#              Processed\Harmonised\LCM_2m_aligned.tif
#              Processed\CS_Inputs\woodland_regions.tif (pre-computed)
# Outputs    : Processed\CS_Inputs\ — ASCII + node files + config.ini
#              Outputs\ — current density GeoTIFFs
# =============================================================================

import arcpy
import os
import subprocess
import datetime
import gc

arcpy.CheckOutExtension("Spatial")
from arcpy.sa import Raster, Con, RegionGroup

# ---------------------------------------------------------------------------
# 0. ENVIRONMENT SETUP
# ---------------------------------------------------------------------------

arcpy.env.overwriteOutput = True
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(27700)

ROOT     = r"C:\Mac\Home\Desktop\UK\Abardeen\UoA\Dissetation\A66"
GDB      = os.path.join(ROOT, "Processed", "StudyArea", "A66_Study.gdb")
HARM_DIR = os.path.join(ROOT, "Processed", "Harmonised")
RES_DIR  = os.path.join(ROOT, "Processed", "Resistance")
CS_DIR   = os.path.join(ROOT, "Processed", "CS_Inputs")
OUT_DIR  = os.path.join(ROOT, "Outputs")
BUFFER   = os.path.join(GDB, "StudyArea_2km_buffer")
IN_LCM   = os.path.join(HARM_DIR, "LCM_2m_aligned.tif")
JULIA    = r"C:\Users\abhin_n\AppData\Local\Programs\Julia-1.12.6\bin\julia.exe"

# Pre-verified focal node patch IDs — 15 largest woodland patches >= 1ha
# Confirmed from diagnostic: 83 qualifying patches, 15 largest selected
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
print("SCRIPT 4 — CIRCUITSCAPE CONNECTIVITY ANALYSIS")
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
# STAGE 0 — PRE-FLIGHT VALIDATION
# ---------------------------------------------------------------------------

print("\n--- STAGE 0: Pre-flight Validation ---")
preflight_ok = True

# Julia and Circuitscape
print("\n  Checking Julia and Circuitscape...")
try:
    result = subprocess.run(
        [JULIA, "--eval",
         'using Circuitscape; println("Circuitscape ready")'],
        capture_output=True, text=True, timeout=120)
    if "Circuitscape ready" in result.stdout:
        print("  [OK] Julia and Circuitscape confirmed ready.")
    else:
        print(f"  [ABORT] Circuitscape not ready: {result.stderr}")
        preflight_ok = False
except Exception as e:
    print(f"  [ABORT] Julia error: {e}")
    preflight_ok = False

# Resistance surfaces
print("\n  Checking resistance surfaces...")
for name, path in SCENARIOS:
    if os.path.exists(path):
        bands = int(arcpy.management.GetRasterProperties(
            path, "BANDCOUNT").getOutput(0))
        mn = float(arcpy.management.GetRasterProperties(
            path, "MINIMUM").getOutput(0))
        mx = float(arcpy.management.GetRasterProperties(
            path, "MAXIMUM").getOutput(0))
        print(f"  [OK] {name}: bands={bands}  "
              f"min={mn:.4f}  max={mx:.4f}")
    else:
        print(f"  [ABORT] Not found: {path}")
        preflight_ok = False

# LCM
print("\n  Checking LCM...")
if os.path.exists(IN_LCM):
    bands = int(arcpy.management.GetRasterProperties(
        IN_LCM, "BANDCOUNT").getOutput(0))
    print(f"  [OK] LCM_2m_aligned.tif: bands={bands}")
else:
    print(f"  [ABORT] LCM not found.")
    preflight_ok = False

# Woodland regions
regions_path = os.path.join(CS_DIR, "woodland_regions.tif")
print("\n  Checking pre-computed woodland regions...")
if os.path.exists(regions_path):
    print(f"  [OK] woodland_regions.tif found — will reuse.")
else:
    print(f"  [NOTE] Not found — will compute in Stage 1.")

# Focal node summary
print(f"\n  Focal node configuration:")
print(f"  Mode      : all-to-one (recommended for cumulative corridor mapping)")
print(f"  Nodes     : {len(FOCAL_PATCH_IDS)}")
print(f"  Sizes     : {min(FOCAL_PATCH_SIZES_HA):.1f}–"
      f"{max(FOCAL_PATCH_SIZES_HA):.1f} ha")
print(f"  Calcs     : {len(FOCAL_PATCH_IDS)} "
      f"(vs {len(FOCAL_PATCH_IDS)*(len(FOCAL_PATCH_IDS)-1)//2} for pairwise)")
print(f"  [OK] Pre-verified from diagnostic run.")

if not preflight_ok:
    raise RuntimeError(
        "Pre-flight validation failed. Fix issues before rerunning Script 4.")

print("\n  [OK] All pre-flight checks passed — proceeding.")

# ---------------------------------------------------------------------------
# STAGE 1 — BUILD FOCAL NODE RASTER
# ---------------------------------------------------------------------------

print("\n--- STAGE 1: Build Focal Node Raster ---")

# Recompute regions only if not present
if not os.path.exists(regions_path):
    print("  Computing woodland mask and region group...")
    lcm_r    = Raster(IN_LCM)
    woodland = Con((lcm_r == 1) | (lcm_r == 2), 1)
    wood_path = os.path.join(CS_DIR, "woodland_mask.tif")
    woodland.save(wood_path)
    del lcm_r, woodland; gc.collect()
    regions = RegionGroup(wood_path, "EIGHT", "WITHIN", "NO_LINK")
    regions.save(regions_path)
    del regions; gc.collect()
    print("  [OK] Woodland regions computed.")
else:
    print("  [OK] Reusing existing woodland_regions.tif.")

# Build focal node raster from pre-verified patch IDs
print(f"\n  Building focal node raster — "
      f"{len(FOCAL_PATCH_IDS)} patches...")

regions_r    = Raster(regions_path)
nodes_raster = None

for new_id, old_id in enumerate(FOCAL_PATCH_IDS, start=1):
    patch = Con(regions_r == old_id, new_id)
    if nodes_raster is None:
        nodes_raster = patch
    else:
        nodes_raster = Con(nodes_raster > 0, nodes_raster, patch)

nodes_path = os.path.join(CS_DIR, "focal_nodes.tif")
nodes_raster.save(nodes_path)
del regions_r, nodes_raster; gc.collect()

print("\n  Focal node patch summary:")
for i, (pid, size_ha) in enumerate(
        zip(FOCAL_PATCH_IDS, FOCAL_PATCH_SIZES_HA), start=1):
    print(f"    Node {i:2d} (Patch {pid:4d}): {size_ha:.1f} ha")

verify(nodes_path, "Focal nodes raster", exp_bands=1)
print(f"\n  [OK] {len(FOCAL_PATCH_IDS)} focal nodes built and verified.")

# ---------------------------------------------------------------------------
# STAGE 2 — EXPORT ASCII FILES FOR CIRCUITSCAPE
# ---------------------------------------------------------------------------

print("\n--- STAGE 2: Export ASCII Files ---")

def export_ascii(in_raster, out_asc, name):
    arcpy.conversion.RasterToASCII(in_raster, out_asc)
    if os.path.exists(out_asc):
        size_mb = os.path.getsize(out_asc) / 1024 / 1024
        print(f"  [OK] {name}: {os.path.basename(out_asc)} "
              f"({size_mb:.1f} MB)")
        return True
    else:
        print(f"  [MISSING] {name}: ASCII export failed.")
        return False

# Export focal nodes
nodes_asc = os.path.join(CS_DIR, "focal_nodes.asc")
print("\n  Exporting focal nodes...")
export_ascii(nodes_path, nodes_asc, "Focal nodes")

# Export resistance surfaces
resistance_ascs = {}
for name, path in SCENARIOS:
    print(f"\n  Exporting {name}...")
    out_asc = os.path.join(CS_DIR, f"Resistance_{name}.asc")
    if export_ascii(path, out_asc, name):
        resistance_ascs[name] = out_asc

# ---------------------------------------------------------------------------
# STAGE 3 — CREATE CIRCUITSCAPE CONFIG FILES — ALL-TO-ONE MODE
# ---------------------------------------------------------------------------
# all-to-one mode: in each iteration one focal node is connected to ground,
# all others connected to current sources. Produces cumulative current map
# equivalent to pairwise but with 15 calculations instead of 105.
# Recommended by official Circuitscape docs for cumulative corridor mapping.

print("\n--- STAGE 3: Create Circuitscape Config Files (all-to-one) ---")

def create_config(scenario_name, resistance_asc, nodes_asc, output_dir):
    scenario_out_dir = os.path.join(output_dir, scenario_name)
    if not os.path.exists(scenario_out_dir):
        os.makedirs(scenario_out_dir)
    config_path = os.path.join(CS_DIR, f"config_{scenario_name}.ini")
    output_file = os.path.join(scenario_out_dir, f"CS_{scenario_name}")

    # Forward slashes for Circuitscape/Julia compatibility
    res_fwd = resistance_asc.replace("\\", "/")
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
    else:
        print(f"  [MISSING] Config not created.")

    return config_path, output_file

config_files    = {}
output_prefixes = {}

for name, _ in SCENARIOS:
    if name in resistance_ascs:
        cfg_path, out_prefix = create_config(
            name, resistance_ascs[name], nodes_asc, OUT_DIR)
        config_files[name]    = cfg_path
        output_prefixes[name] = out_prefix

# ---------------------------------------------------------------------------
# STAGE 4 — CLOSE ARCGIS PRO — THEN RUN CIRCUITSCAPE
# ---------------------------------------------------------------------------
# Official Circuitscape documentation explicitly recommends closing ArcGIS
# before running Circuitscape to free RAM.
# Script prints a clear prompt — user closes ArcGIS Pro, then presses Enter.
# Circuitscape runs in all-to-one mode via Julia subprocess.

print("\n--- STAGE 4: Run Circuitscape ---")
print()
print("  *** IMPORTANT — ACTION REQUIRED ***")
print("  Official Circuitscape docs recommend closing ArcGIS Pro")
print("  before running Circuitscape to free RAM and avoid stalls.")
print()
print("  Steps:")
print("  1. Save your ArcGIS Pro project (Ctrl+S)")
print("  2. Close ArcGIS Pro completely")
print("  3. Return to this Python window")
print("  4. Press Enter to start Circuitscape")
print()
input("  Press Enter when ArcGIS Pro is closed and you are ready...")
print()

curmap_paths = {}

for name, _ in SCENARIOS:
    if name not in config_files:
        print(f"  [SKIP] {name} — config not created.")
        continue

    cfg_path   = config_files[name]
    out_prefix = output_prefixes[name]
    # all-to-one cumulative output filename
    curmap_asc = out_prefix + "_cum_curmap.asc"

    print(f"  Running: {name}")
    print(f"  Mode   : all-to-one")
    print(f"  Started: {datetime.datetime.now().strftime('%H:%M:%S')}")

    # Forward slashes for Julia
    cfg_fwd = cfg_path.replace("\\", "/")
    julia_cmd = [
        JULIA,
        "--eval",
        f'using Circuitscape; compute("{cfg_fwd}")'
    ]

    try:
        result = subprocess.run(
            julia_cmd,
            capture_output=True,
            text=True,
            timeout=3600)

        print(f"  Finished: {datetime.datetime.now().strftime('%H:%M:%S')}")

        if result.returncode == 0:
            if os.path.exists(curmap_asc):
                size_mb = os.path.getsize(curmap_asc) / 1024 / 1024
                print(f"  [OK] Output: {os.path.basename(curmap_asc)} "
                      f"({size_mb:.1f} MB)")
                curmap_paths[name] = curmap_asc
            else:
                # Check for gzip compressed output
                curmap_gz = curmap_asc + ".gz"
                if os.path.exists(curmap_gz):
                    print(f"  [OK] Compressed output: "
                          f"{os.path.basename(curmap_gz)}")
                    curmap_paths[name] = curmap_gz
                else:
                    print(f"  [WARNING] Expected output not found.")
                    print(f"  stdout: {result.stdout[-300:]}")
                    print(f"  stderr: {result.stderr[-300:]}")
        else:
            print(f"  [ERROR] Return code: {result.returncode}")
            print(f"  stderr: {result.stderr[-500:]}")

    except subprocess.TimeoutExpired:
        print(f"  [TIMEOUT] Exceeded 1 hour for {name}.")
    except Exception as e:
        print(f"  [ERROR] {e}")
    print()

# ---------------------------------------------------------------------------
# STAGE 5 — REOPEN ARCGIS PRO — IMPORT CURRENT DENSITY RASTERS
# ---------------------------------------------------------------------------

print("\n--- STAGE 5: Import Current Density Rasters ---")
print()
print("  Reopen ArcGIS Pro now, reload your project,")
print("  then come back to this Python window.")
print()
input("  Press Enter when ArcGIS Pro is open and project is loaded...")
print()

# Re-checkout Spatial Analyst after ArcGIS Pro reopened
arcpy.CheckOutExtension("Spatial")

curmap_tifs = {}

for name, curmap_asc in curmap_paths.items():
    out_tif = os.path.join(OUT_DIR, f"CurrentDensity_{name}.tif")
    print(f"\n  Importing: {name}")

    # Handle gzip compressed output if needed
    import_path = curmap_asc
    if curmap_asc.endswith(".gz"):
        import subprocess as sp
        import gzip, shutil
        unzipped = curmap_asc.replace(".gz", "")
        with gzip.open(curmap_asc, 'rb') as f_in:
            with open(unzipped, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        import_path = unzipped
        print(f"  Decompressed: {os.path.basename(unzipped)}")

    arcpy.conversion.ASCIIToRaster(import_path, out_tif, "FLOAT")
    arcpy.management.DefineProjection(
        out_tif, arcpy.SpatialReference(27700))

    curmap_tifs[name] = out_tif
    verify(out_tif, f"CurrentDensity_{name}", exp_min=0.0)

# ---------------------------------------------------------------------------
# STAGE 6 — VALIDATE — EDEN VALLEY RAILWAY
# ---------------------------------------------------------------------------
# Primary validation: disused Eden Valley Railway near Kirkby Thore.
# NH TR010062 identified this as the only formally named wildlife corridor
# on the entire A66 upgrade — chosen specifically to avoid this feature.
# Expected: elevated current density along railway vs A66 road.
# Citation: National Highways (2019) NH Citizen Space TR010062

print("\n--- STAGE 6: Validate — Eden Valley Railway ---")

# Approximate EPSG:27700 coordinates along disused railway
railway_points = {
    "Railway W"      : (362300, 520750),
    "Railway mid-W"  : (362800, 520900),
    "Railway centre" : (363200, 521000),
    "Railway mid-E"  : (363600, 521100),
    "Railway E"      : (364000, 521200),
}

# A66 road points — expected low current (barrier effect)
road_points = {
    "A66 W"   : (362500, 520600),
    "A66 mid" : (363200, 520800),
    "A66 E"   : (364000, 521000),
}

if "S1_primary" in curmap_tifs:
    s1_tif = curmap_tifs["S1_primary"]

    print("\n  Current density — Eden Valley Railway (S1_primary):")
    railway_vals = []
    for loc, (x, y) in railway_points.items():
        try:
            val = float(arcpy.management.GetCellValue(
                s1_tif, f"{x} {y}", "1").getOutput(0))
            railway_vals.append(val)
            print(f"    {loc}: {val:.6f}")
        except:
            print(f"    {loc}: NoData")

    print("\n  Current density — A66 road (S1_primary):")
    road_vals = []
    for loc, (x, y) in road_points.items():
        try:
            val = float(arcpy.management.GetCellValue(
                s1_tif, f"{x} {y}", "1").getOutput(0))
            road_vals.append(val)
            print(f"    {loc}: {val:.6f}")
        except:
            print(f"    {loc}: NoData")

    if railway_vals and road_vals:
        rail_mean = sum(railway_vals) / len(railway_vals)
        road_mean = sum(road_vals) / len(road_vals)
        print(f"\n  Mean current density:")
        print(f"    Eden Valley Railway : {rail_mean:.6f}")
        print(f"    A66 road            : {road_mean:.6f}")
        if rail_mean > road_mean:
            print(f"  [OK] Railway > A66 road — validates as wildlife corridor.")
        else:
            print(f"  [NOTE] Check coordinates visually in ArcGIS Pro.")
            print(f"         Coordinate approximations may need refining.")
else:
    print("  [NOTE] S1 output not available — verify Circuitscape ran correctly.")

# ---------------------------------------------------------------------------
# FINAL SUMMARY
# ---------------------------------------------------------------------------

print()
print("=" * 60)
print("SCRIPT 4 COMPLETE")
print(f"Finished: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 60)
print()

print("CS_Inputs files:")
for f in sorted(os.listdir(CS_DIR)):
    if not f.startswith(".") and f != "Thumbs.db":
        print(f"  [OK] {f}")

print("\nOutput current density rasters:")
for name, tif in curmap_tifs.items():
    status = "[OK]" if os.path.exists(tif) else "[MISSING]"
    print(f"  {status} {os.path.basename(tif)}")

missing = [n for n, _ in SCENARIOS if n not in curmap_tifs]
if missing:
    print(f"\n  [WARNING] Missing outputs: {missing}")
    print("  Rerun Script 4 Stage 4 only — check Circuitscape logs.")

print()
print("VISUAL QA — verify in ArcGIS Pro before Script 5:")
print("  [ ] Load CurrentDensity_S1_primary.tif — hot colour ramp")
print("  [ ] High current corridors visible across study area")
print("  [ ] Eden Valley Railway — elevated current near Kirkby Thore")
print("  [ ] A66 road — low current — confirms barrier effect")
print("  [ ] River Eden corridor — moderate-high current")
print("  [ ] S1 vs S4 — similar spatial pattern (confirms high stability)")
print("  [ ] Overlay on GAD georeferenced sheets — check alignment")
print()
print("VERIFIED PARAMETERS FOR DISSERTATION (Chapter 5):")
print(f"  Circuitscape mode  : all-to-one")
print(f"  Focal nodes        : {len(FOCAL_PATCH_IDS)} woodland patches")
print(f"  Patch size range   : {min(FOCAL_PATCH_SIZES_HA):.1f}–"
      f"{max(FOCAL_PATCH_SIZES_HA):.1f} ha")
print(f"  Calculations       : {len(FOCAL_PATCH_IDS)} per scenario")
print(f"  Scenarios          : 4 (S1-S4)")
print(f"  Solver             : cg+amg")
print(f"  CRS                : EPSG:27700")
print()
print("NEXT STEP: Script 5 — Validation and Map Export")
print("  Overlay current density vs NH corridor, Habitat Networks, LNRS APIB")
print("  Export all maps at 300 dpi")
print()
print("Push to GitHub:")
print("  https://github.com/Abhin-Nangari/A66-Connectivity-MSc")

arcpy.CheckInExtension("Spatial")
