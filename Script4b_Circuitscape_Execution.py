# =============================================================================
# Script 4b — Circuitscape Execution and Results Import
# Project : MSc GIS Dissertation — A66 Northern Trans-Pennine Corridor
# Author  : Abhin Nangari (Student ID: 52536856)
# University: University of Aberdeen — MSc GIS GG5910/GG5912
# Supervisor: Dr Shaktiman Singh
# GitHub  : https://github.com/Abhin-Nangari/A66-Connectivity-MSc
# License : MIT
#
# Purpose :
#   Runs Circuitscape and imports results — run AFTER Script 4a.
#   1. Pre-flight validation — verify all Script 4a outputs exist
#   2. Run Circuitscape via Julia subprocess — all 4 scenarios sequentially
#   3. Verify Circuitscape outputs — confirm all 4 cumulative current maps
#   4. Import cumulative current density rasters into ArcGIS Pro
#   5. Verify imported GeoTIFFs — band count, extent, statistics
#   6. Validate current density along Eden Valley Railway (primary validation)
#   7. Compare Railway vs A66 road current density — quantitative validation
#
# Prerequisites — Script 4a must be complete and verified before running:
#   Processed\CS_Inputs\focal_nodes.asc
#   Processed\CS_Inputs\Resistance_S*.asc (4 files)
#   Processed\CS_Inputs\config_S*.ini (4 files)
#
# Circuitscape mode — all-to-one:
#   Each focal node connected to ground in turn, all others act as current
#   sources. Produces cumulative current map identifying all movement
#   corridors and pinch points. 15 calculations per scenario.
#   Citation: McRae et al. (2008) Ecology 89(10) pp.2712-2724
#             Wilson et al. (2017) PLoS Computational Biology 13(6) e1005510
#             Dickson et al. (2019) Conservation Biology 33(2) pp.239-249
#
# Eden Valley Railway validation:
#   Primary validation layer per NH TR010062 methodology.
#   National Highways selected the Temple Sowerby–Appleby bypass alignment
#   specifically to avoid the disused Eden Valley Railway near Kirkby Thore
#   — the only formally named wildlife corridor on the entire A66 upgrade.
#   Expected: elevated current density along railway vs A66 road surface.
#   Citation: National Highways (2019) NH Citizen Space TR010062
#             SI 2024/360 Schedule 7 Part 3
#
# IMPORTANT — Official Circuitscape documentation recommends closing
# RAM-intensive applications before running. If ArcGIS Pro is open,
# close it before running Circuitscape (Stage 2), then reopen for
# Stage 3 onwards. Script prompts user at the correct moment.
#
# Run from : ArcGIS Pro Python console (Analysis tab → Python)
#            OR Windows Command Prompt using ArcGIS Python executable
# Julia     : C:\Users\abhin_n\AppData\Local\Programs\Julia-1.12.6\bin\julia.exe
# CRS       : EPSG:27700 British National Grid (all outputs)
# Inputs    : Processed\CS_Inputs\ (from Script 4a)
# Outputs   : Outputs\CurrentDensity_S*.tif (4 GeoTIFFs)
# =============================================================================

import arcpy
import os
import subprocess
import datetime
import gc
import shutil
import gzip

arcpy.CheckOutExtension("Spatial")

# ---------------------------------------------------------------------------
# 0. ENVIRONMENT SETUP
# ---------------------------------------------------------------------------

arcpy.env.overwriteOutput = True
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(27700)

ROOT    = r"C:\Mac\Home\Desktop\UK\Abardeen\UoA\Dissetation\A66"
RES_DIR = os.path.join(ROOT, "Processed", "Resistance")
CS_DIR  = os.path.join(ROOT, "Processed", "CS_Inputs")
OUT_DIR = os.path.join(ROOT, "Outputs")
JULIA   = r"C:\Users\abhin_n\AppData\Local\Programs\Julia-1.12.6\bin\julia.exe"

# Scenario names — must match Script 4a config file names
SCENARIO_NAMES = [
    "S1_primary",
    "S2_LCM_plus20pp",
    "S3_LCM_minus20pp",
    "S4_equal",
]

# Master grid from S1 resistance surface — for extent verification
REF        = os.path.join(RES_DIR, "Resistance_S1_primary.tif")
arcpy.env.snapRaster = REF
arcpy.env.cellSize   = 2
ref_desc   = arcpy.Describe(REF)
ref_ext    = ref_desc.extent
MASTER_EXT = (ref_ext.XMin, ref_ext.YMin, ref_ext.XMax, ref_ext.YMax)

if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)

print("=" * 60)
print("SCRIPT 4b — CIRCUITSCAPE EXECUTION AND RESULTS IMPORT")
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
    print(f"    Min={r_min:.6f}  Max={r_max:.6f}  "
          f"Mean={r_mean:.6f}  STD={r_std:.6f}")
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
# Verifies all Script 4a outputs exist and are correctly sized.
# Aborts immediately if any required file is missing.

print("\n--- STAGE 0: Pre-flight Validation (Script 4a Outputs) ---")
preflight_ok = True

print("\n  Checking focal nodes ASCII...")
nodes_asc = os.path.join(CS_DIR, "focal_nodes.asc")
if os.path.exists(nodes_asc) and os.path.getsize(nodes_asc) > 50*1024*1024:
    size_mb = os.path.getsize(nodes_asc) / 1024 / 1024
    print(f"  [OK] focal_nodes.asc ({size_mb:.1f} MB)")
else:
    print(f"  [ABORT] focal_nodes.asc missing or too small. Run Script 4a first.")
    preflight_ok = False

print("\n  Checking resistance ASCII files...")
for name in SCENARIO_NAMES:
    asc = os.path.join(CS_DIR, f"Resistance_{name}.asc")
    if os.path.exists(asc) and os.path.getsize(asc) > 100*1024*1024:
        size_mb = os.path.getsize(asc) / 1024 / 1024
        print(f"  [OK] Resistance_{name}.asc ({size_mb:.1f} MB)")
    else:
        print(f"  [ABORT] Resistance_{name}.asc missing. Run Script 4a first.")
        preflight_ok = False

print("\n  Checking config files...")
config_files    = {}
output_prefixes = {}
for name in SCENARIO_NAMES:
    cfg = os.path.join(CS_DIR, f"config_{name}.ini")
    if os.path.exists(cfg) and os.path.getsize(cfg) > 0:
        print(f"  [OK] config_{name}.ini")
        config_files[name] = cfg
        # Read output path from config
        with open(cfg) as f:
            for line in f:
                if line.strip().startswith("output_file"):
                    out_path = line.split("=")[1].strip().replace("/", "\\")
                    output_prefixes[name] = out_path
                    break
    else:
        print(f"  [ABORT] config_{name}.ini missing. Run Script 4a first.")
        preflight_ok = False

print("\n  Checking Julia...")
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

print("\n  Checking disk space...")
_, _, free = shutil.disk_usage(r"C:\Mac\Home")
free_gb = free / 1024 / 1024 / 1024
print(f"  Free: {free_gb:.1f} GB "
      f"{'[OK]' if free_gb >= 5 else '[WARNING] Less than 5GB'}")
if free_gb < 5:
    preflight_ok = False

if not preflight_ok:
    raise RuntimeError(
        "Pre-flight validation failed. "
        "Run Script 4a and fix issues before rerunning Script 4b.")

print("\n  [OK] All pre-flight checks passed — proceeding.")

# ---------------------------------------------------------------------------
# STAGE 1 — PROMPT: CLOSE ARCGIS PRO BEFORE CIRCUITSCAPE
# ---------------------------------------------------------------------------
# Official Circuitscape documentation recommends closing ArcGIS Pro
# before running Circuitscape to free RAM and prevent stalls.
# This script runs Circuitscape via subprocess — ArcGIS Pro does not
# need to be open for this stage.

print("\n--- STAGE 1: Close ArcGIS Pro Before Circuitscape ---")
print()
print("  *** ACTION REQUIRED ***")
print("  Close ArcGIS Pro before running Circuitscape to free RAM.")
print("  Circuitscape runs independently — ArcGIS Pro not needed.")
print()
print("  1. Press Ctrl+S in ArcGIS Pro to save project")
print("  2. Close ArcGIS Pro completely")
print("  3. Return here and press Enter")
print()
input("  Press Enter when ArcGIS Pro is closed...")
print()
print(f"  Circuitscape starting: "
      f"{datetime.datetime.now().strftime('%H:%M:%S')}")

# ---------------------------------------------------------------------------
# STAGE 2 — RUN CIRCUITSCAPE VIA JULIA — ALL 4 SCENARIOS
# ---------------------------------------------------------------------------
# All-to-one mode. 15 calculations per scenario.
# Each scenario run sequentially — results saved to Outputs\{scenario_name}\
# Timeout: 1 hour per scenario.
# Both standard (.asc) and gzip-compressed (.asc.gz) outputs handled.

print("\n--- STAGE 2: Run Circuitscape ---")
print(f"  Mode     : all-to-one")
print(f"  Nodes    : 15")
print(f"  Calcs    : 15 per scenario (vs 105 pairwise)")
print(f"  Scenarios: 4")
print(f"  Solver   : cg+amg")
print()

curmap_paths = {}

for name in SCENARIO_NAMES:
    if name not in config_files:
        print(f"  [SKIP] {name} — config not found.")
        continue

    cfg_path   = config_files[name]
    out_prefix = output_prefixes.get(name, "")
    curmap_asc = out_prefix + "_cum_curmap.asc"
    cfg_fwd    = cfg_path.replace("\\", "/")

    print(f"  Running  : {name}")
    print(f"  Started  : {datetime.datetime.now().strftime('%H:%M:%S')}")

    julia_cmd = [
        JULIA, "--eval",
        f'using Circuitscape; compute("{cfg_fwd}")'
    ]

    try:
        result = subprocess.run(
            julia_cmd,
            capture_output=True,
            text=True,
            timeout=3600)

        print(f"  Finished : {datetime.datetime.now().strftime('%H:%M:%S')}")

        if result.returncode == 0:
            if os.path.exists(curmap_asc):
                size_mb = os.path.getsize(curmap_asc) / 1024 / 1024
                print(f"  [OK] {os.path.basename(curmap_asc)} "
                      f"({size_mb:.1f} MB)")
                curmap_paths[name] = curmap_asc
            elif os.path.exists(curmap_asc + ".gz"):
                size_mb = os.path.getsize(
                    curmap_asc + ".gz") / 1024 / 1024
                print(f"  [OK] Compressed output ({size_mb:.1f} MB)")
                curmap_paths[name] = curmap_asc + ".gz"
            else:
                print(f"  [WARNING] Output not found: {curmap_asc}")
                print(f"  Stdout: {result.stdout[-300:]}")
                print(f"  Stderr: {result.stderr[-300:]}")
        else:
            print(f"  [ERROR] Return code: {result.returncode}")
            print(f"  Stderr: {result.stderr[-500:]}")

    except subprocess.TimeoutExpired:
        print(f"  [TIMEOUT] Exceeded 1 hour for {name}.")
        print(f"  Consider reducing study area extent or node count.")
    except Exception as e:
        print(f"  [ERROR] {e}")
    print()

print(f"  Circuitscape complete: "
      f"{datetime.datetime.now().strftime('%H:%M:%S')}")
print(f"  Successful outputs: {len(curmap_paths)} of {len(SCENARIO_NAMES)}")

if len(curmap_paths) == 0:
    raise RuntimeError(
        "No Circuitscape outputs produced. "
        "Check Julia installation and config file paths.")

# ---------------------------------------------------------------------------
# STAGE 3 — VERIFY CIRCUITSCAPE OUTPUTS
# ---------------------------------------------------------------------------

print("\n--- STAGE 3: Verify Circuitscape Outputs ---")
for name, path in curmap_paths.items():
    exists = os.path.exists(path)
    size_mb = os.path.getsize(path) / 1024 / 1024 if exists else 0
    print(f"  {'[OK]' if exists else '[MISSING]'} {name}: "
          f"{os.path.basename(path)} ({size_mb:.1f} MB)")

# ---------------------------------------------------------------------------
# STAGE 4 — PROMPT: REOPEN ARCGIS PRO
# ---------------------------------------------------------------------------

print("\n--- STAGE 4: Reopen ArcGIS Pro ---")
print()
print("  Reopen ArcGIS Pro and reload your project.")
print("  Then return here and press Enter.")
print()
input("  Press Enter when ArcGIS Pro is open and project loaded...")
print()

arcpy.CheckOutExtension("Spatial")

# ---------------------------------------------------------------------------
# STAGE 5 — IMPORT CURRENT DENSITY RASTERS INTO ARCGIS PRO
# ---------------------------------------------------------------------------
# Converts Circuitscape ASCII outputs to GeoTIFF.
# Applies EPSG:27700 coordinate system.
# Handles gzip-compressed outputs automatically.
# Verifies each imported raster — band count, extent, statistics.

print("\n--- STAGE 5: Import Current Density Rasters ---")

curmap_tifs = {}

for name, curmap_src in curmap_paths.items():
    out_tif = os.path.join(OUT_DIR, f"CurrentDensity_{name}.tif")
    print(f"\n  Importing: {name}")

    # Decompress if gzip compressed
    import_path = curmap_src
    if curmap_src.endswith(".gz"):
        unzipped = curmap_src.replace(".gz", "")
        with gzip.open(curmap_src, 'rb') as f_in:
            with open(unzipped, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        import_path = unzipped
        size_mb = os.path.getsize(unzipped) / 1024 / 1024
        print(f"  Decompressed: {os.path.basename(unzipped)} "
              f"({size_mb:.1f} MB)")

    arcpy.conversion.ASCIIToRaster(import_path, out_tif, "FLOAT")
    arcpy.management.DefineProjection(
        out_tif, arcpy.SpatialReference(27700))

    curmap_tifs[name] = out_tif
    verify(out_tif, f"CurrentDensity_{name}", exp_min=0.0, exp_bands=1)

# ---------------------------------------------------------------------------
# STAGE 6 — VALIDATE — EDEN VALLEY RAILWAY
# ---------------------------------------------------------------------------
# Primary validation layer: disused Eden Valley Railway near Kirkby Thore.
# National Highways TR010062 (2019) identified this as the only formally
# named wildlife corridor on the entire A66 upgrade. The DCO alignment
# (SI 2024/360 Schedule 7 Part 3) was specifically routed to avoid this
# feature — making it the most important validation location in the study.
#
# Validation method:
#   Sample current density at 5 points along the railway centreline and
#   3 points along the A66 road surface. Compare mean values.
#   Expected: railway mean > A66 road mean — confirms circuit theory
#   correctly identifies the railway as a preferred movement corridor.
#
# Citation: National Highways (2019) NH Citizen Space TR010062
#           SI 2024/360 Schedule 7 Part 3
#           McRae et al. (2008) Ecology 89(10) pp.2712-2724

print("\n--- STAGE 6: Validate — Eden Valley Railway ---")

# Approximate EPSG:27700 coordinates along disused railway centreline
railway_points = {
    "Railway W end"  : (362300, 520750),
    "Railway mid-W"  : (362800, 520900),
    "Railway centre" : (363200, 521000),
    "Railway mid-E"  : (363600, 521100),
    "Railway E end"  : (364000, 521200),
}

# A66 road points — expected low current (barrier effect)
road_points = {
    "A66 W"   : (362500, 520600),
    "A66 mid" : (363200, 520800),
    "A66 E"   : (364000, 521000),
}

validation_results = {}

for scenario_name, tif_path in curmap_tifs.items():
    print(f"\n  Scenario: {scenario_name}")

    railway_vals = []
    print("  Current density — Eden Valley Railway:")
    for loc, (x, y) in railway_points.items():
        try:
            val = float(arcpy.management.GetCellValue(
                tif_path, f"{x} {y}", "1").getOutput(0))
            railway_vals.append(val)
            print(f"    {loc}: {val:.6f}")
        except:
            print(f"    {loc}: NoData")

    road_vals = []
    print("  Current density — A66 road:")
    for loc, (x, y) in road_points.items():
        try:
            val = float(arcpy.management.GetCellValue(
                tif_path, f"{x} {y}", "1").getOutput(0))
            road_vals.append(val)
            print(f"    {loc}: {val:.6f}")
        except:
            print(f"    {loc}: NoData")

    if railway_vals and road_vals:
        rail_mean = sum(railway_vals) / len(railway_vals)
        road_mean = sum(road_vals) / len(road_vals)
        ratio     = rail_mean / road_mean if road_mean > 0 else 0
        validation_results[scenario_name] = {
            "railway_mean": rail_mean,
            "road_mean"   : road_mean,
            "ratio"       : ratio
        }
        print(f"\n  Mean current density:")
        print(f"    Eden Valley Railway : {rail_mean:.6f}")
        print(f"    A66 road            : {road_mean:.6f}")
        print(f"    Ratio (rail/road)   : {ratio:.2f}x")
        if rail_mean > road_mean:
            print(f"  [OK] Railway > A66 road — validates as "
                  f"wildlife corridor.")
        else:
            print(f"  [NOTE] Railway not > road at sampled points.")
            print(f"         Check visually in ArcGIS Pro — coordinate")
            print(f"         approximations may need refining.")

# Summary across all scenarios
if validation_results:
    print("\n  Validation summary across all scenarios:")
    print(f"  {'Scenario':<25} {'Railway':>10} {'Road':>10} {'Ratio':>8}")
    print(f"  {'-'*55}")
    for name, vals in validation_results.items():
        print(f"  {name:<25} {vals['railway_mean']:>10.6f} "
              f"{vals['road_mean']:>10.6f} {vals['ratio']:>7.2f}x")

# ---------------------------------------------------------------------------
# FINAL SUMMARY
# ---------------------------------------------------------------------------

print()
print("=" * 60)
print("SCRIPT 4b COMPLETE")
print(f"Finished: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 60)
print()

print("Output current density rasters:")
for name in SCENARIO_NAMES:
    tif = os.path.join(OUT_DIR, f"CurrentDensity_{name}.tif")
    if os.path.exists(tif):
        size_mb = os.path.getsize(tif) / 1024 / 1024
        print(f"  [OK] CurrentDensity_{name}.tif ({size_mb:.1f} MB)")
    else:
        print(f"  [MISSING] CurrentDensity_{name}.tif")

missing = [n for n in SCENARIO_NAMES
           if not os.path.exists(
               os.path.join(OUT_DIR, f"CurrentDensity_{n}.tif"))]
if missing:
    print(f"\n  [WARNING] Missing outputs: {missing}")
    print("  Check Circuitscape stdout/stderr logs above.")

print()
print("VISUAL QA — verify in ArcGIS Pro before Script 5:")
print("  [ ] Load CurrentDensity_S1_primary.tif — apply hot colour ramp")
print("  [ ] High current density corridors visible across study area")
print("  [ ] Eden Valley Railway — elevated current near Kirkby Thore")
print("  [ ] A66 road — low current — confirms barrier effect")
print("  [ ] River Eden corridor — moderate-high current")
print("  [ ] S1 vs S4 — similar spatial pattern (confirms high stability)")
print("  [ ] Overlay on GAD georeferenced sheets — confirm alignment")
print("  [ ] British Gypsum site — low current (high resistance confirmed)")
print()
print("VERIFIED PARAMETERS FOR DISSERTATION (Chapter 5):")
print(f"  Mode           : all-to-one")
print(f"  Focal nodes    : 15 woodland patches (8.3–64.1 ha)")
print(f"  Calculations   : 15 per scenario")
print(f"  Scenarios      : 4 (S1–S4)")
print(f"  Solver         : cg+amg")
print(f"  CRS            : EPSG:27700")
print(f"  Outputs        : {len(curmap_tifs)} current density GeoTIFFs")
print()
print("NEXT STEP: Script 5 — Validation and Map Export")
print("  Overlay current density vs:")
print("  - NH TR010062 disused Eden Valley Railway corridor")
print("  - Natural England Habitat Networks Combined Habitats England")
print("  - Cumbria LNRS APIB priority areas")
print("  Export all maps at 300 dpi, EPSG:27700, with scale bar,")
print("  north arrow, legend, and figure caption.")
print()
print("Push to GitHub:")
print("  https://github.com/Abhin-Nangari/A66-Connectivity-MSc")

arcpy.CheckInExtension("Spatial")
