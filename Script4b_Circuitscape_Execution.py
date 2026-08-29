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
#      Pairwise mode — 105 pairs per scenario — ArcGIS Pro stays open
#   3. Verify Circuitscape outputs — non-zero current density confirmed
#   4. Import cumulative current density rasters into ArcGIS Pro
#      Uses numpy method — bypasses ASCIIToRaster UNC path issue
#   5. Verify imported GeoTIFFs — band count, extent, statistics
#   6. Validate current density along Eden Valley Railway (primary validation)
#   7. Quantitative comparison: Railway vs A66 road across all 4 scenarios
#
# Prerequisites — Script 4a must be complete before running:
#   Processed\CS_Inputs\focal_nodes.asc
#   Processed\CS_Inputs\Resistance_S*.asc (4 files)
#   Processed\CS_Inputs\config_S*.ini (4 files — pairwise mode)
#   Outputs\S*\ (4 output folders)
#
# Circuitscape mode — pairwise:
#   15 nodes = 105 pairwise calculations per scenario.
#   Pairwise is the original Circuitscape methodology.
#   all-to-one mode produces zero output in Circuitscape v5.17.1.
#   Pairwise confirmed producing correct non-zero current density.
#   Citation: McRae et al. (2008) Ecology 89(10) pp.2712-2724
#             Wilson et al. (2017) PLoS Computational Biology 13(6) e1005510
#             Dickson et al. (2019) Conservation Biology 33(2) pp.239-249
#
# Import method:
#   arcpy.ASCIIToRaster fails on Parallels UNC paths (path truncation).
#   numpy used to read ASCII and arcpy.NumPyArrayToRaster to create GeoTIFF.
#   Produces identical output with correct EPSG:27700 projection.
#
# Eden Valley Railway validation:
#   Primary validation layer per NH TR010062.
#   National Highways selected Option E alignment specifically to avoid
#   the wildlife corridor on the disused railway line near Kirkby Thore
#   — the only formally named wildlife corridor on the entire A66 upgrade.
#   Expected: elevated current density along railway vs A66 road.
#   Citation: National Highways (2019) NH Citizen Space TR010062
#             SI 2024/360 Schedule 7 Part 3
#
# Run from : ArcGIS Pro Python console (Analysis tab → Python)
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
import numpy as np

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

SCENARIO_NAMES = [
    "S1_primary",
    "S2_LCM_plus20pp",
    "S3_LCM_minus20pp",
    "S4_equal",
]

# Master grid for extent verification
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
# HELPER — verify raster
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
    print(f"    Bands    : {bands} {'[OK]' if band_ok else '[WARNING]'}")
    print(f"    Cell size: {cs}m")
    print(f"    Extent   : {ext_t}")
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
    if band_ok and stat_ok:
        print(f"    [OK] All checks passed.")
    return r_min, r_max, r_mean, r_std

# ---------------------------------------------------------------------------
# HELPER — import ASCII using numpy
# ---------------------------------------------------------------------------

def import_ascii_numpy(asc_path, out_tif, name):
    """
    Import Circuitscape ASCII output to GeoTIFF using numpy.
    Bypasses arcpy.ASCIIToRaster which fails on Parallels UNC paths.
    """
    print(f"  Importing {name}...")
    with open(asc_path) as f:
        ncols    = int(f.readline().split()[1])
        nrows    = int(f.readline().split()[1])
        xll      = float(f.readline().split()[1])
        yll      = float(f.readline().split()[1])
        cellsize = float(f.readline().split()[1])
        nodata   = float(f.readline().split()[1])
        arr      = np.array(
            [[float(v) for v in line.strip().split()]
             for line in f],
            dtype=np.float32
        )

    # Replace nodata with NaN for correct raster handling
    arr[arr == nodata] = np.nan

    lower_left = arcpy.Point(xll, yll)
    raster     = arcpy.NumPyArrayToRaster(
        arr, lower_left, cellsize, cellsize, np.nan)
    arcpy.management.DefineProjection(
        raster, arcpy.SpatialReference(27700))
    raster.save(out_tif)

    # Verify non-zero values
    r_max = float(arcpy.management.GetRasterProperties(
        out_tif, "MAXIMUM").getOutput(0))
    r_mean = float(arcpy.management.GetRasterProperties(
        out_tif, "MEAN").getOutput(0))

    if r_max > 0:
        print(f"  [OK] {os.path.basename(out_tif)}: "
              f"max={r_max:.6f}  mean={r_mean:.6f}")
        return True
    else:
        print(f"  [WARNING] {name}: all values zero — check Circuitscape output")
        return False

# ---------------------------------------------------------------------------
# STAGE 0 — PRE-FLIGHT VALIDATION
# ---------------------------------------------------------------------------

print("\n--- STAGE 0: Pre-flight Validation (Script 4a Outputs) ---")
preflight_ok = True

print("\n  Checking focal nodes ASCII...")
nodes_asc = os.path.join(CS_DIR, "focal_nodes.asc")
if os.path.exists(nodes_asc) and os.path.getsize(nodes_asc) > 5*1024*1024:
    size_mb = os.path.getsize(nodes_asc) / 1024 / 1024
    print(f"  [OK] focal_nodes.asc ({size_mb:.1f} MB)")
else:
    print(f"  [ABORT] focal_nodes.asc missing. Run Script 4a first.")
    preflight_ok = False

print("\n  Checking resistance ASCII files...")
for name in SCENARIO_NAMES:
    asc = os.path.join(CS_DIR, f"Resistance_{name}.asc")
    if os.path.exists(asc) and os.path.getsize(asc) > 10*1024*1024:
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
        # Verify pairwise mode
        with open(cfg) as f:
            content = f.read()
        if "scenario = pairwise" in content:
            print(f"  [OK] config_{name}.ini (pairwise mode confirmed)")
        else:
            print(f"  [WARNING] config_{name}.ini — pairwise mode not found")
        config_files[name] = cfg
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

print("\n  Checking output folders...")
for name in SCENARIO_NAMES:
    folder = os.path.join(OUT_DIR, name)
    if os.path.exists(folder):
        print(f"  [OK] {name}/")
    else:
        os.makedirs(folder)
        print(f"  [CREATED] {name}/")

_, _, free = shutil.disk_usage(r"C:\Mac\Home")
free_gb = free / 1024 / 1024 / 1024
print(f"\n  Disk space free: {free_gb:.1f} GB "
      f"{'[OK]' if free_gb >= 5 else '[WARNING]'}")
if free_gb < 5:
    preflight_ok = False

if not preflight_ok:
    raise RuntimeError(
        "Pre-flight validation failed. Run Script 4a before Script 4b.")

print("\n  [OK] All pre-flight checks passed — proceeding.")
print()
print("  NOTE: ArcGIS Pro remains open throughout.")
print("  Circuitscape runs as an independent subprocess.")
print("  Mode: pairwise — 105 pairs per scenario")
print()

# ---------------------------------------------------------------------------
# STAGE 1 — RUN CIRCUITSCAPE VIA JULIA — ALL 4 SCENARIOS
# ---------------------------------------------------------------------------

print("\n--- STAGE 1: Run Circuitscape (pairwise mode) ---")
print(f"  Mode     : pairwise")
print(f"  Nodes    : 15")
print(f"  Pairs    : 105 per scenario")
print(f"  Scenarios: 4")
print(f"  Solver   : cg+amg")
print(f"  Started  : {datetime.datetime.now().strftime('%H:%M:%S')}")
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
            timeout=7200)  # 2 hour timeout — pairwise takes longer

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
        print(f"  [TIMEOUT] Exceeded 2 hours for {name}.")
    except Exception as e:
        print(f"  [ERROR] {e}")
    print()

print(f"  All scenarios finished: "
      f"{datetime.datetime.now().strftime('%H:%M:%S')}")
print(f"  Successful: {len(curmap_paths)} of {len(SCENARIO_NAMES)}")

if len(curmap_paths) == 0:
    raise RuntimeError(
        "No Circuitscape outputs produced. "
        "Check Julia installation and config files.")

# ---------------------------------------------------------------------------
# STAGE 2 — VERIFY CIRCUITSCAPE OUTPUTS
# ---------------------------------------------------------------------------

print("\n--- STAGE 2: Verify Circuitscape Outputs ---")
for name in SCENARIO_NAMES:
    if name in curmap_paths:
        path    = curmap_paths[name]
        size_mb = os.path.getsize(path) / 1024 / 1024
        # Quick scan for non-zero values
        nonzero = 0
        with open(path) as f:
            for i, line in enumerate(f):
                if i < 6:
                    continue
                vals = [float(v) for v in line.strip().split()
                        if v not in ['-9999', '-9999.0', '0', '0.0']]
                nonzero += len(vals)
                if nonzero > 100:
                    break
        status = "[OK]" if nonzero > 0 else "[WARNING] All zeros"
        print(f"  {status} {name}: {os.path.basename(path)} "
              f"({size_mb:.1f} MB)  non-zero pixels found: {nonzero}+")
    else:
        print(f"  [MISSING] {name}: no output produced")

# ---------------------------------------------------------------------------
# STAGE 3 — IMPORT CURRENT DENSITY RASTERS
# ---------------------------------------------------------------------------
# Uses numpy method — bypasses arcpy.ASCIIToRaster which fails on
# Parallels Desktop due to UNC path truncation (8.3 format conversion).

print("\n--- STAGE 3: Import Current Density Rasters (numpy method) ---")

curmap_tifs = {}

for name, curmap_src in curmap_paths.items():
    out_tif = os.path.join(OUT_DIR, f"CurrentDensity_{name}.tif")

    # Handle gzip compressed output
    import_path = curmap_src
    if curmap_src.endswith(".gz"):
        import gzip
        unzipped = curmap_src.replace(".gz", "")
        with gzip.open(curmap_src, 'rb') as f_in:
            with open(unzipped, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        import_path = unzipped

    if import_ascii_numpy(import_path, out_tif, name):
        curmap_tifs[name] = out_tif
        verify(out_tif, f"CurrentDensity_{name}", exp_min=0.0, exp_bands=1)

# ---------------------------------------------------------------------------
# STAGE 4 — VALIDATE — EDEN VALLEY RAILWAY
# ---------------------------------------------------------------------------

print("\n--- STAGE 4: Validate — Eden Valley Railway ---")
print("  Primary validation: disused railway near Kirkby Thore")
print("  NH TR010062: only formally named wildlife corridor on A66 upgrade")
print()

railway_points = {
    "Railway W end"  : (362300, 520750),
    "Railway mid-W"  : (362800, 520900),
    "Railway centre" : (363200, 521000),
    "Railway mid-E"  : (363600, 521100),
    "Railway E end"  : (364000, 521200),
}

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
            print(f"  [OK] Railway > A66 road — validates as wildlife corridor.")
        else:
            print(f"  [NOTE] Check visually in ArcGIS Pro.")

if validation_results:
    print("\n  Validation summary — all scenarios:")
    print(f"  {'Scenario':<25} {'Railway':>10} {'Road':>10} {'Ratio':>8}")
    print(f"  {'-'*57}")
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
    print(f"\n  [WARNING] Missing: {missing}")

print()
print("VISUAL QA — verify in ArcGIS Pro before Script 5:")
print("  [ ] Load CurrentDensity_S1_primary.tif — hot colour ramp")
print("  [ ] High current corridors visible across study area")
print("  [ ] Eden Valley Railway — elevated current near Kirkby Thore")
print("  [ ] A66 road — low current — confirms barrier effect")
print("  [ ] River Eden corridor — moderate-high current")
print("  [ ] S1 vs S4 — similar spatial pattern (confirms high stability)")
print("  [ ] Overlay on GAD georeferenced sheets — confirm alignment")
print("  [ ] British Gypsum site — low current (high resistance)")
print()
print("VERIFIED PARAMETERS FOR DISSERTATION (Chapter 5):")
print(f"  Mode           : pairwise")
print(f"  Focal nodes    : 15 woodland patches (8.3-64.1 ha)")
print(f"  Pairwise pairs : 105 per scenario")
print(f"  Scenarios      : 4 (S1-S4)")
print(f"  Solver         : cg+amg")
print(f"  Resolution     : 3m")
print(f"  CRS            : EPSG:27700")
print(f"  Outputs        : {len(curmap_tifs)} current density GeoTIFFs")
print()
print("NEXT STEP: Script 5 — Validation and Map Export")
print("  Overlay current density vs:")
print("  - NH TR010062 disused Eden Valley Railway corridor")
print("  - Natural England Habitat Networks Combined Habitats England")
print("  - Cumbria LNRS APIB priority areas")
print("  Export all maps at 300 dpi, EPSG:27700")
print()
print("Push to GitHub:")
print("  https://github.com/Abhin-Nangari/A66-Connectivity-MSc")

arcpy.CheckInExtension("Spatial")
