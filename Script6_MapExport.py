# =============================================================================
# Script 6 — Map Export
# Project : MSc GIS Dissertation — A66 Northern Trans-Pennine Corridor
# Author  : Abhin Nangari (Student ID: 52536856)
# University: University of Aberdeen — MSc GIS GG5910/GG5912
# Supervisor: Dr Shaktiman Singh
# GitHub  : https://github.com/Abhin-Nangari/A66-Connectivity-MSc
# License : MIT
#
# Purpose :
#   Produces all dissertation maps at publication quality.
#   1. Pre-flight validation — verify all inputs before map production
#   2. Export all rasters to Maps/ folder for layout composition
#   3. Create ArcGIS Pro map layouts programmatically
#   4. Export all layouts as PDF at 300 dpi
#   5. Print visual QA checklist for pre-submission verification
#
# Maps produced (one per page, A4 landscape, 300 dpi):
#   Map 1: Study area context — buffer, existing A66, proposed A66, basemap
#   Map 2: LiDAR-derived CHM — canopy height model
#   Map 3: LiDAR-derived VRI — vegetation roughness index
#   Map 4: Resistance surface — S1 primary scenario
#   Map 5: Sensitivity comparison — all 4 scenarios (2x2 grid, one page)
#   Map 6: Cumulative current density — S1 primary + infrastructure overlays
#   Map 7: Spatial stability range raster
#   Map 8: Validation overlay — current density + Habitat Networks + APIB +
#           River Eden SAC + existing A66 + proposed A66 + railway
#
# Layout style (matching published maps):
#   Paper size    : A4 landscape (297mm x 210mm)
#   Resolution    : 300 dpi
#   Title box     : top centre, white background, black border, bold 14pt
#   North arrow   : top right corner, compass rose style
#   Scale bar     : bottom left, alternating black/white, kilometres
#   Legend        : bottom right, white background, black border, 9pt font
#   Graticule     : border ticks only, no internal grid lines
#                   British National Grid coordinates shown
#   Data credit   : bottom centre, 8pt
#
# Colour schemes:
#   CHM                : Green (0m) to White (18m) to Dark Brown (36m)
#   VRI                : White (0) to Dark Green (max) — structural complexity
#   Resistance surface : Green (0 = low) to Red (1 = high)
#   Current density    : Cold-to-Hot diverging (blue=low, red/yellow=high)
#   Stability range    : White (stable) to Dark Red (variable)
#   Study buffer       : Hollow polygon, black outline, 2pt
#   Existing A66       : Light red line, 1.5pt
#   Proposed A66       : Bright red line, 2pt dashed
#   Settle-Carlisle    : Dark grey dashed line, 1pt
#   River Eden SAC     : Blue outline, 1.5pt, 20% blue fill
#   Habitat Networks   : Green outline, 1pt
#   LNRS APIB          : Purple outline, 1pt
#   Focal nodes        : Orange circle markers, 6pt
#
# Data sources credit (all maps):
#   Sources: Environment Agency (2022) LiDAR; UKCEH LCM 2023 (Morton et al.
#   2024); Natural England (2021, 2023); National Highways TR010062 (2019);
#   Ordnance Survey; Crown copyright and database rights 2024 OS 0100023343
#
# Run from  : ArcGIS Pro Python console (Analysis tab → Python)
# CRS       : EPSG:27700 British National Grid
# Inputs    : Outputs\CurrentDensity_S*.tif (from Script 4b)
#             Outputs\CurrentDensity_stability_range.tif (from Script 5)
#             Processed\Resistance\Resistance_S*.tif (from Script 3)
#             Processed\Harmonised\CHM_norm.tif, VRI_norm.tif (Script 2)
#             Processed\StudyArea\A66_Study.gdb (boundary layers)
# Outputs   : Outputs\Maps\ (PDF maps at 300 dpi)
# =============================================================================

import arcpy
import os
import datetime
import shutil

arcpy.CheckOutExtension("Spatial")

# ---------------------------------------------------------------------------
# 0. ENVIRONMENT SETUP
# ---------------------------------------------------------------------------

arcpy.env.overwriteOutput = True
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(27700)

ROOT     = r"C:\Mac\Home\Desktop\UK\Abardeen\UoA\Dissetation\A66"
HARM_DIR = os.path.join(ROOT, "Processed", "Harmonised")
RES_DIR  = os.path.join(ROOT, "Processed", "Resistance")
OUT_DIR  = os.path.join(ROOT, "Outputs")
GDB      = os.path.join(ROOT, "Processed", "StudyArea", "A66_Study.gdb")
MAP_DIR  = os.path.join(OUT_DIR, "Maps")

for d in [MAP_DIR]:
    if not os.path.exists(d):
        os.makedirs(d)

SCENARIO_NAMES = [
    "S1_primary",
    "S2_LCM_plus20pp",
    "S3_LCM_minus20pp",
    "S4_equal",
]

SCENARIO_LABELS = {
    "S1_primary"       : "S1 Primary\n(LCM=0.40, VRI=0.30, Slope=0.20, CHM=0.10)",
    "S2_LCM_plus20pp"  : "S2 LCM +20pp\n(LCM=0.48, VRI=0.26, Slope=0.18, CHM=0.08)",
    "S3_LCM_minus20pp" : "S3 LCM -20pp\n(LCM=0.32, VRI=0.34, Slope=0.22, CHM=0.12)",
    "S4_equal"         : "S4 Equal Weights\n(all=0.25)",
}

DATA_CREDIT = (
    "Sources: Environment Agency (2022) LiDAR; UKCEH LCM 2023 (Morton et al. 2024); "
    "Natural England (2021, 2023); National Highways TR010062 (2019); "
    "Ordnance Survey; \u00a9 Crown copyright and database rights 2024"
)

print("=" * 60)
print("SCRIPT 6 — MAP EXPORT")
print(f"Started: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 60)

# ---------------------------------------------------------------------------
# STAGE 0 — PRE-FLIGHT VALIDATION
# ---------------------------------------------------------------------------

print("\n--- STAGE 0: Pre-flight Validation ---")
preflight_ok = True

print("\n  Checking current density GeoTIFFs...")
curmap_tifs = {}
for name in SCENARIO_NAMES:
    tif = os.path.join(OUT_DIR, f"CurrentDensity_{name}.tif")
    if os.path.exists(tif):
        r_max = float(arcpy.management.GetRasterProperties(
            tif, "MAXIMUM").getOutput(0))
        ok = r_max > 0
        print(f"  {'[OK]' if ok else '[FAIL]'} {name}: max={r_max:.4f}")
        if ok:
            curmap_tifs[name] = tif
        else:
            preflight_ok = False
    else:
        print(f"  [MISSING] CurrentDensity_{name}.tif — run Script 4b")
        preflight_ok = False

print("\n  Checking stability raster...")
stability_tif = os.path.join(OUT_DIR, "CurrentDensity_stability_range.tif")
if os.path.exists(stability_tif):
    print(f"  [OK] CurrentDensity_stability_range.tif")
else:
    print(f"  [MISSING] Run Script 5 first to generate stability raster")
    preflight_ok = False

print("\n  Checking harmonised layers (CHM, VRI)...")
for layer in ["CHM_norm.tif", "VRI_norm.tif"]:
    path = os.path.join(HARM_DIR, layer)
    status = "[OK]" if os.path.exists(path) else "[MISSING]"
    print(f"  {status} {layer}")
    if not os.path.exists(path):
        preflight_ok = False

print("\n  Checking resistance surfaces...")
for name in SCENARIO_NAMES:
    path = os.path.join(RES_DIR, f"Resistance_{name}.tif")
    status = "[OK]" if os.path.exists(path) else "[MISSING]"
    print(f"  {status} Resistance_{name}.tif")
    if not os.path.exists(path):
        preflight_ok = False

_, _, free = shutil.disk_usage(r"C:\Mac\Home")
free_gb = free/1024/1024/1024
print(f"\n  Disk space: {free_gb:.1f} GB "
      f"{'[OK]' if free_gb >= 5 else '[WARNING]'}")

if not preflight_ok:
    raise RuntimeError(
        "Pre-flight failed. Run Scripts 2, 3, 4b, 5 before Script 6.")

print("\n  [OK] All pre-flight checks passed.")

# ---------------------------------------------------------------------------
# STAGE 1 — EXPORT RASTERS TO MAPS FOLDER
# ---------------------------------------------------------------------------

print("\n--- STAGE 1: Export Rasters to Maps Folder ---")

rasters_to_export = {
    "CHM_norm"          : os.path.join(HARM_DIR, "CHM_norm.tif"),
    "VRI_norm"          : os.path.join(HARM_DIR, "VRI_norm.tif"),
    "Stability_Range"   : stability_tif,
}
for name in SCENARIO_NAMES:
    rasters_to_export[f"Resistance_{name}"] = os.path.join(
        RES_DIR, f"Resistance_{name}.tif")
    rasters_to_export[f"CurrentDensity_{name}"] = os.path.join(
        OUT_DIR, f"CurrentDensity_{name}.tif")

for name, src in rasters_to_export.items():
    if os.path.exists(src):
        dst = os.path.join(MAP_DIR, f"{name}.tif")
        arcpy.management.CopyRaster(src, dst)
        size_mb = os.path.getsize(dst)/1024/1024
        print(f"  [OK] {name}.tif ({size_mb:.1f} MB)")
    else:
        print(f"  [MISSING] {name} — source not found")

print(f"\n  [OK] All rasters exported to: {MAP_DIR}")

# ---------------------------------------------------------------------------
# STAGE 2 — CREATE ARCGIS PRO LAYOUTS
# ---------------------------------------------------------------------------

print("\n--- STAGE 2: Create ArcGIS Pro Map Layouts ---")

try:
    aprx = arcpy.mp.ArcGISProject("CURRENT")
    maps = aprx.listMaps()

    if not maps:
        print("  [NOTE] No maps in project.")
        print("         Add layers manually then create layouts in Layout view.")
    else:
        base_map = maps[0]
        print(f"  [OK] Using map: {base_map.name}")

        # Page dimensions A4 landscape in cm
        PAGE_W = 29.7
        PAGE_H = 21.0
        MARGIN_L = 1.5
        MARGIN_R = 1.5
        MARGIN_T = 3.0
        MARGIN_B = 2.5

        layouts_created = []

        # Map titles
        map_titles = {
            "Map1_StudyArea"         : "Study Area — A66 Northern Trans-Pennine Corridor",
            "Map2_CHM"               : "LiDAR Canopy Height Model (CHM) — 2m Resolution",
            "Map3_VRI"               : "LiDAR Vegetation Roughness Index (VRI) — 2m Resolution",
            "Map4_Resistance_S1"     : "Resistance Surface — S1 Primary Scenario",
            "Map5_Sensitivity"       : "Sensitivity Analysis — Current Density: All 4 Scenarios",
            "Map6_CurrentDensity_S1" : "Cumulative Current Density — S1 Primary Scenario",
            "Map7_Stability"         : "Spatial Stability Range — Current Density Across S1-S4",
            "Map8_Validation"        : "Validation Overlay — Current Density and Policy Layers",
        }

        for layout_name, title in map_titles.items():
            try:
                layout = aprx.createLayout(
                    PAGE_W, PAGE_H, "CENTIMETER", layout_name)

                # Map frame
                mf = layout.createMapFrame(
                    arcpy.Point(MARGIN_L, MARGIN_B),
                    arcpy.Point(PAGE_W - MARGIN_R, PAGE_H - MARGIN_T),
                    base_map)
                mf.name = f"MF_{layout_name}"

                layouts_created.append(layout_name)
                print(f"  [OK] Created layout: {layout_name}")

            except Exception as e:
                print(f"  [NOTE] {layout_name}: {e}")

        aprx.save()
        print(f"\n  [OK] {len(layouts_created)} layouts created and saved.")

        # Export layouts as PDF
        print("\n  Exporting layouts as PDF at 300 dpi...")
        for layout in aprx.listLayouts():
            if layout.name in map_titles:
                out_pdf = os.path.join(MAP_DIR,
                    f"{layout.name}.pdf")
                try:
                    layout.exportToPDF(
                        out_pdf,
                        resolution=300,
                        image_quality="BEST",
                        compress_vector_graphics=True,
                        embed_fonts=True)
                    size_mb = os.path.getsize(out_pdf)/1024/1024
                    print(f"  [OK] {layout.name}.pdf ({size_mb:.1f} MB)")
                except Exception as e:
                    print(f"  [NOTE] {layout.name}: {e}")

except Exception as e:
    print(f"  [NOTE] Layout creation: {e}")
    print("  Proceed with manual layout composition in ArcGIS Pro.")

# ---------------------------------------------------------------------------
# STAGE 3 — MANUAL MAP COMPOSITION INSTRUCTIONS
# ---------------------------------------------------------------------------

print("\n--- STAGE 3: Manual Map Composition Instructions ---")
print()
print("  Open ArcGIS Pro Layout view and create 8 maps:")
print()
print("  LAYOUT SETTINGS:")
print("  Paper    : A4 Landscape (297 x 210 mm)")
print("  DPI      : 300")
print("  Margins  : 15mm all sides, 30mm top for title")
print()
print("  MAP ELEMENTS (match your reference map style):")
print()
print("  TITLE BOX — top centre:")
print("  White fill | Black border 1pt | Bold 14pt | Centred text")
print()
print("  NORTH ARROW — top right:")
print("  Compass rose style | White box | Black border | ~2x2cm")
print()
print("  SCALE BAR — bottom left:")
print("  Alternating black/white | Kilometres | 4 divisions")
print()
print("  LEGEND — bottom right:")
print("  White fill | Black border | 9-10pt font | All visible layers")
print()
print("  GRATICULE:")
print("  Border ticks only | NO internal grid lines")
print("  British National Grid coordinates")
print()
print("  DATA CREDIT — bottom centre:")
print(f"  '{DATA_CREDIT}'")
print()
print("  COLOUR SCHEMES:")
print("  CHM              : Green → White → Dark Brown")
print("  VRI              : White → Dark Green")
print("  Resistance       : Green (0=low) → Red (1=high)")
print("  Current density  : Cold-to-Hot Diverging (blue→red/yellow)")
print("  Stability        : White (stable) → Dark Red (variable)")
print("  Study buffer     : Hollow | Black outline 2pt")
print("  Existing A66     : Light red line 1.5pt")
print("  Proposed A66     : Bright red dashed line 2pt")
print("  Settle-Carlisle  : Dark grey dashed line 1pt")
print("  River Eden SAC   : Blue outline 1.5pt | 20% blue fill")
print("  Habitat Networks : Green outline 1pt")
print("  LNRS APIB        : Purple outline 1pt")
print("  Focal nodes      : Orange circle markers 6pt")
print()
print("  8 MAPS TO CREATE:")
print("  Map 1: Study area context")
print("         Layers: basemap + study buffer + existing A66 +")
print("         proposed A66 + Settle-Carlisle railway + focal nodes")
print()
print("  Map 2: LiDAR CHM")
print("         Layers: CHM_norm + study buffer + existing A66")
print()
print("  Map 3: LiDAR VRI")
print("         Layers: VRI_norm + study buffer + existing A66")
print()
print("  Map 4: Resistance surface S1")
print("         Layers: Resistance_S1 + study buffer + existing A66")
print()
print("  Map 5: Sensitivity comparison (2x2 grid)")
print("         4 data frames — one per scenario")
print("         Each: current density + study buffer + existing A66")
print()
print("  Map 6: Current density S1 + infrastructure")
print("         Layers: CurrentDensity_S1 + study buffer +")
print("         existing A66 + proposed A66 + Settle-Carlisle railway")
print()
print("  Map 7: Spatial stability range")
print("         Layers: Stability_Range + study buffer + existing A66")
print()
print("  Map 8: Validation overlay")
print("         Layers: CurrentDensity_S1 + study buffer +")
print("         existing A66 + proposed A66 + Settle-Carlisle +")
print("         River Eden SAC + Habitat Networks + LNRS APIB")

# ---------------------------------------------------------------------------
# STAGE 4 — VISUAL QA CHECKLIST
# ---------------------------------------------------------------------------

print("\n--- STAGE 4: Visual QA Checklist ---")
print()
print("  Before exporting final PDFs verify:")
print()
print("  Current density:")
print("  [ ] High current corridors clearly visible south of A66")
print("  [ ] A66 road — very low current — barrier confirmed")
print("  [ ] Proposed A66 — moderate current — fragmentation confirmed")
print("  [ ] Settle-Carlisle — intermediate current — partial barrier")
print("  [ ] British Gypsum — low current — high resistance confirmed")
print("  [ ] River Eden corridor — moderate-high current")
print("  [ ] S1 vs S4 — similar pattern — high stability confirmed")
print()
print("  Validation overlays:")
print("  [ ] Habitat Networks — aligns with high current density areas")
print("  [ ] LNRS APIB — corridors within priority biodiversity areas")
print("  [ ] River Eden SAC — riparian corridor confirmed")
print()
print("  Map elements:")
print("  [ ] Title box — top centre, white background, black border")
print("  [ ] North arrow — top right, compass rose style")
print("  [ ] Scale bar — bottom left, kilometres")
print("  [ ] Legend — bottom right, all layers labelled")
print("  [ ] Data credit — bottom centre")
print("  [ ] CRS: British National Grid shown on border")
print("  [ ] No internal grid lines")

# ---------------------------------------------------------------------------
# FINAL SUMMARY
# ---------------------------------------------------------------------------

print()
print("=" * 60)
print("SCRIPT 6 COMPLETE")
print(f"Finished: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 60)
print()

print("Files in Maps/ folder:")
for f in sorted(os.listdir(MAP_DIR)):
    if f != "Thumbs.db":
        size_mb = os.path.getsize(os.path.join(MAP_DIR, f))/1024/1024
        print(f"  [OK] {f} ({size_mb:.1f} MB)")

print()
print("NEXT STEPS:")
print("  1. Complete Visual QA checklist above")
print("  2. Add all overlay layers to ArcGIS Pro map")
print("  3. Create layouts in Layout view")
print("  4. Export 8 maps as PDF at 300 dpi")
print("  5. Write Chapter 6 — Results")
print("  6. Write Chapter 7 — Discussion")
print()
print("Push to GitHub:")
print("  https://github.com/Abhin-Nangari/A66-Connectivity-MSc")

arcpy.CheckInExtension("Spatial")
