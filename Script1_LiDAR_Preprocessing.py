# =============================================================================
# Script 1 — LiDAR Preprocessing
# Project : MSc GIS Dissertation — A66 Northern Trans-Pennine Corridor
# Author  : Abhin Nangari (Student ID: 52536856)
# University: University of Aberdeen — MSc GIS GG5910/GG5912
# Supervisor: Dr Shaktiman Singh
# GitHub  : https://github.com/Abhin-Nangari/A66-Connectivity-MSc
# License : MIT
#
# Purpose : 
#   1. Mosaic all EA LiDAR DTM 1m tiles  → DTM_mosaic.tif
#   2. Mosaic all EA LiDAR DSM 1m tiles  → DSM_mosaic.tif
#   3. CHM = DSM − DTM                   → CHM_1m.tif
#   4. VRI = Focal Statistics (3×3 STD)  → VRI_1m.tif
#   5. Survey Index QA — log acquisition dates, flag pre-2010 tiles
#
# Run from: ArcGIS Pro Python console (Analysis tab → Python)
# Python  : C:\ArcGIS\Sem1\GG5567\Assessment2\DroughtAnalysis_Scotland.venv\Scripts\python.exe
# CRS     : EPSG:27700 British National Grid (all outputs)
# Outputs : C:\Mac\Home\Desktop\UK\Abardeen\UoA\Dissetation\A66\Processed\LiDAR\
# =============================================================================

import arcpy
import os
import glob
import datetime

# ---------------------------------------------------------------------------
# 0. ENVIRONMENT SETUP
# ---------------------------------------------------------------------------

arcpy.env.overwriteOutput = True

# Root paths — all Windows-side paths via Parallels C:\Mac\Home\... convention
ROOT        = r"C:\Mac\Home\Desktop\UK\Abardeen\UoA\Dissetation\A66"
RAW_DTM     = os.path.join(ROOT, "Raw", "LiDaR", "LiDAR_DTM")
RAW_DSM     = os.path.join(ROOT, "Raw", "LiDaR", "LiDAR_DSM")
RAW_INDEX   = os.path.join(ROOT, "Raw", "LiDaR", "LIDAR_DSM_Time_Stamped_Extents1",
                           "LIDAR_DSM_Time_Stamped_Extents")
OUT_DIR     = os.path.join(ROOT, "Processed", "LiDAR")
GDB         = os.path.join(ROOT, "Processed", "StudyArea", "A66_Study.gdb")
STUDY_BUFFER = os.path.join(GDB, "StudyArea_2km_buffer")

# Output file paths
DTM_MOSAIC  = os.path.join(OUT_DIR, "DTM_mosaic.tif")
DSM_MOSAIC  = os.path.join(OUT_DIR, "DSM_mosaic.tif")
CHM_1M      = os.path.join(OUT_DIR, "CHM_1m.tif")
VRI_1M      = os.path.join(OUT_DIR, "VRI_1m.tif")

# Ensure output directory exists
if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)
    print(f"Created output directory: {OUT_DIR}")

# Set coordinate system for all outputs — EPSG:27700 British National Grid
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(27700)

print("=" * 60)
print("SCRIPT 1 — LiDAR PREPROCESSING")
print(f"Started: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 60)

# ---------------------------------------------------------------------------
# 1. SURVEY INDEX QA — log tile dates, flag pre-2010 tiles
# ---------------------------------------------------------------------------
# Reads the EA Survey Index shapefile to log acquisition dates for all tiles.
# Flags any tiles acquired before 2010 as potentially lower quality.

print("\n--- STAGE 1: Survey Index QA ---")

# Find the survey index shapefile (only one .shp expected in this folder)
index_shps = glob.glob(os.path.join(RAW_INDEX, "*.shp"))

if not index_shps:
    print("  WARNING: No Survey Index shapefile found in:")
    print(f"  {RAW_INDEX}")
    print("  Skipping Survey Index QA — proceed with caution.")
else:
    index_shp = index_shps[0]
    print(f"  Survey Index shapefile: {os.path.basename(index_shp)}")

    # Field names confirmed from pre-flight inspection of the shapefile:
    #   lidar_name — String — tile identifier
    #   year       — String — acquisition year (already clean, no parsing needed)
    #   sdflown    — Date   — survey start date (available if year is ever null)
    TILE_FIELD = "lidar_name"
    YEAR_FIELD = "year"

    print(f"  Tile field : {TILE_FIELD}")
    print(f"  Year field : {YEAR_FIELD}")
    print()

    early_tiles = []
    total_tiles = 0
    min_year = 9999
    max_year = 0

    with arcpy.da.SearchCursor(index_shp, [TILE_FIELD, YEAR_FIELD]) as cursor:
        for row in cursor:
            total_tiles += 1
            tile_name = row[0] if row[0] else f"Tile_{total_tiles}"
            year_val  = row[1]

            year = None
            if year_val:
                try:
                    year = int(str(year_val).strip()[:4])
                except ValueError:
                    pass

            if year:
                min_year = min(min_year, year)
                max_year = max(max_year, year)
                if year < 2010:
                    early_tiles.append((tile_name, year))
                    print(f"  [FLAG] Pre-2010 tile: {tile_name} — year {year}")
                else:
                    print(f"  [OK]   {tile_name} — year {year}")

    print()
    print(f"  Total tiles in Survey Index : {total_tiles}")
    print(f"  Acquisition year range      : {min_year}–{max_year}")
    print(f"  Pre-2010 tiles flagged      : {len(early_tiles)}")
    if early_tiles:
        print("  NOTE: Pre-2010 tiles may have lower point density.")
        print("        Note these tile names in dissertation limitations.")

print("\n  Survey Index QA complete.")

# ---------------------------------------------------------------------------
# 2. MOSAIC DTM TILES → DTM_mosaic.tif
# ---------------------------------------------------------------------------
# Collects all .tif files in the DTM folder and mosaics them to a single
# raster. Uses MEAN mosaic operator (appropriate where tiles overlap slightly
# at edges) and FLOAT32 pixel type to preserve sub-metre elevation precision.

print("\n--- STAGE 2: Mosaic DTM tiles ---")

dtm_tiles = glob.glob(os.path.join(RAW_DTM, "*.tif"))
if not dtm_tiles:
    dtm_tiles = glob.glob(os.path.join(RAW_DTM, "**", "*.tif"), recursive=True)

print(f"  DTM tiles found: {len(dtm_tiles)}")
if len(dtm_tiles) == 0:
    raise FileNotFoundError(
        f"No DTM .tif files found in {RAW_DTM}. "
        "Check that all tiles are extracted."
    )

# Build semicolon-delimited input string for MosaicToNewRaster
dtm_input_str = ";".join(dtm_tiles)

print(f"  Mosaicking {len(dtm_tiles)} DTM tiles...")
print(f"  Output: {DTM_MOSAIC}")

arcpy.management.MosaicToNewRaster(
    input_rasters        = dtm_input_str,
    output_location      = OUT_DIR,
    raster_dataset_name_with_extension = "DTM_mosaic.tif",
    coordinate_system_for_the_raster  = arcpy.SpatialReference(27700),
    pixel_type           = "32_BIT_FLOAT",
    cellsize             = 1,          # 1m native resolution
    number_of_bands      = 1,
    mosaic_method        = "MEAN",
    mosaic_colormap_mode = "FIRST"
)

print("  DTM mosaic complete.")

# Quick check — describe the output
dtm_desc = arcpy.Describe(DTM_MOSAIC)
print(f"  DTM mosaic extent   : {dtm_desc.extent}")
print(f"  DTM mosaic cell size: {dtm_desc.meanCellWidth}m × {dtm_desc.meanCellHeight}m")

# ---------------------------------------------------------------------------
# 3. MOSAIC DSM TILES → DSM_mosaic.tif
# ---------------------------------------------------------------------------
# Same approach as DTM. DSM First Return tiles include NY22ne which was
# confirmed extracted. MEAN operator used — appropriate for edge overlaps.

print("\n--- STAGE 3: Mosaic DSM tiles ---")

dsm_tiles = glob.glob(os.path.join(RAW_DSM, "*.tif"))
if not dsm_tiles:
    dsm_tiles = glob.glob(os.path.join(RAW_DSM, "**", "*.tif"), recursive=True)

print(f"  DSM tiles found: {len(dsm_tiles)}")
if len(dsm_tiles) == 0:
    raise FileNotFoundError(
        f"No DSM .tif files found in {RAW_DSM}. "
        "Check that all tiles including NY22ne are extracted."
    )

# Verify NY22ne is present — it was the last tile confirmed extracted
ny22ne_present = any("NY22ne" in os.path.basename(t) or "ny22ne" in os.path.basename(t).lower()
                     for t in dsm_tiles)
if ny22ne_present:
    print("  [OK] NY22ne tile confirmed present in DSM folder.")
else:
    print("  [WARNING] NY22ne tile NOT found in DSM folder — check extraction.")

dsm_input_str = ";".join(dsm_tiles)

print(f"  Mosaicking {len(dsm_tiles)} DSM tiles...")
print(f"  Output: {DSM_MOSAIC}")

arcpy.management.MosaicToNewRaster(
    input_rasters        = dsm_input_str,
    output_location      = OUT_DIR,
    raster_dataset_name_with_extension = "DSM_mosaic.tif",
    coordinate_system_for_the_raster  = arcpy.SpatialReference(27700),
    pixel_type           = "32_BIT_FLOAT",
    cellsize             = 1,
    number_of_bands      = 1,
    mosaic_method        = "MEAN",
    mosaic_colormap_mode = "FIRST"
)

print("  DSM mosaic complete.")

dsm_desc = arcpy.Describe(DSM_MOSAIC)
print(f"  DSM mosaic extent   : {dsm_desc.extent}")
print(f"  DSM mosaic cell size: {dsm_desc.meanCellWidth}m × {dsm_desc.meanCellHeight}m")

# ---------------------------------------------------------------------------
# 4. CANOPY HEIGHT MODEL — CHM = DSM − DTM → CHM_1m.tif
# ---------------------------------------------------------------------------
# Raster Calculator: CHM = DSM_mosaic − DTM_mosaic
# Represents above-ground vegetation/structure height.
# Expected range: 0m (bare ground/water) to ~30m (mature broadleaf woodland).
# Negative values can occur where DSM < DTM due to processing artefacts —
# these are clamped to 0 using Con().

print("\n--- STAGE 4: Canopy Height Model (CHM = DSM − DTM) ---")
print(f"  Output: {CHM_1M}")

# Enable Spatial Analyst extension
arcpy.CheckOutExtension("Spatial")

from arcpy.sa import (
    Raster, FocalStatistics, NbrRectangle, Con, IsNull, Float
)

dtm_raster = Raster(DTM_MOSAIC)
dsm_raster = Raster(DSM_MOSAIC)

# Calculate raw CHM
chm_raw = dsm_raster - dtm_raster

# Clamp negative values to 0 — negative CHM is physically meaningless
# and typically arises from LiDAR processing artefacts in water bodies.
chm_clamped = Con(chm_raw < 0, 0, chm_raw)

chm_clamped.save(CHM_1M)

print("  CHM calculation complete.")

# Report basic statistics
chm_result = arcpy.management.GetRasterProperties(CHM_1M, "MINIMUM")
chm_min = float(chm_result.getOutput(0))
chm_result = arcpy.management.GetRasterProperties(CHM_1M, "MAXIMUM")
chm_max = float(chm_result.getOutput(0))
chm_result = arcpy.management.GetRasterProperties(CHM_1M, "MEAN")
chm_mean = float(chm_result.getOutput(0))

print(f"  CHM statistics — Min: {chm_min:.2f}m  Max: {chm_max:.2f}m  Mean: {chm_mean:.2f}m")

# Visual QA guidance
print()
print("  VISUAL QA — CHM:")
print("  Expected: open fields = 0–2m, scrub = 2–5m, woodland = 5–20m+")
if chm_max < 3:
    print("  [WARNING] CHM max < 3m — possible DTM/DSM mismatch or wrong tile set.")
elif chm_max > 40:
    print("  [WARNING] CHM max > 40m — possible building artefacts or outliers.")
else:
    print("  [OK] CHM range looks plausible for upland farmland/woodland mosaic.")

# ---------------------------------------------------------------------------
# 5. VEGETATION ROUGHNESS INDEX — VRI = Focal Statistics 3×3 STD → VRI_1m.tif
# ---------------------------------------------------------------------------
# VRI captures structural complexity at the local scale — high values at
# woodland edges, heterogeneous scrub, and riparian margins. Low values in
# homogeneous open fields or uniform closed canopy.
# Neighbourhood: 3×3 cell rectangle (~3m × 3m at 1m resolution).
# Statistic: Standard Deviation.

print("\n--- STAGE 5: Vegetation Roughness Index (VRI = Focal STD 3×3) ---")
print(f"  Output: {VRI_1M}")

chm_raster = Raster(CHM_1M)

# 3×3 rectangular neighbourhood — 3 cells width and height
neighbourhood = NbrRectangle(3, 3, "CELL")

vri_raster = FocalStatistics(
    in_raster          = chm_raster,
    neighborhood       = neighbourhood,
    statistics_type    = "STD",
    ignore_nodata      = "DATA"
)

vri_raster.save(VRI_1M)

print("  VRI calculation complete.")

# Report basic statistics
vri_result = arcpy.management.GetRasterProperties(VRI_1M, "MINIMUM")
vri_min = float(vri_result.getOutput(0))
vri_result = arcpy.management.GetRasterProperties(VRI_1M, "MAXIMUM")
vri_max = float(vri_result.getOutput(0))
vri_result = arcpy.management.GetRasterProperties(VRI_1M, "MEAN")
vri_mean = float(vri_result.getOutput(0))

print(f"  VRI statistics — Min: {vri_min:.4f}  Max: {vri_max:.4f}  Mean: {vri_mean:.4f}")

print()
print("  VISUAL QA — VRI:")
print("  Expected: highest values at woodland edges and riparian margins.")
print("  Open improved grassland: near 0. Dense uniform canopy: low-moderate.")
print("  Hedgerows and field margins: moderate. Woodland edges: high.")
if vri_max < 0.5:
    print("  [WARNING] VRI max very low — check CHM values are correct first.")
else:
    print("  [OK] VRI range looks plausible.")

# Return Spatial Analyst extension
arcpy.CheckInExtension("Spatial")

# ---------------------------------------------------------------------------
# 6. FINAL SUMMARY AND QA CHECKLIST
# ---------------------------------------------------------------------------

print()
print("=" * 60)
print("SCRIPT 1 COMPLETE")
print(f"Finished: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 60)
print()
print("Outputs written to:")
print(f"  {OUT_DIR}")
print()
print("Files created:")
for fpath in [DTM_MOSAIC, DSM_MOSAIC, CHM_1M, VRI_1M]:
    exists = os.path.exists(fpath)
    status = "[OK]" if exists else "[MISSING]"
    print(f"  {status} {os.path.basename(fpath)}")

print()
print("VISUAL QA CHECKLIST — open each output in ArcGIS Pro and verify:")
print("  [ ] DTM_mosaic.tif — continuous elevation surface, no visible tile seams")
print("  [ ] DSM_mosaic.tif — slightly elevated over DTM where vegetation/buildings exist")
print("  [ ] CHM_1m.tif    — 0m fields, 0.5–5m hedgerows, 5–20m woodland, A66 = near 0")
print("  [ ] VRI_1m.tif    — high values at woodland edges and riparian margins")
print("  [ ] Overlay CHM on GAD georeferenced sheets — spot-check Kirkby Thore area")
print("  [ ] Eden Valley Railway corridor should show linear low-CHM feature")
print()
print("NEXT STEP: Run Script 2 — Raster Harmonisation")
print("  Reproject → Resample 2m → Clip to StudyArea_2km_buffer → Normalise 0–1")
print()
print("Push to GitHub at end of session:")
print("  https://github.com/Abhin-Nangari/A66-Connectivity-MSc")
