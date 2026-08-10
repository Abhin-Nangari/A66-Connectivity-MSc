# A66-Connectivity-MSc
Automated ArcPy pipeline for wildlife connectivity modelling along the A66 Northern Trans-Pennine Corridor — MSc GIS Dissertation, University of Aberdeen 2026

# A66 Northern Trans-Pennine Corridor — Wildlife Connectivity Modelling

## MSc GIS Dissertation — University of Aberdeen 2026
**Student:** Abhin Nangari  
**Supervisor:** Dr Shaktiman Singh  
**Programme:** MSc Geographical Information Systems GG5910/GG5912  
**Deadline:** 25 August 2026  

## Study Area
DCO Scheme 0405 — Temple Sowerby to Appleby bypass junction, Westmorland, Cumbria  
8.2km corridor confirmed from SI 2024/360 Schedule 7 Part 3  

## Dissertation Title
Automating Wildlife Connectivity Modelling Along the A66 Northern Trans-Pennine Corridor: An ArcPy-Based GIS Workflow Integrating Environment Agency LiDAR Structural Metrics, Fuzzy Resistance Modelling, and Circuit Theory

## Pipeline Overview
1. EA LiDAR preprocessing — CHM and VRI derivation
2. Raster harmonisation — reproject, resample, clip, normalise
3. Fuzzy MCE resistance surface — 4 sensitivity scenarios
4. Circuitscape 5.0 pairwise connectivity modelling
5. 3-layer structural validation
6. GitHub toolbox — Sphinx documented, redeployable

## Data Sources
- EA LiDAR Composite DTM and DSM 1m (2022)
- UKCEH Land Cover Map 2023 10m
- Natural England Living England 2022-23 10m
- OS Terrain 50 DTM
- OS Open Roads and OS Open Rivers
- National Highways TR010062 DCO documents
- Natural England Habitat Networks England
- Cumbria LNRS spatial layers

## Focal Species
Badger (Meles meles), Roe deer (Capreolus capreolus), Pine marten (Martes martes)
All confirmed present by Arup PEI 2021 (TR010062)

## Repository Structure
- 03_Scripts/ — all ArcPy Python scripts
- README.md — project documentation
