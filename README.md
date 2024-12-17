# Snakemake workflow: Tropical Cyclone and Storm Surge generation with CLIMADA

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)

## Usage

If you use this workflow in a paper please contact the author.

# Dataset Overview

This dataset provides historical and synthetic tropical cyclone (TC) tracks, wind-field calculations, and storm surge data under various climate scenarios and sea level rise (SLR) conditions. It is designed for studying the impacts of tropical cyclones and associated storm surges, with a focus on historical data and projections under climate change scenarios.

## Data Description

### Tropical Cyclone Tracks
   - Historical Tracks: Based on IBTrACS dataset IBTrACS.ALL.v04r00.nc for the period 1980–2023.
   - Synthetic Tracks: For each historical track, 25 synthetic tracks were generated, resulting in 26 members per cyclone.

### Wind-Field Computation
   - Computed using CLIMADA 5.0.0 and its petals (Holland 2008 model).
   - Resolution: 150 arcseconds (~4.5 km).

### Climate Change Scenarios
   - Climate scenarios applied using improved Knutson scaling, which updates event frequency:
        * RCP2.6 (low emissions scenario)
        * RCP8.5 (high emissions scenario)
    Time Horizon: 2060.

### Storm Surge Data
  - Computed using a bathtub model in CLIMADA petals.
  - Downscaling of wind-field to ~28 arcseconds (~0.008°) resolution via spatial interpolation.
  - Elevation data from SRTM15+V2 Digital Elevation Model (SRTM15+V2).
  - SLR scenarios are based on AR6 regional sea level rise data (AR6 SLR dataset):
      * No SLR
      * SLR under SSP1-2.6 (low SLR scenario) for 2060
      * SLR under SSP5-8.5 (high SLR scenario) for 2060

### File Naming Convention

Surge data files follow the format:

`surge/<Climate Scenario for TC>/<SLR Scenario>/surge_28arcsec_25synth_<country>_<SLR Scenario>_<Climate Scenario for TC>.hdf5`

Note that files with 0 size indicate no damaging events for the corresponding country.

### Tools and Dependencies
  * CLIMADA 5.0.0 and CLIMADA petals (GitHub Repository).
  * Elevation data from SRTM15+V2.
  * SLR data from AR6 SLR dataset.

# Contact

For any questions or further information, feel free to contact:
    Samuel Juhel
        Email: sjuhel[at]ethz.ch
        Email: pro[at]sjuhel.org

