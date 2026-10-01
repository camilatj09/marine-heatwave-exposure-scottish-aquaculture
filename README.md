# Analysis workflow

This repository contains the R code supporting the manuscript 
**“Marine heatwave exposure is changing unevenly across Scottish aquaculture sites.”**

## Script order

Scripts should be run from the repository root in numerical order.

1. `R/01_prepare_and_detect_regional_MHWs.R`  
   Defines the seven aquaculture-relevant regions and detects MHWs across 
   OSTIA grid cells for the regional analysis.

2. `R/02_summarise_regional_MHWs.R`  
   Calculates annual regional MHW metrics and produces the regional time-series figure.

3. `R/03_prepare_aquaculture_sites.R`  
   Prepares active marine finfish and shellfish locations for OSTIA extraction.

4. `R/04_extract_and_detect_cell_MHWs.R`  
   Extracts OSTIA SST at farm locations, identifies the 218 unique OSTIA cells 
   represented by the 361 farms, and detects MHWs at cell level.

5. `R/05_PCA_and_clustering.R`  
   Calculates long-term MHW metrics, performs PCA and hierarchical clustering, 
   and assigns the resulting exposure profiles to aquaculture locations.

6. `R/06_continuous_year_models.R`  
   Fits the continuous-year models for MHW frequency, duration, mean intensity, 
   cumulative intensity and severity, including model diagnostics and the 
   negative-binomial sensitivity check for frequency.

7. `R/07_period_models_and_breakpoint_sensitivity.R`  
   Fits the pre-/post-2020 models and repeats the analysis using 2017, 2018 and 
   2019 as alternative temporal boundaries.

8. `R/08_in_situ_QC_and_validation.R`  
   Implements the quality-control and OSTIA–in situ validation workflow. 
   The original farm-temperature observations are confidential and are not 
   included in the repository; when they are unavailable, the script uses 
   synthetic data to demonstrate the workflow.

9. `R/09_create_study_area_map.R`  
   Produces the study-area map using aquaculture locations, 
   aquaculture-relevant regions and 2025 OSTIA SST.

`sessionInfo.txt` records the R environment used for the analysis.

## Required public inputs

- `data/raw/OSTIA_1990_2025.nc`
- `data/raw/scottish_marine_regions/administrative_units_scottish_marine_regions.shp`
- `data/raw/scottish_aquaculture_sites/marine_freshwater_aquaculture_sites.shp`

Data sources and access information are provided in the manuscript Data Availability statement.

## In situ validation

The farm-temperature observations used for validation are not publicly available. 
For an authorised local run, `R/08_in_situ_QC_and_validation.R` can use an anonymised 
file at `data/confidential/in_situ_OSTIA_joined.csv` containing `site_id`, `date`, 
`site_temperature` and `OSTIA_temperature`. An optional 
`data/confidential/in_situ_manual_exclusions.csv` file can be used to identify known 
erroneous observations by `site_id` and `date`.

The cell-level workflow includes checks for the analytical dataset used in the study: 
361 aquaculture locations, 218 unique OSTIA cells and 21,099 detected MHW events.
