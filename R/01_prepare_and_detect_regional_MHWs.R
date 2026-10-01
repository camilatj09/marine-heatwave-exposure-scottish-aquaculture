# ==============================================================================
# 01 - Prepare regions and detect regional MHWs
# ==============================================================================
#
# Purpose:
#   Construct the seven aquaculture-relevant marine regions used in the
#   analysis and detect marine heatwaves independently in each OSTIA ocean
#   grid cell from 1990 to 2025.
#
# Data sources:
#
#   OSTIA reprocessed Level 4 daily sea surface temperature
#   Copernicus Marine Service
#   Product: SST_GLO_SST_L4_REP_OBSERVATIONS_010_011
#   Dataset: METOFFICE-GLO-SST-L4-REP-OBS-SST
#   Variable: analysed_sst
#   Spatial resolution: 0.05°
#   Spatial extent: -9.5, 2.7, 48.85, 61.2  (xmin, xmax, ymin, ymax)
#   DOI: https://doi.org/10.48670/moi-00168
#  
#   Administrative Units - Scottish Marine Regions
#   Scottish Government Marine Directorate
#   Metadata:
#   https://spatialdata.gov.scot/geonetwork/srv/eng/catalog.search#/
#   metadata/Marine_Scotland_FishDAC_1609
#
#   Contains information from Scottish Government licensed under the
#   Open Government Licence v3.0.
#
# Inputs:
#   data/raw/OSTIA_1990_2025.nc
#   data/raw/scottish_marine_regions/
#     administrative_units_scottish_marine_regions.shp
#
# Outputs:
#   data/derived/regional_mhw_events/
#     all_regional_mhw_events_1990_2025.rds
#
#   data/derived/regional_mhw_events/
#     regional_valid_cell_counts.csv 
#
#    data/derived/spatial/aquaculture_relevant_regions.gpkg
#
# Notes:
#   - OSTIA SST is converted from kelvin to degrees Celsius.
#   - A fixed 1990-01-01 to 2020-12-31 climatology is used.
# ==============================================================================


# ------------------------------------------------------------------------------
# Packages
# ------------------------------------------------------------------------------

library(terra)
library(sf)
library(dplyr)
library(heatwaveR)


# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------

sst_file <- file.path(
  "data",
  "raw",
  "OSTIA_1990_2025.nc"
)

smr_file <- file.path(
  "data",
  "raw",
  "scottish_marine_regions",
  "administrative_units_scottish_marine_regions.shp"
)

output_dir <- file.path(
  "data",
  "derived",
  "regional_mhw_events"
)

events_file <- file.path(
  output_dir,
  "all_regional_mhw_events_1990_2025.rds"
)

cell_counts_file <- file.path(
  output_dir,
  "regional_valid_cell_counts.csv"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

regions_file <- file.path(
  "data",
  "derived",
  "spatial",
  "aquaculture_relevant_regions.gpkg"
)

dir.create(
  dirname(regions_file),
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# Analysis settings
# ------------------------------------------------------------------------------

sst_variable <- "analysed_sst"

study_start <- as.Date("1990-01-01")
study_end   <- as.Date("2025-12-31")

clim_start <- as.Date("1990-01-01")
clim_end   <- as.Date("2020-12-31")

clim_window_half_width <- 5
clim_smoothing_width   <- 31
mhw_percentile         <- 90
min_duration           <- 5
max_gap                <- 2

region_levels <- c(
  "Shetland Isles",
  "Orkney Islands",
  "North Coast & West Highlands",
  "Outer Hebrides",
  "Argyll & Clyde",
  "Moray Firth & North East",
  "Forth & Tay"
)


# ------------------------------------------------------------------------------
# Read OSTIA SST
# ------------------------------------------------------------------------------

stopifnot(
  file.exists(sst_file),
  file.exists(smr_file)
)

# Read the analysed SST variable from the Copernicus NetCDF file.
sst <- terra::rast(
  sst_file,
  subds = sst_variable,
  drivers = "NETCDF"
)

# The original analysed_sst variable is supplied in kelvin.
sst <- sst - 273.15

sst_dates <- as.Date(terra::time(sst))

stopifnot(
  length(sst_dates) == terra::nlyr(sst),
  !anyNA(sst_dates),
  !anyDuplicated(sst_dates)
)

# Retain the exact study period.
keep_layers <- (
  sst_dates >= study_start &
    sst_dates <= study_end
)

sst <- sst[[which(keep_layers)]]
sst_dates <- sst_dates[keep_layers]

expected_dates <- seq.Date(
  study_start,
  study_end,
  by = "day"
)

stopifnot(
  length(sst_dates) == length(expected_dates),
  all(sst_dates == expected_dates)
)


# ------------------------------------------------------------------------------
# Construct aquaculture-relevant regions
# ------------------------------------------------------------------------------

smr <- sf::st_read(
  smr_file,
  quiet = TRUE
)

stopifnot(
  "objnam" %in% names(smr)
)

# Standardise an alternative spelling that may occur in the source file.
smr <- smr %>%
  mutate(
    objnam = recode(
      objnam,
      "Forth and Tay" = "Forth & Tay"
    )
  )

arr_definition <- tibble::tribble(
  ~source_region,    ~region,
  "Shetland Isles",  "Shetland Isles",
  "Orkney Islands",  "Orkney Islands",
  "North Coast",     "North Coast & West Highlands",
  "West Highlands",  "North Coast & West Highlands",
  "Outer Hebrides",  "Outer Hebrides",
  "Argyll",          "Argyll & Clyde",
  "Clyde",           "Argyll & Clyde",
  "Moray Firth",     "Moray Firth & North East",
  "North East",      "Moray Firth & North East",
  "Forth & Tay",     "Forth & Tay"
)

stopifnot(
  all(arr_definition$source_region %in% smr$objnam)
)

# Attach the regional grouping and dissolve boundaries between source regions.
arr <- smr %>%
  filter(
    objnam %in% arr_definition$source_region
  ) %>%
  left_join(
    arr_definition,
    by = c("objnam" = "source_region")
  ) %>%
  group_by(region) %>%
  summarise(
    .groups = "drop"
  ) %>%
  mutate(
    region = factor(
      region,
      levels = region_levels
    )
  ) %>%
  arrange(region)

stopifnot(
  nrow(arr) == length(region_levels),
  setequal(as.character(arr$region), region_levels)
)

# Match the polygons to the OSTIA coordinate reference system.
arr <- sf::st_transform(
  arr,
  crs = terra::crs(sst)
)

arr_vect <- terra::vect(arr)

# Export
sf::st_write(
  arr,
  regions_file,
  delete_dsn = TRUE,
  quiet = TRUE
)

# ------------------------------------------------------------------------------
# Detect MHW events in one SST time series
# ------------------------------------------------------------------------------

detect_pixel_events <- function(temp, dates) {
  
  pixel_data <- data.frame(
    t = dates,
    temp = as.numeric(temp)
  )
  
  climatology <- heatwaveR::ts2clm(
    pixel_data,
    climatologyPeriod = format(
      c(clim_start, clim_end),
      "%Y-%m-%d"
    ),
    maxPadLength = FALSE,
    windowHalfWidth = clim_window_half_width,
    pctile = mhw_percentile,
    smoothPercentile = TRUE,
    smoothPercentileWidth = clim_smoothing_width,
    roundClm = 4
  )
  
  heatwaveR::detect_event(
    climatology,
    minDuration = min_duration,
    joinAcrossGaps = TRUE,
    maxGap = max_gap,
    categories = FALSE,
    roundRes = 4
  )$event
}


# ------------------------------------------------------------------------------
# Detect MHW events in one region
# ------------------------------------------------------------------------------

process_region <- function(region_name) {
  
  region_polygon <- arr_vect[
    arr_vect$region == region_name,
  ]
  
  # Crop to the regional extent and mask cells outside the polygon.
  region_sst <- terra::crop(
    sst,
    region_polygon,
    mask = TRUE,
    touches = TRUE
  )
  
  # OSTIA is gap-free over ocean cells, so the first daily layer is sufficient
  # to identify the valid ocean-cell mask.
  valid_local_cells <- which(
    !is.na(
      terra::values(
        region_sst[[1]],
        mat = FALSE
      )
    )
  )
  
  stopifnot(
    length(valid_local_cells) > 0
  )
  
  cell_coordinates <- terra::xyFromCell(
    region_sst,
    valid_local_cells
  )
  
  # Record cell identifiers relative to the original OSTIA raster rather than
  # the cropped regional raster.
  global_cell_ids <- terra::cellFromXY(
    sst[[1]],
    cell_coordinates
  )
  
  temperature_matrix <- terra::values(
    region_sst,
    mat = TRUE
  )[valid_local_cells, , drop = FALSE]
  
  stopifnot(
    nrow(temperature_matrix) == length(valid_local_cells),
    ncol(temperature_matrix) == length(sst_dates),
    !anyNA(global_cell_ids)
  )
  
  event_list <- lapply(
    seq_len(nrow(temperature_matrix)),
    function(i) {
      
      events <- detect_pixel_events(
        temp = temperature_matrix[i, ],
        dates = sst_dates
      )
      
      if (nrow(events) == 0) {
        return(NULL)
      }
      
      events %>%
        mutate(
          region = region_name,
          cell_id = global_cell_ids[i],
          x = cell_coordinates[i, 1],
          y = cell_coordinates[i, 2],
          .before = 1
        )
    }
  )
  
  list(
    events = bind_rows(event_list),
    cell_count = tibble::tibble(
      region = region_name,
      n_valid_cells = length(valid_local_cells)
    )
  )
}


# ------------------------------------------------------------------------------
# Run detection across all seven regions
# ------------------------------------------------------------------------------

regional_results <- lapply(
  region_levels,
  process_region
)

regional_events <- bind_rows(
  lapply(
    regional_results,
    `[[`,
    "events"
  )
)

regional_cell_counts <- bind_rows(
  lapply(
    regional_results,
    `[[`,
    "cell_count"
  )
)


# ------------------------------------------------------------------------------
# Final checks
# ------------------------------------------------------------------------------

stopifnot(
  setequal(
    unique(regional_events$region),
    region_levels
  ),
  nrow(regional_cell_counts) == length(region_levels),
  !anyDuplicated(regional_cell_counts$region),
  all(regional_cell_counts$n_valid_cells > 0),
  all(regional_events$duration >= min_duration),
  all(
    regional_events$date_peak >= regional_events$date_start &
      regional_events$date_peak <= regional_events$date_end
  )
)

regional_cell_counts <- regional_cell_counts %>%
  mutate(
    region = factor(
      region,
      levels = region_levels
    )
  ) %>%
  arrange(region) %>%
  mutate(
    region = as.character(region)
  )


# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------

saveRDS(
  regional_events,
  events_file
)

write.csv(
  regional_cell_counts,
  cell_counts_file,
  row.names = FALSE
)