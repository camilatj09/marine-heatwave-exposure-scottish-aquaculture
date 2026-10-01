# ==============================================================================
# 04 - Extract OSTIA SST and detect MHWs at unique aquaculture grid cells
# ==============================================================================
# Purpose:
#   Extract daily OSTIA SST at active Scottish aquaculture locations, identify
#   farms sharing the same OSTIA cell, retain one SST history per unique cell,
#   and detect marine heatwaves independently for those cells.
#
# Analysis units:
#   361 farm locations represented by 218 independent OSTIA grid cells.
#
# Inputs:
#   data/raw/OSTIA_1990_2025.nc
#   data/derived/aquaculture_sites/
#     scotland_active_marine_aquaculture_sites_unique.rds
#
# Main outputs:
#   data/derived/aquaculture_sites/aquaculture_sites_with_OSTIA.csv
#   data/derived/aquaculture_cells/scotland_farm_OSTIA_cell_lookup.rds
#   data/derived/aquaculture_cells/OSTIA_unique_cell_SST_1990_2025.rds
#   data/derived/aquaculture_cells/OSTIA_cell_representation_summary.csv
#   data/derived/cell_mhw_events/
#     scotland_OSTIA_cell_mhw_events_1990_2025.rds
# ==============================================================================

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(heatwaveR)

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------

sst_file <- file.path("data", "raw", "OSTIA_1990_2025.nc")

sites_file <- file.path(
  "data", "derived", "aquaculture_sites",
  "scotland_active_marine_aquaculture_sites_unique.rds"
)

sites_outside_mask_file <- file.path(
  "data", "derived", "aquaculture_sites",
  "aquaculture_sites_outside_OSTIA_mask.csv"
)

sites_with_ostia_file <- file.path(
  "data", "derived", "aquaculture_sites",
  "aquaculture_sites_with_OSTIA.csv"
)

cell_dir <- file.path("data", "derived", "aquaculture_cells")
events_dir <- file.path("data", "derived", "cell_mhw_events")

farm_cell_lookup_file <- file.path(
  cell_dir, "scotland_farm_OSTIA_cell_lookup.rds"
)

cell_sst_file <- file.path(
  cell_dir, "OSTIA_unique_cell_SST_1990_2025.rds"
)

cell_representation_file <- file.path(
  cell_dir, "OSTIA_cell_representation_summary.csv"
)

cell_counts_file <- file.path(
  cell_dir, "OSTIA_cell_farm_counts.csv"
)

events_file <- file.path(
  events_dir, "scotland_OSTIA_cell_mhw_events_1990_2025.rds"
)

dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(events_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Analysis settings
# ------------------------------------------------------------------------------

sst_variable <- "analysed_sst"

study_start <- as.Date("1990-01-01")
study_end   <- as.Date("2025-12-31")

clim_start <- "1990-01-01"
clim_end   <- "2020-12-31"

window_half_width       <- 5
threshold_percentile    <- 90
smooth_percentile_width <- 31
minimum_duration        <- 5
maximum_gap             <- 2

expected_n_farms  <- 361L
expected_n_cells  <- 218L
expected_n_events <- 21099L

# ------------------------------------------------------------------------------
# Read OSTIA SST
# ------------------------------------------------------------------------------

stopifnot(file.exists(sst_file), file.exists(sites_file))

sst <- terra::rast(
  sst_file,
  subds = sst_variable,
  drivers = "NETCDF"
)

# analysed_sst is supplied in kelvin.
sst <- sst - 273.15

sst_dates <- as.Date(terra::time(sst))
keep <- sst_dates >= study_start & sst_dates <= study_end

sst <- sst[[which(keep)]]
sst_dates <- sst_dates[keep]

expected_dates <- seq.Date(study_start, study_end, by = "day")

stopifnot(
  length(sst_dates) == length(expected_dates),
  all(sst_dates == expected_dates),
  !anyDuplicated(sst_dates)
)

names(sst) <- format(sst_dates, "%Y-%m-%d")
sst_columns <- names(sst)

# ------------------------------------------------------------------------------
# Read aquaculture locations and assign OSTIA cells
# ------------------------------------------------------------------------------

aquaculture_sites <- readRDS(sites_file)

required_site_columns <- c(
  "site_id", "source_FID", "siteno", "sitename", "sitetype", "speciescd",
  "n_species", "mgtarea", "regionname", "longitude", "latitude"
)

stopifnot(
  all(required_site_columns %in% names(aquaculture_sites)),
  !anyDuplicated(aquaculture_sites$site_id)
)

sites_sf <- sf::st_as_sf(
  aquaculture_sites,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

sites_vect <- terra::vect(sites_sf)

if (!terra::same.crs(sites_vect, sst)) {
  sites_vect <- terra::project(sites_vect, terra::crs(sst))
}

site_xy <- terra::crds(sites_vect)
site_cells <- terra::cellFromXY(sst[[1]], site_xy)

# ------------------------------------------------------------------------------
# Extract daily SST at farm locations and remove cells outside the ocean mask
# ------------------------------------------------------------------------------

sst_extracted <- terra::extract(sst, sites_vect, ID = TRUE)

stopifnot(
  nrow(sst_extracted) == nrow(aquaculture_sites),
  all(sst_extracted$ID == seq_len(nrow(aquaculture_sites)))
)

site_metadata <- aquaculture_sites %>%
  select(all_of(required_site_columns)) %>%
  mutate(ostia_cell = site_cells)

sst_extracted <- sst_extracted %>%
  mutate(site_id = aquaculture_sites$site_id[ID], .after = ID) %>%
  left_join(site_metadata, by = "site_id")

all_na <- rowSums(
  !is.na(as.matrix(sst_extracted[, sst_columns, drop = FALSE]))
) == 0

sites_outside_mask <- sst_extracted[all_na, ]
sites_with_sst <- sst_extracted[!all_na, ]

sites_outside_mask %>%
  select(any_of(c(required_site_columns, "ostia_cell"))) %>%
  write.csv(sites_outside_mask_file, row.names = FALSE, na = "")

# ------------------------------------------------------------------------------
# Create farm-to-cell lookup and unique-cell metadata
# ------------------------------------------------------------------------------

cell_counts <- sites_with_sst %>%
  count(ostia_cell, name = "n_farms_in_cell") %>%
  arrange(desc(n_farms_in_cell), ostia_cell)

unique_cells <- sort(unique(sites_with_sst$ostia_cell))
cell_xy <- terra::xyFromCell(sst[[1]], unique_cells)

cell_points <- terra::vect(
  data.frame(x = cell_xy[, 1], y = cell_xy[, 2]),
  geom = c("x", "y"),
  crs = terra::crs(sst)
)

cell_points_wgs84 <- terra::project(cell_points, "EPSG:4326")
cell_lonlat <- terra::crds(cell_points_wgs84)

cell_metadata <- data.frame(
  ostia_cell = unique_cells,
  ostia_cell_lon = cell_lonlat[, 1],
  ostia_cell_lat = cell_lonlat[, 2]
) %>%
  left_join(cell_counts, by = "ostia_cell")

farm_cell_lookup <- sites_with_sst %>%
  select(any_of(required_site_columns), ostia_cell) %>%
  left_join(cell_metadata, by = "ostia_cell", relationship = "many-to-one") %>%
  arrange(site_id)

sites_with_ostia <- farm_cell_lookup

write.csv(
  sites_with_ostia,
  sites_with_ostia_file,
  row.names = FALSE,
  na = ""
)

cell_representation_summary <- data.frame(
  metric = c(
    "Farm locations with valid OSTIA SST",
    "Unique OSTIA grid cells",
    "OSTIA cells containing more than one farm",
    "Farms located in shared OSTIA cells"
  ),
  value = c(
    nrow(farm_cell_lookup),
    nrow(cell_metadata),
    sum(cell_counts$n_farms_in_cell > 1),
    sum(cell_counts$n_farms_in_cell[cell_counts$n_farms_in_cell > 1])
  )
)

stopifnot(
  nrow(farm_cell_lookup) == expected_n_farms,
  nrow(cell_metadata) == expected_n_cells,
  n_distinct(farm_cell_lookup$ostia_cell) == expected_n_cells
)

saveRDS(farm_cell_lookup, farm_cell_lookup_file)
write.csv(cell_representation_summary, cell_representation_file, row.names = FALSE)
write.csv(cell_counts, cell_counts_file, row.names = FALSE)

# ------------------------------------------------------------------------------
# Retain one SST series per unique OSTIA cell
# ------------------------------------------------------------------------------

cell_sst_wide <- sites_with_sst %>%
  arrange(ostia_cell, site_id) %>%
  distinct(ostia_cell, .keep_all = TRUE) %>%
  select(ostia_cell, all_of(sst_columns))

cell_sst <- cell_sst_wide %>%
  pivot_longer(
    cols = all_of(sst_columns),
    names_to = "date",
    values_to = "sst"
  ) %>%
  mutate(date = as.Date(date)) %>%
  left_join(cell_metadata, by = "ostia_cell", relationship = "many-to-one") %>%
  arrange(ostia_cell, date)

stopifnot(
  n_distinct(cell_sst$ostia_cell) == expected_n_cells,
  nrow(cell_sst) == expected_n_cells * length(expected_dates),
  !anyDuplicated(cell_sst[c("ostia_cell", "date")]),
  !anyNA(cell_sst$sst)
)

saveRDS(cell_sst, cell_sst_file)

# ------------------------------------------------------------------------------
# Marine heatwave detection
# ------------------------------------------------------------------------------

run_mhw <- function(cell_data) {

  ts <- cell_data %>%
    transmute(t = date, temp = sst) %>%
    arrange(t)

  climatology <- heatwaveR::ts2clm(
    data = ts,
    x = t,
    y = temp,
    climatologyPeriod = c(clim_start, clim_end),
    windowHalfWidth = window_half_width,
    pctile = threshold_percentile,
    smoothPercentile = TRUE,
    smoothPercentileWidth = smooth_percentile_width,
    returnDF = TRUE
  )

  detected <- heatwaveR::detect_event(
    data = climatology,
    x = t,
    y = temp,
    minDuration = minimum_duration,
    joinAcrossGaps = TRUE,
    maxGap = maximum_gap,
    coldSpells = FALSE,
    categories = TRUE,
    returnDF = TRUE
  )

  if (is.list(detected) &&
      !inherits(detected, "data.frame") &&
      "event" %in% names(detected)) {
    detected <- detected$event
  }

  detected
}

cell_mhw_events <- cell_sst %>%
  select(ostia_cell, date, sst) %>%
  group_by(ostia_cell) %>%
  nest() %>%
  mutate(events = map(data, run_mhw)) %>%
  select(ostia_cell, events) %>%
  unnest(events) %>%
  left_join(cell_metadata, by = "ostia_cell", relationship = "many-to-one") %>%
  arrange(ostia_cell, date_start)

required_event_columns <- c(
  "ostia_cell", "date_start", "date_peak", "date_end", "duration",
  "intensity_mean", "intensity_max", "intensity_cumulative", "category"
)

stopifnot(
  all(required_event_columns %in% names(cell_mhw_events)),
  n_distinct(cell_mhw_events$ostia_cell) == expected_n_cells,
  nrow(cell_mhw_events) == expected_n_events,
  all(cell_mhw_events$duration >= minimum_duration)
)

saveRDS(cell_mhw_events, events_file)

print(cell_representation_summary)
cat("Detected MHW events:", nrow(cell_mhw_events), "\n")
