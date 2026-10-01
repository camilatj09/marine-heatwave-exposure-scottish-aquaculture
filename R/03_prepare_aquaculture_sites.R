# ==============================================================================
# 03 - Prepare Scottish aquaculture site locations
# ==============================================================================
# Purpose:
#   Select active marine finfish and shellfish records, combine duplicate source
#   records at identical coordinates, and create one geographic farm location
#   per site for subsequent OSTIA extraction.
#
# Input:
#   data/raw/scottish_aquaculture_sites/
#     marine_freshwater_aquaculture_sites.shp
#
# Outputs:
#   data/derived/aquaculture_sites/
#     scotland_active_marine_aquaculture_records.rds
#     scotland_active_marine_aquaculture_sites_unique.rds
# ==============================================================================

library(sf)
library(dplyr)
library(stringr)

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------

aquaculture_file <- file.path(
  "data", "raw", "scottish_aquaculture_sites",
  "marine_freshwater_aquaculture_sites.shp"
)

output_dir <- file.path("data", "derived", "aquaculture_sites")

active_records_file <- file.path(
  output_dir, "scotland_active_marine_aquaculture_records.rds"
)

unique_sites_file <- file.path(
  output_dir, "scotland_active_marine_aquaculture_sites_unique.rds"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Read and filter source records
# ------------------------------------------------------------------------------

stopifnot(file.exists(aquaculture_file))

aquaculture_all <- sf::st_read(aquaculture_file, quiet = TRUE)

required_columns <- c(
  "FID", "sitename", "siteno", "sitestatcd", "watertype", "sitetype",
  "easting", "northing", "speciescd", "stagecode"
)

stopifnot(all(required_columns %in% names(aquaculture_all)))

aquaculture_marine_active <- aquaculture_all %>%
  filter(
    watertype %in% c("S", "B"),
    sitestatcd == "A",
    sitetype %in% c("Fish", "Shellfish")
  )

if (inherits(aquaculture_marine_active, "sf")) {
  aquaculture_marine_active <- sf::st_drop_geometry(aquaculture_marine_active)
}

aquaculture_marine_active <- aquaculture_marine_active %>%
  mutate(
    across(where(is.character), stringr::str_trim),
    across(where(is.character), ~ na_if(.x, ""))
  )

stopifnot(
  nrow(aquaculture_marine_active) > 0,
  !anyNA(aquaculture_marine_active$easting),
  !anyNA(aquaculture_marine_active$northing)
)

# ------------------------------------------------------------------------------
# Add WGS84 coordinates to the selected source records
# ------------------------------------------------------------------------------

active_sf <- sf::st_as_sf(
  aquaculture_marine_active,
  coords = c("easting", "northing"),
  crs = 27700,
  remove = FALSE
) %>%
  sf::st_transform(4326)

active_xy <- sf::st_coordinates(active_sf)

aquaculture_marine_active <- active_sf %>%
  mutate(
    longitude = active_xy[, 1],
    latitude = active_xy[, 2]
  ) %>%
  sf::st_drop_geometry() %>%
  select(-any_of("scotaqurl"))

saveRDS(aquaculture_marine_active, active_records_file)

# ------------------------------------------------------------------------------
# Combine records at identical coordinates
# ------------------------------------------------------------------------------

collapse_values <- function(x) {
  x <- sort(unique(na.omit(x)))
  paste(x, collapse = ", ")
}

aquaculture_sites_unique <- aquaculture_marine_active %>%
  group_by(easting, northing) %>%
  summarise(
    source_FID = first(FID),
    sitename = first(sitename),
    siteno = first(siteno),
    n_records = n(),
    n_species = n_distinct(speciescd[!is.na(speciescd)]),
    n_stages = n_distinct(stagecode[!is.na(stagecode)]),
    n_types = n_distinct(sitetype[!is.na(sitetype)]),
    sitecatcd = collapse_values(sitecatcd),
    sitestatcd = collapse_values(sitestatcd),
    sitetype = collapse_values(sitetype),
    mgtarea = collapse_values(mgtarea),
    unitauth = collapse_values(unitauth),
    regionname = collapse_values(regionname),
    watertype = collapse_values(watertype),
    osgridref = collapse_values(osgridref),
    facilitytp = collapse_values(facilitytp),
    stagecode = collapse_values(stagecode),
    speciescd = collapse_values(speciescd),
    .groups = "drop"
  )

unique_sf <- sf::st_as_sf(
  aquaculture_sites_unique,
  coords = c("easting", "northing"),
  crs = 27700,
  remove = FALSE
) %>%
  sf::st_transform(4326)

unique_xy <- sf::st_coordinates(unique_sf)

aquaculture_sites_unique <- unique_sf %>%
  mutate(
    longitude = unique_xy[, 1],
    latitude = unique_xy[, 2]
  ) %>%
  sf::st_drop_geometry() %>%
  arrange(northing, easting) %>%
  mutate(site_id = row_number(), .before = 1)

# ------------------------------------------------------------------------------
# Checks and output
# ------------------------------------------------------------------------------

stopifnot(
  !anyDuplicated(aquaculture_sites_unique$site_id),
  !anyDuplicated(aquaculture_sites_unique[c("easting", "northing")]),
  !anyNA(aquaculture_sites_unique$longitude),
  !anyNA(aquaculture_sites_unique$latitude),
  sum(aquaculture_sites_unique$n_records) == nrow(aquaculture_marine_active)
)

saveRDS(aquaculture_sites_unique, unique_sites_file)
