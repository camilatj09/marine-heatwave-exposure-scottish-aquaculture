# ==============================================================================
# 09 - Create study-area map
# ==============================================================================
#
# Purpose:
#   Create a map showing the seven aquaculture-relevant marine regions,
#   aquaculture sites retained within the OSTIA ocean mask, and mean OSTIA sea
#   surface temperature during 2025.
#
# Inputs:
#   data/raw/OSTIA_1990_2025.nc
#
#   data/derived/spatial/
#     aquaculture_relevant_regions.gpkg
#
#   data/derived/aquaculture_sites/
#     aquaculture_sites_with_OSTIA.csv
#
# Output:
#   outputs/figures/
#     aquaculture_regions_sites_and_2025_SST.png
#
# Notes:
#   - Mean SST is calculated from the 2025 layers of the same reprocessed OSTIA
#     dataset used in the main analyses.
#   - The map is an illustrative example and does not reproduce the final
#     publication figure.
# ==============================================================================


# ------------------------------------------------------------------------------
# Packages
# ------------------------------------------------------------------------------

library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)


# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------

sst_file <- file.path(
  "data",
  "raw",
  "OSTIA_1990_2025.nc"
)

regions_file <- file.path(
  "data",
  "derived",
  "spatial",
  "aquaculture_relevant_regions.gpkg"
)

sites_file <- file.path(
  "data",
  "derived",
  "aquaculture_sites",
  "aquaculture_sites_with_OSTIA.csv"
)

figure_file <- file.path(
  "outputs",
  "figures",
  "aquaculture_regions_sites_and_2025_SST.png"
)

dir.create(
  dirname(figure_file),
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# Analysis settings
# ------------------------------------------------------------------------------

sst_variable <- "analysed_sst"

map_year <- 2025

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
# Read spatial inputs
# ------------------------------------------------------------------------------

stopifnot(
  file.exists(sst_file),
  file.exists(regions_file),
  file.exists(sites_file)
)

aquaculture_regions <- sf::st_read(
  regions_file,
  quiet = TRUE
)

aquaculture_sites <- read.csv(
  sites_file,
  stringsAsFactors = FALSE
)

stopifnot(
  all(c("region", "geometry") %in% names(aquaculture_regions)),
  all(
    c(
      "site_id",
      "longitude",
      "latitude"
    ) %in% names(aquaculture_sites)
  ),
  !anyNA(aquaculture_sites$longitude),
  !anyNA(aquaculture_sites$latitude)
)

aquaculture_regions <- aquaculture_regions %>%
  mutate(
    region = factor(
      as.character(region),
      levels = region_levels
    )
  ) %>%
  arrange(region)

aquaculture_sites_sf <- sf::st_as_sf(
  aquaculture_sites,
  coords = c(
    "longitude",
    "latitude"
  ),
  crs = 4326,
  remove = FALSE
)


# ------------------------------------------------------------------------------
# Read and summarise 2025 OSTIA SST
# ------------------------------------------------------------------------------

sst <- terra::rast(
  sst_file,
  subds = sst_variable,
  drivers = "NETCDF"
)

sst_dates <- as.Date(
  terra::time(sst)
)

stopifnot(
  length(sst_dates) == terra::nlyr(sst),
  !anyNA(sst_dates)
)

keep_2025 <- as.integer(
  format(sst_dates, "%Y")
) == map_year

stopifnot(
  any(keep_2025)
)

sst_2025 <- sst[[
  which(keep_2025)
]]

# The original analysed_sst variable is supplied in kelvin.
sst_2025 <- sst_2025 - 273.15

mean_sst_2025 <- terra::app(
  sst_2025,
  mean,
  na.rm = TRUE
)

names(mean_sst_2025) <- "mean_sst"


# ------------------------------------------------------------------------------
# Crop SST to the study area
# ------------------------------------------------------------------------------

regions_vect <- terra::vect(
  aquaculture_regions
)

regions_vect <- terra::project(
  regions_vect,
  terra::crs(mean_sst_2025)
)

study_extent <- terra::ext(
  regions_vect
)

# Add a small margin around the seven regions.
study_extent <- study_extent + 0.4

mean_sst_crop <- terra::crop(
  mean_sst_2025,
  study_extent
)

sst_plot_data <- as.data.frame(
  mean_sst_crop,
  xy = TRUE,
  na.rm = TRUE
)


# ------------------------------------------------------------------------------
# Prepare map layers
# ------------------------------------------------------------------------------

map_crs <- sf::st_crs(
  terra::crs(mean_sst_crop)
)

aquaculture_regions <- sf::st_transform(
  aquaculture_regions,
  map_crs
)

aquaculture_sites_sf <- sf::st_transform(
  aquaculture_sites_sf,
  map_crs
)

land <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
) %>%
  sf::st_transform(
    map_crs
  )


# ------------------------------------------------------------------------------
# Simple study-area visualisation
# ------------------------------------------------------------------------------

study_area_map <- ggplot() +
  geom_raster(
    data = sst_plot_data,
    aes(
      x = x,
      y = y,
      fill = mean_sst
    )
  ) +
  geom_sf(
    data = land,
    fill = "grey90",
    colour = "grey60",
    linewidth = 0.2
  ) +
  geom_sf(
    data = aquaculture_regions,
    fill = NA,
    colour = "black",
    linewidth = 0.5
  ) +
  geom_sf(
    data = aquaculture_sites_sf,
    colour = "black",
    size = 0.8,
    alpha = 0.7
  ) +
  scale_fill_viridis_c(
    name = "Mean SST\n2025 (°C)"
  ) +
  coord_sf(
    xlim = c(
      terra::xmin(study_extent),
      terra::xmax(study_extent)
    ),
    ylim = c(
      terra::ymin(study_extent),
      terra::ymax(study_extent)
    ),
    expand = FALSE
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(study_area_map)

ggsave(
  filename = figure_file,
  plot = study_area_map,
  width = 8,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)