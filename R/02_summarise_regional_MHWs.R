# ==============================================================================
# 02 - Summarise regional MHW metrics
# ==============================================================================
#
# Purpose:
#   Summarise marine heatwave events by aquaculture-relevant region and year
#   for 1990–2025.
#
# Inputs:
#   data/derived/regional_mhw_events/
#     all_regional_mhw_events_1990_2025.rds
#
#   data/derived/regional_mhw_events/
#     regional_valid_cell_counts.csv
#
# Outputs:
#   data/derived/regional_annual_mhw_metrics.csv
#   data/derived/regional_annual_mhw_metrics_long.csv
#   outputs/figures/regional_mhw_metrics.png
#
# Notes:
#   - Events are assigned to the year of peak intensity.
#   - Annual frequency is divided by the number of valid OSTIA ocean cells
#     in each region.
#   - Duration, mean intensity and cumulative intensity are averaged across
#     all events detected within each region and year.
# ==============================================================================


# ------------------------------------------------------------------------------
# Packages
# ------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)


# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------

events_file <- file.path(
  "data",
  "derived",
  "regional_mhw_events",
  "all_regional_mhw_events_1990_2025.rds"
)

cell_counts_file <- file.path(
  "data",
  "derived",
  "regional_mhw_events",
  "regional_valid_cell_counts.csv"
)

annual_metrics_file <- file.path(
  "data",
  "derived",
  "regional_annual_mhw_metrics.csv"
)

annual_metrics_long_file <- file.path(
  "data",
  "derived",
  "regional_annual_mhw_metrics_long.csv"
)

regional_plot_file <- file.path(
  "outputs",
  "figures",
  "regional_mhw_metrics.png"
)

dir.create(
  dirname(annual_metrics_file),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dirname(regional_plot_file),
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# Analysis settings
# ------------------------------------------------------------------------------

study_years <- 1990:2025

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
# Read input data
# ------------------------------------------------------------------------------

stopifnot(
  file.exists(events_file),
  file.exists(cell_counts_file)
)

regional_events <- readRDS(
  events_file
)

regional_cell_counts <- read.csv(
  cell_counts_file,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------------------------
# Check input structure
# ------------------------------------------------------------------------------

required_event_columns <- c(
  "region",
  "cell_id",
  "date_start",
  "date_peak",
  "date_end",
  "duration",
  "intensity_mean",
  "intensity_cumulative"
)

required_count_columns <- c(
  "region",
  "n_valid_cells"
)

stopifnot(
  all(
    required_event_columns %in%
      names(regional_events)
  ),
  all(
    required_count_columns %in%
      names(regional_cell_counts)
  ),
  !anyDuplicated(regional_cell_counts$region),
  setequal(
    regional_cell_counts$region,
    region_levels
  ),
  all(regional_cell_counts$n_valid_cells > 0)
)


# ------------------------------------------------------------------------------
# Prepare event data
# ------------------------------------------------------------------------------

regional_events <- regional_events %>%
  mutate(
    region = as.character(region),
    date_start = as.Date(date_start),
    date_peak = as.Date(date_peak),
    date_end = as.Date(date_end),
    year = as.integer(
      format(date_peak, "%Y")
    )
  )

regional_cell_counts <- regional_cell_counts %>%
  transmute(
    region = as.character(region),
    n_valid_cells = as.integer(n_valid_cells)
  )

stopifnot(
  setequal(
    unique(regional_events$region),
    region_levels
  ),
  !anyNA(regional_events$date_peak),
  all(regional_events$year %in% study_years),
  all(
    regional_events$date_peak >= regional_events$date_start &
      regional_events$date_peak <= regional_events$date_end
  ),
  all(regional_events$duration >= 5)
)


# ------------------------------------------------------------------------------
# Summarise events by region and year
# ------------------------------------------------------------------------------

annual_event_summary <- regional_events %>%
  group_by(
    region,
    year
  ) %>%
  summarise(
    n_events = n(),
    mean_duration = mean(
      duration,
      na.rm = TRUE
    ),
    mean_intensity = mean(
      intensity_mean,
      na.rm = TRUE
    ),
    mean_cumulative_intensity = mean(
      intensity_cumulative,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ------------------------------------------------------------------------------
# Construct the complete region-year table
# ------------------------------------------------------------------------------

regional_annual_metrics <- tidyr::expand_grid(
  region = region_levels,
  year = study_years
) %>%
  left_join(
    annual_event_summary,
    by = c("region", "year")
  ) %>%
  left_join(
    regional_cell_counts,
    by = "region"
  ) %>%
  mutate(
    n_events = replace_na(
      n_events,
      0L
    ),
    frequency_per_cell = (
      n_events / n_valid_cells
    ),
    region = factor(
      region,
      levels = region_levels
    )
  ) %>%
  select(
    region,
    year,
    n_valid_cells,
    n_events,
    frequency_per_cell,
    mean_duration,
    mean_intensity,
    mean_cumulative_intensity
  ) %>%
  arrange(
    region,
    year
  )


# ------------------------------------------------------------------------------
# Check annual summary
# ------------------------------------------------------------------------------

expected_rows <- (
  length(region_levels) *
    length(study_years)
)

stopifnot(
  nrow(regional_annual_metrics) == expected_rows,
  !anyDuplicated(
    regional_annual_metrics[
      c("region", "year")
    ]
  ),
  !anyNA(
    regional_annual_metrics$frequency_per_cell
  ),
  all(
    regional_annual_metrics$frequency_per_cell >= 0
  )
)


# ------------------------------------------------------------------------------
# Convert to long format for Figure 3
# ------------------------------------------------------------------------------

metric_labels <- c(
  frequency_per_cell =
    "Frequency\n(events pixel\u207B\u00B9 yr\u207B\u00B9)",
  mean_duration =
    "Annual mean duration\n(days)",
  mean_intensity =
    "Annual mean intensity\n(\u00B0C)",
  mean_cumulative_intensity =
    "Annual mean cumulative intensity\n(\u00B0C days)"
)

regional_annual_metrics_long <- regional_annual_metrics %>%
  select(
    region,
    year,
    frequency_per_cell,
    mean_duration,
    mean_intensity,
    mean_cumulative_intensity
  ) %>%
  pivot_longer(
    cols = c(
      frequency_per_cell,
      mean_duration,
      mean_intensity,
      mean_cumulative_intensity
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric_label = factor(
      metric,
      levels = names(metric_labels),
      labels = unname(metric_labels)
    )
  )


# ------------------------------------------------------------------------------
# Save annual summary tables
# ------------------------------------------------------------------------------

write.csv(
  regional_annual_metrics,
  annual_metrics_file,
  row.names = FALSE,
  na = ""
)

write.csv(
  regional_annual_metrics_long,
  annual_metrics_long_file,
  row.names = FALSE,
  na = ""
)


# ------------------------------------------------------------------------------
# Visualisation of regional MHW metrics
# ------------------------------------------------------------------------------

regional_metrics_plot <- ggplot(
  regional_annual_metrics_long,
  aes(
    x = value,
    y = region,
    colour = year
  )
) +
  geom_boxplot(
    aes(group = region),
    colour = "grey40",
    fill = "grey92",
    outlier.shape = NA,
    width = 0.6,
    na.rm = TRUE
  ) +
  geom_jitter(
    height = 0.12,
    width = 0,
    size = 1.2,
    alpha = 0.8,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~metric_label,
    scales = "free_x",
    ncol = 2
  ) +
  scale_y_discrete(
    limits = rev(region_levels),
    drop = FALSE
  ) +
  scale_colour_viridis_c(
    name = "Year",
    limits = range(study_years)
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
    legend.position = "right",
    strip.text = element_text(
      face = "bold"
    )
  )

print(regional_metrics_plot)

ggsave(
  filename = regional_plot_file,
  plot = regional_metrics_plot,
  width = 9,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)