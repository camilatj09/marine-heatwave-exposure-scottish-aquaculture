# ==============================================================================
# 08 - In situ temperature QC and OSTIA validation workflow
# ==============================================================================
# Purpose:
#   QC and validation workflow applied to farm temperature observations. 
#
# Real-data input (not distributed):
#   data/confidential/in_situ_OSTIA_joined.csv
#
# Required columns:
#   site_id, date, site_temperature, OSTIA_temperature
#
# Optional manual exclusion file (stable site/date keys):
#   data/confidential/in_situ_manual_exclusions.csv
#
# QC:
#   - Primary isolated-spike threshold: 4.5 degrees C
#   - Sensitivity thresholds: 3.5, 4.0, 4.5, 5.0, 5.5 degrees C
#   - At least four valid neighbouring observations
#   - At least one valid observation before and after
#   - Median before/after difference <= 2 degrees C
#
# Validation metrics:
#   Pearson r, R2, mean difference (OSTIA - in situ), RMSE and cRMSD.
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------
# File paths and mode
# ------------------------------------------------------------------------------

in_situ_file <- file.path(
  "data", "confidential", "in_situ_OSTIA_joined.csv"
)

manual_exclusions_file <- file.path(
  "data", "confidential", "in_situ_manual_exclusions.csv"
)

real_data_available <- file.exists(in_situ_file)

if (real_data_available) {
  output_dir <- file.path("data", "derived", "in_situ_validation")
  figure_dir <- file.path("outputs", "figures")
  analysis_mode <- "confidential data"
} else {
  output_dir <- file.path("outputs", "in_situ_demo")
  figure_dir <- output_dir
  analysis_mode <- "synthetic demonstration"
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Create synthetic demonstration data
# ------------------------------------------------------------------------------

if (real_data_available) {

  in_situ_ostia <- read.csv(
    in_situ_file,
    stringsAsFactors = FALSE
  ) %>%
    mutate(date = as.Date(date))

} else {

  set.seed(812)

  mock_grid <- expand_grid(
    site_id = paste("Site", 1:8),
    date = seq.Date(
      as.Date("2024-01-01"),
      as.Date("2025-12-31"),
      by = "day"
    )
  ) %>%
    mutate(
      day = as.integer(format(date, "%j")),
      site_offset = rep(seq(-0.7, 0.7, length.out = 8), each = 731),
      site_temperature =
        10 +
        4.5 * sin(2 * pi * (day - 170) / 365.25) +
        site_offset +
        rnorm(n(), 0, 0.35),
      OSTIA_temperature =
        site_temperature +
        0.03 +
        rnorm(n(), 0, 0.55)
    ) %>%
    select(
      site_id,
      date,
      site_temperature,
      OSTIA_temperature
    )

  # Add missing periods and a small number of obvious spikes.
  mock_grid$site_temperature[
    runif(nrow(mock_grid)) < 0.28
  ] <- NA_real_

  spike_rows <- c(350, 1760, 4210)
  mock_grid$site_temperature[spike_rows] <- c(72, -12, 26)

  in_situ_ostia <- mock_grid

  message(
    "Confidential input not found. Running synthetic demonstration only; ",
    "outputs do not reproduce manuscript validation results."
  )
}

required_columns <- c(
  "site_id", "date", "site_temperature", "OSTIA_temperature"
)

stopifnot(
  all(required_columns %in% names(in_situ_ostia)),
  !anyDuplicated(in_situ_ostia[c("site_id", "date")])
)

# Optional stable site/date exclusions for known impossible source values.
if (file.exists(manual_exclusions_file)) {

  manual_exclusions <- read.csv(
    manual_exclusions_file,
    stringsAsFactors = FALSE
  ) %>%
    transmute(
      site_id = as.character(site_id),
      date = as.Date(date)
    )

} else {

  manual_exclusions <- data.frame(
    site_id = character(),
    date = as.Date(character())
  )
}

# ------------------------------------------------------------------------------
# QC helper functions
# ------------------------------------------------------------------------------

safe_median <- function(...) {
  x <- c(...)
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else median(x)
}

n_nonmissing <- function(...) {
  sum(!is.na(c(...)))
}

qc_in_situ_temperature <- function(
    data,
    spike_threshold = 4.5,
    temp_min = -2,
    temp_max = 30,
    context_threshold = 2,
    min_neighbours = 4,
    manual_exclusions = data.frame(
      site_id = character(),
      date = as.Date(character())
    )
) {

  exclusion_keys <- paste(
    manual_exclusions$site_id,
    manual_exclusions$date
  )

  data %>%
    mutate(
      site_id = as.character(site_id),
      date = as.Date(date),
      temp_original = site_temperature,
      observation_key = paste(site_id, date),
      flag_physical =
        observation_key %in% exclusion_keys |
        (
          !is.na(temp_original) &
          (temp_original < temp_min | temp_original > temp_max)
        ),
      temp_after_physical = if_else(
        flag_physical,
        NA_real_,
        temp_original
      )
    ) %>%
    arrange(site_id, date) %>%
    group_by(site_id) %>%
    mutate(
      lag1 = lag(temp_after_physical, 1),
      lag2 = lag(temp_after_physical, 2),
      lag3 = lag(temp_after_physical, 3),
      lead1 = lead(temp_after_physical, 1),
      lead2 = lead(temp_after_physical, 2),
      lead3 = lead(temp_after_physical, 3),

      before_median = pmap_dbl(
        list(lag1, lag2, lag3),
        safe_median
      ),

      after_median = pmap_dbl(
        list(lead1, lead2, lead3),
        safe_median
      ),

      local_median = pmap_dbl(
        list(lag1, lag2, lag3, lead1, lead2, lead3),
        safe_median
      ),

      n_before = pmap_int(
        list(lag1, lag2, lag3),
        n_nonmissing
      ),

      n_after = pmap_int(
        list(lead1, lead2, lead3),
        n_nonmissing
      ),

      n_neighbours = n_before + n_after,

      difference_from_context =
        temp_after_physical - local_median,

      flag_local_spike =
        !flag_physical &
        !is.na(temp_after_physical) &
        !is.na(local_median) &
        n_before >= 1 &
        n_after >= 1 &
        n_neighbours >= min_neighbours &
        abs(difference_from_context) > spike_threshold &
        abs(before_median - after_median) <= context_threshold,

      qc_flag = case_when(
        is.na(temp_original) ~ "missing_original",
        flag_physical ~ "physical_impossible",
        flag_local_spike ~ "isolated_spike",
        TRUE ~ "retained"
      ),

      site_temperature_qc = if_else(
        qc_flag == "retained",
        temp_after_physical,
        NA_real_
      )
    ) %>%
    ungroup() %>%
    select(
      -observation_key,
      -lag1, -lag2, -lag3,
      -lead1, -lead2, -lead3
    )
}

# ------------------------------------------------------------------------------
# Validation metrics
# ------------------------------------------------------------------------------

calculate_validation_metrics <- function(data, threshold) {

  paired <- data %>%
    filter(
      !is.na(site_temperature_qc),
      !is.na(OSTIA_temperature)
    )

  metric_summary <- function(df) {

    r <- cor(
      df$site_temperature_qc,
      df$OSTIA_temperature,
      use = "complete.obs"
    )

    bias <- mean(
      df$OSTIA_temperature - df$site_temperature_qc
    )

    rmse <- sqrt(mean(
      (df$OSTIA_temperature - df$site_temperature_qc)^2
    ))

    crmsd <- sqrt(mean(
      (
        (df$site_temperature_qc - mean(df$site_temperature_qc)) -
        (df$OSTIA_temperature - mean(df$OSTIA_temperature))
      )^2
    ))

    data.frame(
      n = nrow(df),
      pearson_r = r,
      r_squared = r^2,
      bias_degC = bias,
      rmse_degC = rmse,
      crmsd_degC = crmsd
    )
  }

  site_metrics <- paired %>%
    group_by(site_id) %>%
    group_modify(~ metric_summary(.x)) %>%
    ungroup() %>%
    mutate(spike_threshold = threshold, .before = 1)

  overall_metrics <- metric_summary(paired) %>%
    mutate(spike_threshold = threshold, .before = 1)

  list(
    site = site_metrics,
    overall = overall_metrics
  )
}

# ------------------------------------------------------------------------------
# Threshold sensitivity
# ------------------------------------------------------------------------------

spike_thresholds <- c(3.5, 4.0, 4.5, 5.0, 5.5)

qc_results <- map(
  spike_thresholds,
  ~ qc_in_situ_temperature(
    in_situ_ostia,
    spike_threshold = .x,
    manual_exclusions = manual_exclusions
  )
)

qc_summary <- map2_dfr(
  qc_results,
  spike_thresholds,
  ~ .x %>%
    summarise(
      spike_threshold = .y,
      n_nonmissing_original = sum(!is.na(temp_original)),
      n_removed_physical = sum(flag_physical, na.rm = TRUE),
      n_removed_local_spikes = sum(flag_local_spike, na.rm = TRUE),
      n_retained = sum(qc_flag == "retained", na.rm = TRUE)
    )
)

qc_summary_by_site <- map2_dfr(
  qc_results,
  spike_thresholds,
  ~ .x %>%
    group_by(site_id) %>%
    summarise(
      spike_threshold = .y,
      n_nonmissing_original = sum(!is.na(temp_original)),
      n_removed_physical = sum(flag_physical, na.rm = TRUE),
      n_removed_local_spikes = sum(flag_local_spike, na.rm = TRUE),
      n_retained = sum(qc_flag == "retained", na.rm = TRUE),
      .groups = "drop"
    )
)

validation_results <- map2(
  qc_results,
  spike_thresholds,
  calculate_validation_metrics
)

validation_site <- map_dfr(validation_results, "site")
validation_overall <- map_dfr(validation_results, "overall")

primary_threshold <- 4.5
primary_index <- match(primary_threshold, spike_thresholds)
primary_qc <- qc_results[[primary_index]]

primary_validation <- validation_site %>%
  filter(spike_threshold == primary_threshold) %>%
  select(-spike_threshold)

# ------------------------------------------------------------------------------
# Save tabular outputs
# ------------------------------------------------------------------------------

write.csv(
  qc_summary,
  file.path(output_dir, "QC_threshold_summary.csv"),
  row.names = FALSE
)

write.csv(
  qc_summary_by_site,
  file.path(output_dir, "QC_threshold_summary_by_site.csv"),
  row.names = FALSE
)

write.csv(
  validation_site,
  file.path(
    output_dir,
    "OSTIA_in_situ_validation_metrics_by_threshold.csv"
  ),
  row.names = FALSE
)

write.csv(
  validation_overall,
  file.path(
    output_dir,
    "OSTIA_in_situ_validation_overall_by_threshold.csv"
  ),
  row.names = FALSE
)

write.csv(
  primary_validation,
  file.path(
    output_dir,
    "OSTIA_in_situ_validation_metrics_4.5C.csv"
  ),
  row.names = FALSE
)

saveRDS(
  primary_qc,
  file.path(output_dir, "OSTIA_in_situ_joined_primary_QC.rds")
)

# ------------------------------------------------------------------------------
# OSTIA-versus-in-situ comparison
# ------------------------------------------------------------------------------

paired_plot_data <- primary_qc %>%
  select(site_id, date, site_temperature_qc, OSTIA_temperature) %>%
  pivot_longer(
    cols = c(site_temperature_qc, OSTIA_temperature),
    names_to = "source",
    values_to = "temperature"
  ) %>%
  mutate(
    source = recode(
      source,
      site_temperature_qc = "In situ",
      OSTIA_temperature = "OSTIA"
    )
  )

paired_plot <- ggplot(
  paired_plot_data,
  aes(date, temperature, colour = source)
) +
  geom_line(linewidth = 0.35, alpha = 0.75, na.rm = TRUE) +
  facet_wrap(~ site_id, ncol = 2, scales = "free_x") +
  labs(
    x = "Date",
    y = "Temperature (degrees C)",
    colour = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top")

ggsave(
  file.path(figure_dir, "in_situ_OSTIA_temperature_comparison.png"),
  paired_plot,
  width = 9,
  height = 8,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------------------------
# Sensitivity figure
# ------------------------------------------------------------------------------

p_removed <- ggplot(
  qc_summary,
  aes(spike_threshold, n_removed_local_spikes)
) +
  geom_vline(xintercept = primary_threshold, linetype = "dashed") +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = spike_thresholds) +
  labs(
    x = "Isolated-spike threshold (degrees C)",
    y = "Flagged isolated spikes (n)"
  ) +
  theme_classic(base_size = 11)

validation_long <- validation_overall %>%
  pivot_longer(
    cols = c(r_squared, bias_degC, rmse_degC, crmsd_degC),
    names_to = "metric",
    values_to = "value"
  )

p_metrics <- ggplot(
  validation_long,
  aes(spike_threshold, value)
) +
  geom_vline(xintercept = primary_threshold, linetype = "dashed") +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = spike_thresholds) +
  labs(
    x = "Isolated-spike threshold (degrees C)",
    y = NULL
  ) +
  theme_classic(base_size = 11)

sensitivity_figure <- p_removed + p_metrics +
  patchwork::plot_layout(widths = c(1, 1.8)) +
  patchwork::plot_annotation(tag_levels = "A")

ggsave(
  file.path(figure_dir, "in_situ_QC_sensitivity.png"),
  sensitivity_figure,
  width = 10,
  height = 5,
  dpi = 300,
  bg = "white"
)

cat("Mode:", analysis_mode, "\n")
print(qc_summary)
print(validation_overall)

if (!real_data_available) {
  message(
    "Synthetic demo outputs are illustrative only and must not be used as ",
    "manuscript results."
  )
}
