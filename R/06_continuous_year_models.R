# ==============================================================================
# 06 - Continuous-year mixed models and temporal trends
# ==============================================================================
# Purpose:
#   Test temporal change in MHW frequency, duration, mean intensity, cumulative
#   intensity and severity across the two clusters.
#
# Analysis unit:
#   - Frequency: OSTIA cell-year, including zero-event years.
#   - Event characteristics/severity: individual MHW events.
#
# Model structure:
#   response ~ cluster * continuous year + (1 | OSTIA cell) + (1 | year)
#
# Inputs:
#   data/derived/cell_exposure/
#     scotland_OSTIA_cell_mhw_events_clustered.rds
#
# Main outputs:
#   data/derived/cell_models/
#     continuous_year_trends.csv
#     continuous_year_slope_contrasts.csv
#     frequency_poisson_negative_binomial_sensitivity.csv
#     MHW_continuous_year_models.rds
#   outputs/figures/
#     figure_continuous_year_trends.png
#     continuous_year_LMM_diagnostics.png
#     continuous_year_GLMM_diagnostics.png
# ==============================================================================

library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
library(DHARMa)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------

events_file <- file.path(
  "data", "derived", "cell_exposure",
  "scotland_OSTIA_cell_mhw_events_clustered.rds"
)

output_dir <- file.path("data", "derived", "cell_models")
figure_dir <- file.path("outputs", "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

models_file <- file.path(output_dir, "MHW_continuous_year_models.rds")
trends_file <- file.path(output_dir, "continuous_year_trends.csv")
contrasts_file <- file.path(output_dir, "continuous_year_slope_contrasts.csv")
interaction_tests_file <- file.path(output_dir, "continuous_year_interaction_tests.csv")
frequency_sensitivity_file <- file.path(
  output_dir, "frequency_poisson_negative_binomial_sensitivity.csv"
)
annual_summary_file <- file.path(
  output_dir, "cluster_annual_MHW_summary.csv"
)
prediction_file <- file.path(
  output_dir, "continuous_year_model_predictions.csv"
)

trend_figure_file <- file.path(
  figure_dir, "figure_continuous_year_trends.png"
)
lmm_diag_file <- file.path(
  figure_dir, "continuous_year_LMM_diagnostics.png"
)
glmm_diag_file <- file.path(
  figure_dir, "continuous_year_GLMM_diagnostics.png"
)

# ------------------------------------------------------------------------------
# Prepare modelling datasets
# ------------------------------------------------------------------------------

events <- readRDS(events_file) %>%
  ungroup() %>%
  mutate(
    date_peak = as.Date(date_peak),
    year = as.integer(format(date_peak, "%Y")),
    year_c = year - 2000,
    year_factor = factor(year),
    cluster = factor(cluster, levels = c("Cluster 1", "Cluster 2")),
    ostia_cell = factor(ostia_cell),
    log_duration = log1p(duration),
    log_cumulative_intensity = log1p(intensity_cumulative),
    severe_extreme = as.integer(
      category %in% c("III Severe", "IV Extreme")
    )
  ) %>%
  filter(
    !is.na(year),
    !is.na(cluster),
    !is.na(ostia_cell)
  )

stopifnot(
  n_distinct(events$ostia_cell) == 218,
  nrow(events) == 21099,
  range(events$year) == c(1990, 2025)
)

cell_cluster_lookup <- events %>%
  distinct(ostia_cell, cluster)

frequency_data <- expand_grid(
  cell_cluster_lookup,
  year = 1990:2025
) %>%
  mutate(
    year_dec = (year - 2000) / 10,
    year_factor = factor(year),
    ostia_cell = factor(ostia_cell),
    cluster = factor(cluster, levels = c("Cluster 1", "Cluster 2"))
  ) %>%
  left_join(
    events %>% count(ostia_cell, year, name = "n_events"),
    by = c("ostia_cell", "year")
  ) %>%
  mutate(n_events = replace_na(n_events, 0L))

stopifnot(nrow(frequency_data) == 218 * 36)

# ------------------------------------------------------------------------------
# Fit continuous-year models
# ------------------------------------------------------------------------------

model_frequency <- glmer(
  n_events ~ cluster * year_dec +
    (1 | ostia_cell) +
    (1 | year_factor),
  data = frequency_data,
  family = poisson(link = "log"),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

# Negative-binomial (sensitivity model check)
model_frequency_nb <- glmer.nb(
  n_events ~ cluster * year_dec +
    (1 | ostia_cell) +
    (1 | year_factor),
  data = frequency_data,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

model_duration <- lmer(
  log_duration ~ cluster * year_c +
    (1 | ostia_cell) +
    (1 | year_factor),
  data = events,
  REML = TRUE
)

model_intensity <- lmer(
  intensity_mean ~ cluster * year_c +
    (1 | ostia_cell) +
    (1 | year_factor),
  data = events,
  REML = TRUE
)

model_cumulative_intensity <- lmer(
  log_cumulative_intensity ~ cluster * year_c +
    (1 | ostia_cell) +
    (1 | year_factor),
  data = events,
  REML = TRUE
)

model_severity <- glmer(
  severe_extreme ~ cluster * year_c +
    (1 | ostia_cell) +
    (1 | year_factor),
  data = events,
  family = binomial(link = "logit"),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

# ------------------------------------------------------------------------------
# Fixed-effect tests
# ------------------------------------------------------------------------------

anova_to_df <- function(x, response) {
  out <- as.data.frame(x)
  data.frame(
    response = response,
    term = rownames(out),
    out,
    row.names = NULL,
    check.names = FALSE
  )
}

frequency_drop <- drop1(model_frequency, test = "Chisq")
frequency_nb_drop <- drop1(model_frequency_nb, test = "Chisq")

anova_duration <- anova(
  model_duration,
  type = 3,
  ddf = "Satterthwaite"
)

anova_intensity <- anova(
  model_intensity,
  type = 3,
  ddf = "Satterthwaite"
)

anova_cumulative <- anova(
  model_cumulative_intensity,
  type = 3,
  ddf = "Satterthwaite"
)

frequency_p_col <- grep("^Pr", names(frequency_drop), value = TRUE)[1]
severity_coef <- coef(summary(model_severity))
severity_row <- "clusterCluster 2:year_c"

interaction_tests <- bind_rows(
  data.frame(
    response = "Frequency",
    statistic_type = "likelihood-ratio Chi-square",
    statistic = frequency_drop["cluster:year_dec", "LRT"],
    df1 = NA_real_,
    df2 = NA_real_,
    p_value = frequency_drop["cluster:year_dec", frequency_p_col]
  ),
  data.frame(
    response = "Duration",
    statistic_type = "F",
    statistic = anova_duration["cluster:year_c", "F value"],
    df1 = anova_duration["cluster:year_c", "NumDF"],
    df2 = anova_duration["cluster:year_c", "DenDF"],
    p_value = anova_duration["cluster:year_c", "Pr(>F)"]
  ),
  data.frame(
    response = "Mean intensity",
    statistic_type = "F",
    statistic = anova_intensity["cluster:year_c", "F value"],
    df1 = anova_intensity["cluster:year_c", "NumDF"],
    df2 = anova_intensity["cluster:year_c", "DenDF"],
    p_value = anova_intensity["cluster:year_c", "Pr(>F)"]
  ),
  data.frame(
    response = "Cumulative intensity",
    statistic_type = "F",
    statistic = anova_cumulative["cluster:year_c", "F value"],
    df1 = anova_cumulative["cluster:year_c", "NumDF"],
    df2 = anova_cumulative["cluster:year_c", "DenDF"],
    p_value = anova_cumulative["cluster:year_c", "Pr(>F)"]
  ),
  data.frame(
    response = "Severity",
    statistic_type = "Wald z",
    statistic = severity_coef[severity_row, "z value"],
    df1 = NA_real_,
    df2 = NA_real_,
    p_value = severity_coef[severity_row, "Pr(>|z|)"]
  )
)

write.csv(
  interaction_tests,
  interaction_tests_file,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Cluster-specific slopes
# ------------------------------------------------------------------------------

normalise_trend <- function(emm_object, metric) {

  x <- as.data.frame(summary(emm_object, infer = c(TRUE, TRUE)))

  trend_col <- grep("\\.trend$", names(x), value = TRUE)[1]
  lower_col <- intersect(c("lower.CL", "asymp.LCL"), names(x))[1]
  upper_col <- intersect(c("upper.CL", "asymp.UCL"), names(x))[1]
  stat_col <- intersect(c("t.ratio", "z.ratio"), names(x))[1]

  data.frame(
    metric = metric,
    cluster = x$cluster,
    trend = x[[trend_col]],
    SE = x$SE,
    df = if ("df" %in% names(x)) x$df else NA_real_,
    lower = x[[lower_col]],
    upper = x[[upper_col]],
    statistic = x[[stat_col]],
    p_value = x$p.value
  )
}

normalise_contrast <- function(contrast_object, metric) {

  x <- as.data.frame(summary(contrast_object, infer = c(TRUE, TRUE)))

  lower_col <- intersect(c("lower.CL", "asymp.LCL"), names(x))[1]
  upper_col <- intersect(c("upper.CL", "asymp.UCL"), names(x))[1]
  stat_col <- intersect(c("t.ratio", "z.ratio"), names(x))[1]

  data.frame(
    metric = metric,
    contrast = x$contrast,
    estimate = x$estimate,
    SE = x$SE,
    df = if ("df" %in% names(x)) x$df else NA_real_,
    lower = x[[lower_col]],
    upper = x[[upper_col]],
    statistic = x[[stat_col]],
    p_value = x$p.value
  )
}

frequency_trends <- emtrends(model_frequency, ~ cluster, var = "year_dec")
duration_trends <- emtrends(
  model_duration, ~ cluster, var = "year_c", lmer.df = "asymptotic"
)
intensity_trends <- emtrends(
  model_intensity, ~ cluster, var = "year_c", lmer.df = "asymptotic"
)
cumulative_trends <- emtrends(
  model_cumulative_intensity, ~ cluster, var = "year_c",
  lmer.df = "asymptotic"
)
severity_trends <- emtrends(model_severity, ~ cluster, var = "year_c")

trend_table <- bind_rows(
  normalise_trend(frequency_trends, "Frequency"),
  normalise_trend(duration_trends, "Duration"),
  normalise_trend(intensity_trends, "Mean intensity"),
  normalise_trend(cumulative_trends, "Cumulative intensity"),
  normalise_trend(severity_trends, "Severity")
) %>%
  mutate(
    annual_multiplier = case_when(
      metric == "Frequency" ~ exp(trend / 10),
      metric == "Duration" ~ exp(trend),
      metric == "Cumulative intensity" ~ exp(trend),
      metric == "Severity" ~ exp(trend),
      TRUE ~ NA_real_
    ),
    annual_percent_change = case_when(
      metric %in% c("Frequency", "Duration", "Cumulative intensity") ~
        100 * (annual_multiplier - 1),
      TRUE ~ NA_real_
    ),
    interpretation = case_when(
      metric == "Frequency" ~ "annual rate ratio",
      metric == "Duration" ~ "multiplicative change in duration + 1",
      metric == "Mean intensity" ~ "degrees C per year",
      metric == "Cumulative intensity" ~
        "multiplicative change in cumulative intensity + 1",
      metric == "Severity" ~ "annual odds ratio"
    )
  )

slope_contrasts <- bind_rows(
  normalise_contrast(
    contrast(frequency_trends, method = "revpairwise"),
    "Frequency"
  ),
  normalise_contrast(
    contrast(duration_trends, method = "revpairwise"),
    "Duration"
  ),
  normalise_contrast(
    contrast(intensity_trends, method = "revpairwise"),
    "Mean intensity"
  ),
  normalise_contrast(
    contrast(cumulative_trends, method = "revpairwise"),
    "Cumulative intensity"
  ),
  normalise_contrast(
    contrast(severity_trends, method = "revpairwise"),
    "Severity"
  )
)

write.csv(trend_table, trends_file, row.names = FALSE)
write.csv(slope_contrasts, contrasts_file, row.names = FALSE)

# ------------------------------------------------------------------------------
# Poisson versus negative-binomial sensitivity for frequency
# ------------------------------------------------------------------------------

frequency_nb_trends <- emtrends(
  model_frequency_nb,
  ~ cluster,
  var = "year_dec"
)

extract_interaction_lrt <- function(drop_table) {
  p_col <- grep("^Pr", names(drop_table), value = TRUE)[1]
  data.frame(
    interaction_LRT = drop_table["cluster:year_dec", "LRT"],
    interaction_p = drop_table["cluster:year_dec", p_col]
  )
}

poisson_lrt <- extract_interaction_lrt(frequency_drop)
nb_lrt <- extract_interaction_lrt(frequency_nb_drop)

poisson_slopes <- normalise_trend(frequency_trends, "Frequency") %>%
  transmute(
    model = "Poisson",
    cluster,
    annual_rate_ratio = exp(trend / 10),
    annual_percent_change = 100 * (annual_rate_ratio - 1)
  )

nb_slopes <- normalise_trend(frequency_nb_trends, "Frequency") %>%
  transmute(
    model = "Negative binomial",
    cluster,
    annual_rate_ratio = exp(trend / 10),
    annual_percent_change = 100 * (annual_rate_ratio - 1)
  )

frequency_sensitivity <- bind_rows(poisson_slopes, nb_slopes) %>%
  left_join(
    bind_rows(
      data.frame(model = "Poisson", poisson_lrt),
      data.frame(model = "Negative binomial", nb_lrt)
    ),
    by = "model"
  ) %>%
  mutate(AIC = if_else(
    model == "Poisson",
    AIC(model_frequency),
    AIC(model_frequency_nb)
  ))

write.csv(
  frequency_sensitivity,
  frequency_sensitivity_file,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Observed annual summaries and fixed-effect model trends
# ------------------------------------------------------------------------------

cluster_sizes <- events %>%
  distinct(ostia_cell, cluster) %>%
  count(cluster, name = "n_total_cells")

year_cluster_grid <- expand_grid(
  year = 1990:2025,
  cluster = factor(
    c("Cluster 1", "Cluster 2"),
    levels = c("Cluster 1", "Cluster 2")
  )
)

annual_summary <- events %>%
  group_by(year, cluster) %>%
  summarise(
    n_events = n(),
    mean_duration = mean(duration),
    mean_intensity = mean(intensity_mean),
    mean_cumulative_intensity = mean(intensity_cumulative),
    severe_extreme_probability = mean(severe_extreme),
    .groups = "drop"
  ) %>%
  right_join(year_cluster_grid, by = c("year", "cluster")) %>%
  left_join(cluster_sizes, by = "cluster") %>%
  mutate(
    n_events = replace_na(n_events, 0L),
    events_per_cell = n_events / n_total_cells
  ) %>%
  arrange(cluster, year)

write.csv(annual_summary, annual_summary_file, row.names = FALSE)

find_column <- function(x, candidates) {
  out <- intersect(candidates, names(x))
  if (length(out) == 0) {
    stop("Expected column not found: ", paste(candidates, collapse = ", "))
  }
  out[1]
}

make_lmm_predictions <- function(model, back_transform = FALSE) {

  x <- as.data.frame(
    emmeans(
      model,
      ~ cluster | year_c,
      at = list(year_c = 1990:2025 - 2000),
      lmer.df = "asymptotic"
    )
  )

  lower_col <- find_column(x, c("lower.CL", "asymp.LCL"))
  upper_col <- find_column(x, c("upper.CL", "asymp.UCL"))

  out <- x %>%
    transmute(
      cluster,
      year = year_c + 2000,
      estimate = emmean,
      lower = .data[[lower_col]],
      upper = .data[[upper_col]]
    )

  if (back_transform) {
    out <- out %>%
      mutate(
        estimate = exp(estimate) - 1,
        lower = exp(lower) - 1,
        upper = exp(upper) - 1
      )
  }

  out
}

make_glmm_predictions <- function(
    model,
    time_variable,
    at_values
) {

  formula_text <- paste0("~ cluster | ", time_variable)

  x <- as.data.frame(
    emmeans(
      model,
      specs = as.formula(formula_text),
      at = setNames(list(at_values), time_variable),
      type = "response"
    )
  )

  value_col <- find_column(x, c("rate", "prob", "response", "emmean"))
  lower_col <- find_column(x, c("lower.CL", "asymp.LCL"))
  upper_col <- find_column(x, c("upper.CL", "asymp.UCL"))

  time_values <- x[[time_variable]]

  data.frame(
    cluster = x$cluster,
    year = if (time_variable == "year_dec") {
      time_values * 10 + 2000
    } else {
      time_values + 2000
    },
    estimate = x[[value_col]],
    lower = x[[lower_col]],
    upper = x[[upper_col]]
  )
}

trend_years <- 1990:2025

frequency_predictions <- make_glmm_predictions(
  model_frequency,
  "year_dec",
  (trend_years - 2000) / 10
) %>%
  mutate(metric = "Frequency")

duration_predictions <- make_lmm_predictions(
  model_duration,
  back_transform = TRUE
) %>%
  mutate(metric = "Duration")

intensity_predictions <- make_lmm_predictions(
  model_intensity,
  back_transform = FALSE
) %>%
  mutate(metric = "Mean intensity")

cumulative_predictions <- make_lmm_predictions(
  model_cumulative_intensity,
  back_transform = TRUE
) %>%
  mutate(metric = "Cumulative intensity")

severity_predictions <- make_glmm_predictions(
  model_severity,
  "year_c",
  trend_years - 2000
) %>%
  mutate(metric = "Severity")

model_predictions <- bind_rows(
  frequency_predictions,
  duration_predictions,
  intensity_predictions,
  cumulative_predictions,
  severity_predictions
)

write.csv(model_predictions, prediction_file, row.names = FALSE)

observed_long <- bind_rows(
  annual_summary %>%
    transmute(year, cluster, metric = "Frequency", value = events_per_cell),
  annual_summary %>%
    transmute(year, cluster, metric = "Duration", value = mean_duration),
  annual_summary %>%
    transmute(year, cluster, metric = "Mean intensity", value = mean_intensity),
  annual_summary %>%
    transmute(
      year, cluster, metric = "Cumulative intensity",
      value = mean_cumulative_intensity
    ),
  annual_summary %>%
    transmute(
      year, cluster, metric = "Severity",
      value = severe_extreme_probability
    )
)

metric_levels <- c(
  "Frequency",
  "Mean intensity",
  "Duration",
  "Cumulative intensity",
  "Severity"
)

observed_long$metric <- factor(observed_long$metric, levels = metric_levels)
model_predictions$metric <- factor(
  model_predictions$metric,
  levels = metric_levels
)

trend_plot <- ggplot(
  model_predictions,
  aes(year, estimate, colour = cluster, fill = cluster)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(
    data = observed_long,
    aes(year, value, colour = cluster),
    inherit.aes = FALSE,
    size = 0.8,
    alpha = 0.45
  ) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(
    x = "Year",
    y = NULL,
    colour = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_blank()
  )

ggsave(
  trend_figure_file,
  trend_plot,
  width = 9,
  height = 9,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------------------------
# Model diagnostics
# ------------------------------------------------------------------------------

png(lmm_diag_file, width = 1800, height = 2400, res = 250)
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))

for (obj in list(model_intensity, model_duration, model_cumulative_intensity)) {
  plot(
    fitted(obj),
    residuals(obj),
    xlab = "Fitted values",
    ylab = "Residuals",
    pch = 16,
    cex = 0.45
  )
  abline(h = 0, lty = 2)

  qqnorm(residuals(obj), pch = 16, cex = 0.45)
  qqline(residuals(obj))
}

dev.off()

set.seed(812)

sim_frequency <- DHARMa::simulateResiduals(model_frequency, n = 1000)
sim_severity <- DHARMa::simulateResiduals(model_severity, n = 1000)

png(glmm_diag_file, width = 1800, height = 1800, res = 250)
par(mfrow = c(2, 2))

DHARMa::plotQQunif(
  sim_frequency,
  testUniformity = FALSE,
  testOutliers = FALSE,
  testDispersion = FALSE
)
DHARMa::plotResiduals(sim_frequency)

DHARMa::plotQQunif(
  sim_severity,
  testUniformity = FALSE,
  testOutliers = FALSE,
  testDispersion = FALSE
)
DHARMa::plotResiduals(sim_severity)

dev.off()

print(performance::check_overdispersion(model_frequency))
print(performance::check_singularity(model_frequency))
print(performance::check_singularity(model_duration))
print(performance::check_singularity(model_intensity))
print(performance::check_singularity(model_cumulative_intensity))
print(performance::check_singularity(model_severity))

print(DHARMa::testUniformity(sim_frequency))
print(DHARMa::testDispersion(sim_frequency))
print(DHARMa::testUniformity(sim_severity))
print(DHARMa::testDispersion(sim_severity))

# ------------------------------------------------------------------------------
# Save fitted models
# ------------------------------------------------------------------------------

saveRDS(
  list(
    frequency = model_frequency,
    frequency_negative_binomial = model_frequency_nb,
    duration = model_duration,
    intensity = model_intensity,
    cumulative_intensity = model_cumulative_intensity,
    severity = model_severity
  ),
  models_file
)

print(trend_table)
print(slope_contrasts)
print(frequency_sensitivity)
