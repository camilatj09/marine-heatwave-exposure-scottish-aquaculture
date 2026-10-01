# ==============================================================================
# 07 - Pre-/post-2020 models and alternative-breakpoint sensitivity
# ==============================================================================
# Purpose:
#   Quantify recent changes in event characteristics using a 2020 boundary and
#   test whether the cluster contrasts are robust to alternative boundaries in
#   2017, 2018 and 2019.
#
# Input:
#   data/derived/cell_exposure/
#     scotland_OSTIA_cell_mhw_events_clustered.rds
#
# Main outputs:
#   data/derived/cell_models/
#     period_2020_estimates.csv
#     period_2020_cluster_changes.csv
#     period_2020_interaction_tests.csv
#     breakpoint_sensitivity_interactions_2017_2020.csv
#     breakpoint_sensitivity_cluster_changes_2017_2020.csv
#     MHW_period_2020_models.rds
#   outputs/figures/
#     figure_period_2020_estimates.png
#     supplementary_breakpoint_sensitivity.png
#     period_2020_LMM_diagnostics.png
#     period_2020_GLMM_diagnostics.png
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

models_file <- file.path(output_dir, "MHW_period_2020_models.rds")
estimates_file <- file.path(output_dir, "period_2020_estimates.csv")
changes_file <- file.path(output_dir, "period_2020_cluster_changes.csv")
interactions_file <- file.path(output_dir, "period_2020_interaction_tests.csv")

breakpoint_interactions_file <- file.path(
  output_dir, "breakpoint_sensitivity_interactions_2017_2020.csv"
)
breakpoint_changes_file <- file.path(
  output_dir, "breakpoint_sensitivity_cluster_changes_2017_2020.csv"
)

period_figure_file <- file.path(
  figure_dir, "figure_period_2020_estimates.png"
)
breakpoint_figure_file <- file.path(
  figure_dir, "supplementary_breakpoint_sensitivity.png"
)
lmm_diag_file <- file.path(
  figure_dir, "period_2020_LMM_diagnostics.png"
)
glmm_diag_file <- file.path(
  figure_dir, "period_2020_GLMM_diagnostics.png"
)

# ------------------------------------------------------------------------------
# Prepare event-level dataset
# ------------------------------------------------------------------------------

events <- readRDS(events_file) %>%
  ungroup() %>%
  mutate(
    date_peak = as.Date(date_peak),
    year = as.integer(format(date_peak, "%Y")),
    year_factor = factor(year),
    cluster = factor(cluster, levels = c("Cluster 1", "Cluster 2")),
    ostia_cell = factor(ostia_cell),
    log_duration = log1p(duration),
    log_cumulative_intensity = log1p(intensity_cumulative),
    severe_extreme = as.integer(
      category %in% c("III Severe", "IV Extreme")
    )
  )

stopifnot(
  n_distinct(events$ostia_cell) == 218,
  nrow(events) == 21099
)

# ------------------------------------------------------------------------------
# Model and summary helpers
# ------------------------------------------------------------------------------

fit_period_models <- function(data, breakpoint) {

  dat <- data %>%
    mutate(
      period = if_else(year < breakpoint, "pre", "post"),
      period = factor(period, levels = c("pre", "post"))
    )

  list(
    data = dat,

    intensity = lmer(
      intensity_mean ~ cluster * period +
        (1 | ostia_cell) +
        (1 | year_factor),
      data = dat,
      REML = TRUE
    ),

    duration = lmer(
      log_duration ~ cluster * period +
        (1 | ostia_cell) +
        (1 | year_factor),
      data = dat,
      REML = TRUE
    ),

    cumulative_intensity = lmer(
      log_cumulative_intensity ~ cluster * period +
        (1 | ostia_cell) +
        (1 | year_factor),
      data = dat,
      REML = TRUE
    ),

    severity = glmer(
      severe_extreme ~ cluster * period +
        (1 | ostia_cell) +
        (1 | year_factor),
      data = dat,
      family = binomial(link = "logit"),
      control = glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    )
  )
}

find_column <- function(x, candidates) {
  out <- intersect(candidates, names(x))
  if (length(out) == 0) {
    stop("Expected column not found: ", paste(candidates, collapse = ", "))
  }
  out[1]
}

summarise_lmm_period <- function(model, response, back_transform = FALSE) {

  emm <- emmeans(
    model,
    ~ cluster * period,
    lmer.df = "asymptotic"
  )

  estimates <- as.data.frame(emm)

  lower_col <- find_column(estimates, c("lower.CL", "asymp.LCL"))
  upper_col <- find_column(estimates, c("upper.CL", "asymp.UCL"))

  estimates <- estimates %>%
    transmute(
      response = response,
      cluster,
      period,
      estimate = emmean,
      lower = .data[[lower_col]],
      upper = .data[[upper_col]]
    )

  if (back_transform) {
    estimates <- estimates %>%
      mutate(
        estimate = exp(estimate) - 1,
        lower = exp(lower) - 1,
        upper = exp(upper) - 1
      )
  }

  contrast_raw <- as.data.frame(
    contrast(
      emm,
      method = "revpairwise",
      by = "cluster"
    )
  )

  statistic_col <- find_column(
    contrast_raw,
    c("t.ratio", "z.ratio")
  )

  contrast_df <- data.frame(
    response = response,
    cluster = contrast_raw$cluster,
    contrast = contrast_raw$contrast,
    estimate_model_scale = contrast_raw$estimate,
    SE = contrast_raw$SE,
    df = if ("df" %in% names(contrast_raw)) {
      contrast_raw$df
    } else {
      NA_real_
    },
    statistic = contrast_raw[[statistic_col]],
    p_value = contrast_raw$p.value
  )

  change_df <- estimates %>%
    select(response, cluster, period, estimate) %>%
    pivot_wider(names_from = period, values_from = estimate) %>%
    mutate(
      percent_change = 100 * (post - pre) / pre
    ) %>%
    left_join(
      contrast_df %>% select(response, cluster, p_value),
      by = c("response", "cluster")
    )

  list(
    estimates = estimates,
    contrasts = contrast_df,
    changes = change_df
  )
}

summarise_severity_period <- function(model) {

  emm_response <- emmeans(
    model,
    ~ cluster * period,
    type = "response"
  )

  estimates_raw <- as.data.frame(emm_response)

  value_col <- find_column(estimates_raw, c("prob", "response"))
  lower_col <- find_column(estimates_raw, c("lower.CL", "asymp.LCL"))
  upper_col <- find_column(estimates_raw, c("upper.CL", "asymp.UCL"))

  estimates <- estimates_raw %>%
    transmute(
      response = "Severity",
      cluster,
      period,
      estimate = .data[[value_col]],
      lower = .data[[lower_col]],
      upper = .data[[upper_col]]
    )

  contrasts <- as.data.frame(
    contrast(
      emmeans(model, ~ period | cluster),
      method = "revpairwise",
      type = "response"
    )
  ) %>%
    transmute(
      response = "Severity",
      cluster,
      contrast,
      odds_ratio = odds.ratio,
      SE,
      df = NA_real_,
      statistic = z.ratio,
      p_value = p.value
    )

  changes <- estimates %>%
    select(response, cluster, period, estimate) %>%
    pivot_wider(names_from = period, values_from = estimate) %>%
    left_join(
      contrasts %>%
        select(response, cluster, odds_ratio, p_value),
      by = c("response", "cluster")
    ) %>%
    mutate(percent_change = NA_real_)

  list(
    estimates = estimates,
    contrasts = contrasts,
    changes = changes
  )
}

extract_interaction_test <- function(model, response, model_type) {

  if (model_type == "lmm") {

    x <- anova(model, type = 3, ddf = "Satterthwaite")

    data.frame(
      response = response,
      estimate = fixef(model)[["clusterCluster 2:periodpost"]],
      statistic = x["cluster:period", "F value"],
      df1 = x["cluster:period", "NumDF"],
      df2 = x["cluster:period", "DenDF"],
      p_value = x["cluster:period", "Pr(>F)"]
    )

  } else {

    x <- coef(summary(model))
    row_name <- "clusterCluster 2:periodpost"

    data.frame(
      response = response,
      estimate = x[row_name, "Estimate"],
      statistic = x[row_name, "z value"],
      df1 = NA_real_,
      df2 = NA_real_,
      p_value = x[row_name, "Pr(>|z|)"]
    )
  }
}

# ------------------------------------------------------------------------------
# Primary pre-/post-2020 analysis
# ------------------------------------------------------------------------------

models_2020 <- fit_period_models(events, 2020)

intensity_2020 <- summarise_lmm_period(
  models_2020$intensity,
  "Mean intensity",
  back_transform = FALSE
)

duration_2020 <- summarise_lmm_period(
  models_2020$duration,
  "Duration",
  back_transform = TRUE
)

cumulative_2020 <- summarise_lmm_period(
  models_2020$cumulative_intensity,
  "Cumulative intensity",
  back_transform = TRUE
)

severity_2020 <- summarise_severity_period(models_2020$severity)

period_estimates <- bind_rows(
  intensity_2020$estimates,
  duration_2020$estimates,
  cumulative_2020$estimates,
  severity_2020$estimates
)

period_changes <- bind_rows(
  intensity_2020$changes %>% mutate(odds_ratio = NA_real_),
  duration_2020$changes %>% mutate(odds_ratio = NA_real_),
  cumulative_2020$changes %>% mutate(odds_ratio = NA_real_),
  severity_2020$changes
) %>%
  select(
    response, cluster, pre, post,
    percent_change, odds_ratio, p_value
  )

period_interactions <- bind_rows(
  extract_interaction_test(
    models_2020$intensity,
    "Mean intensity",
    "lmm"
  ),
  extract_interaction_test(
    models_2020$duration,
    "Duration",
    "lmm"
  ),
  extract_interaction_test(
    models_2020$cumulative_intensity,
    "Cumulative intensity",
    "lmm"
  ),
  extract_interaction_test(
    models_2020$severity,
    "Severity",
    "glmm"
  )
)

write.csv(period_estimates, estimates_file, row.names = FALSE)
write.csv(period_changes, changes_file, row.names = FALSE)
write.csv(period_interactions, interactions_file, row.names = FALSE)

saveRDS(
  list(
    intensity = models_2020$intensity,
    duration = models_2020$duration,
    cumulative_intensity = models_2020$cumulative_intensity,
    severity = models_2020$severity
  ),
  models_file
)

# ------------------------------------------------------------------------------
#  Model estimates figure
# ------------------------------------------------------------------------------

response_levels <- c(
  "Mean intensity",
  "Duration",
  "Cumulative intensity",
  "Severity"
)

period_estimates$response <- factor(
  period_estimates$response,
  levels = response_levels
)

period_plot <- ggplot(
  period_estimates,
  aes(period, estimate, colour = cluster, group = cluster)
) +
  geom_line(linewidth = 0.6, position = position_dodge(width = 0.15)) +
  geom_point(size = 2.3, position = position_dodge(width = 0.15)) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.08,
    position = position_dodge(width = 0.15)
  ) +
  facet_wrap(~ response, scales = "free_y", ncol = 2) +
  scale_x_discrete(labels = c(pre = "1990-2019", post = "2020-2025")) +
  labs(x = NULL, y = NULL, colour = NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_blank()
  )

ggsave(
  period_figure_file,
  period_plot,
  width = 8.5,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------------------------
# Alternative-breakpoint sensitivity: 2017-2020
# ------------------------------------------------------------------------------

candidate_breakpoints <- 2017:2020

breakpoint_interactions <- list()
breakpoint_changes <- list()

for (bp in candidate_breakpoints) {

  fitted <- fit_period_models(events, bp)

  intensity <- summarise_lmm_period(
    fitted$intensity,
    "Mean intensity",
    back_transform = FALSE
  )

  duration <- summarise_lmm_period(
    fitted$duration,
    "Duration",
    back_transform = TRUE
  )

  cumulative <- summarise_lmm_period(
    fitted$cumulative_intensity,
    "Cumulative intensity",
    back_transform = TRUE
  )

  severity <- summarise_severity_period(fitted$severity)

  breakpoint_interactions[[as.character(bp)]] <- bind_rows(
    extract_interaction_test(fitted$intensity, "Mean intensity", "lmm"),
    extract_interaction_test(fitted$duration, "Duration", "lmm"),
    extract_interaction_test(
      fitted$cumulative_intensity,
      "Cumulative intensity",
      "lmm"
    ),
    extract_interaction_test(fitted$severity, "Severity", "glmm")
  ) %>%
    mutate(
      breakpoint = bp,
      post_years = 2025 - bp + 1,
      .before = response
    )

  breakpoint_changes[[as.character(bp)]] <- bind_rows(
    intensity$changes %>% mutate(odds_ratio = NA_real_),
    duration$changes %>% mutate(odds_ratio = NA_real_),
    cumulative$changes %>% mutate(odds_ratio = NA_real_),
    severity$changes
  ) %>%
    mutate(
      breakpoint = bp,
      post_years = 2025 - bp + 1,
      .before = response
    ) %>%
    select(
      breakpoint, post_years, response, cluster,
      pre, post, percent_change, odds_ratio, p_value
    )
}

breakpoint_interactions_df <- bind_rows(breakpoint_interactions) %>%
  arrange(response, breakpoint)

breakpoint_changes_df <- bind_rows(breakpoint_changes) %>%
  arrange(response, cluster, breakpoint)

write.csv(
  breakpoint_interactions_df,
  breakpoint_interactions_file,
  row.names = FALSE
)

write.csv(
  breakpoint_changes_df,
  breakpoint_changes_file,
  row.names = FALSE
)

# Simple visual summary of robustness across the four boundaries.
breakpoint_plot_data <- breakpoint_changes_df %>%
  mutate(
    change = if_else(
      response == "Severity",
      odds_ratio,
      percent_change
    ),
    response_label = if_else(
      response == "Severity",
      "Severity (pre/post odds ratio)",
      paste0(response, " (% change)")
    )
  )

breakpoint_plot <- ggplot(
  breakpoint_plot_data,
  aes(breakpoint, change, colour = cluster)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  facet_wrap(~ response_label, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = candidate_breakpoints) +
  labs(
    x = "First year of post period",
    y = NULL,
    colour = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    strip.background = element_blank()
  )

ggsave(
  breakpoint_figure_file,
  breakpoint_plot,
  width = 8.5,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------------------------
# Diagnostics for the primary 2020 models
# ------------------------------------------------------------------------------

png(lmm_diag_file, width = 1800, height = 2400, res = 250)
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))

for (obj in list(
  models_2020$intensity,
  models_2020$duration,
  models_2020$cumulative_intensity
)) {
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
sim_severity_2020 <- DHARMa::simulateResiduals(
  models_2020$severity,
  n = 1000
)

png(glmm_diag_file, width = 1800, height = 900, res = 250)
par(mfrow = c(1, 2))

DHARMa::plotQQunif(
  sim_severity_2020,
  testUniformity = FALSE,
  testOutliers = FALSE,
  testDispersion = FALSE
)
DHARMa::plotResiduals(sim_severity_2020)

dev.off()

print(DHARMa::testUniformity(sim_severity_2020))
print(DHARMa::testDispersion(sim_severity_2020))
print(performance::check_overdispersion(models_2020$severity))

for (obj in list(
  models_2020$intensity,
  models_2020$duration,
  models_2020$cumulative_intensity,
  models_2020$severity
)) {
  print(performance::check_singularity(obj))
}

print(period_changes)
print(period_interactions)
print(breakpoint_interactions_df)
