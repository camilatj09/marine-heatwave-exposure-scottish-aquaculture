# ==============================================================================
# 05 - PCA and hierarchical clustering of aquaculture-associated OSTIA cells
# ==============================================================================
# Purpose:
#   Characterise long-term MHW exposure across the 218 independent OSTIA cells
#   using five event metrics, PCA, and Ward.D2 hierarchical clustering.
#
# Inputs:
#   data/derived/cell_mhw_events/
#     scotland_OSTIA_cell_mhw_events_1990_2025.rds
#   data/derived/aquaculture_cells/
#     scotland_farm_OSTIA_cell_lookup.rds
#
# Main outputs:
#   data/derived/cell_exposure/
#     OSTIA_cell_long_term_MHW_metrics.csv
#     OSTIA_cell_PCA_scores_clusters.csv
#     OSTIA_cell_PCA_loadings.csv
#     OSTIA_cell_cluster_MHW_profiles.csv
#     scotland_OSTIA_cell_mhw_events_clustered.rds
#   outputs/figures/
#     figure_PCA_clusters_and_map.png
#     supplementary_cluster_diagnostics.png
# ==============================================================================

library(dplyr)
library(ggplot2)
library(cluster)
library(sf)
library(rnaturalearth)
library(patchwork)

# ------------------------------------------------------------------------------
# File paths
# ------------------------------------------------------------------------------

events_file <- file.path(
  "data", "derived", "cell_mhw_events",
  "scotland_OSTIA_cell_mhw_events_1990_2025.rds"
)

farm_lookup_file <- file.path(
  "data", "derived", "aquaculture_cells",
  "scotland_farm_OSTIA_cell_lookup.rds"
)

output_dir <- file.path("data", "derived", "cell_exposure")
figure_dir <- file.path("outputs", "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

metrics_file <- file.path(output_dir, "OSTIA_cell_long_term_MHW_metrics.csv")
scores_file <- file.path(output_dir, "OSTIA_cell_PCA_scores_clusters.csv")
loadings_file <- file.path(output_dir, "OSTIA_cell_PCA_loadings.csv")
variance_file <- file.path(output_dir, "OSTIA_cell_PCA_variance.csv")
profiles_file <- file.path(output_dir, "OSTIA_cell_cluster_MHW_profiles.csv")
sizes_file <- file.path(output_dir, "OSTIA_cell_cluster_sizes.csv")
diagnostics_file <- file.path(output_dir, "OSTIA_cell_cluster_diagnostics.csv")

pca_model_file <- file.path(output_dir, "OSTIA_cell_PCA_model.rds")
clustering_file <- file.path(output_dir, "OSTIA_cell_hierarchical_clustering.rds")
clustered_events_file <- file.path(
  output_dir, "scotland_OSTIA_cell_mhw_events_clustered.rds"
)

main_figure_file <- file.path(
  figure_dir, "figure_PCA_clusters_and_map.png"
)

diagnostic_figure_file <- file.path(
  figure_dir, "supplementary_cluster_diagnostics.png"
)

# ------------------------------------------------------------------------------
# Long-term MHW metrics
# ------------------------------------------------------------------------------

cell_events <- readRDS(events_file)
farm_lookup <- readRDS(farm_lookup_file)

stopifnot(n_distinct(cell_events$ostia_cell) == 218)

pca_variables <- c(
  "n_events",
  "mean_duration",
  "mean_intensity",
  "mean_max_intensity",
  "mean_cumulative_intensity"
)

cell_metrics <- cell_events %>%
  group_by(ostia_cell, ostia_cell_lon, ostia_cell_lat) %>%
  summarise(
    n_farms_in_cell = first(n_farms_in_cell),
    n_events = n(),
    mean_duration = mean(duration, na.rm = TRUE),
    mean_intensity = mean(intensity_mean, na.rm = TRUE),
    mean_max_intensity = mean(intensity_max, na.rm = TRUE),
    mean_cumulative_intensity = mean(intensity_cumulative, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(ostia_cell)

stopifnot(
  nrow(cell_metrics) == 218,
  all(complete.cases(cell_metrics[pca_variables]))
)

# ------------------------------------------------------------------------------
# PCA
# ------------------------------------------------------------------------------

pca_matrix <- as.matrix(cell_metrics[pca_variables])
rownames(pca_matrix) <- cell_metrics$ostia_cell

cell_pca <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)

variance <- cell_pca$sdev^2 / sum(cell_pca$sdev^2)

variance_table <- data.frame(
  PC = paste0("PC", seq_along(variance)),
  variance_explained_pct = 100 * variance,
  cumulative_variance_pct = 100 * cumsum(variance)
)

loadings <- data.frame(
  metric = rownames(cell_pca$rotation),
  cell_pca$rotation,
  row.names = NULL
) %>%
  mutate(
    PC1_contribution_pct = 100 * PC1^2,
    PC2_contribution_pct = 100 * PC2^2
  )

scores <- data.frame(
  ostia_cell = as.numeric(rownames(cell_pca$x)),
  cell_pca$x,
  row.names = NULL
) %>%
  left_join(cell_metrics, by = "ostia_cell")

stopifnot(sum(variance[1:2]) > 0.95)

# ------------------------------------------------------------------------------
# Ward.D2 clustering on PC1-PC2
# ------------------------------------------------------------------------------

cluster_data <- scores[, c("PC1", "PC2")]
cell_distance <- dist(cluster_data, method = "euclidean")
cell_tree <- hclust(cell_distance, method = "ward.D2")

scores$cluster_raw <- cutree(cell_tree, k = 2)

# Standardise labels: Cluster 1 is the profile with higher long-term intensity.
cluster_labels <- scores %>%
  group_by(cluster_raw) %>%
  summarise(
    mean_intensity = mean(mean_intensity),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_intensity)) %>%
  mutate(cluster = paste("Cluster", row_number()))

scores <- scores %>%
  left_join(cluster_labels, by = "cluster_raw") %>%
  mutate(cluster = factor(cluster, levels = c("Cluster 1", "Cluster 2")))

cluster_sizes <- scores %>%
  count(cluster, name = "n_OSTIA_cells")

stopifnot(
  sum(cluster_sizes$n_OSTIA_cells) == 218,
  setequal(cluster_sizes$n_OSTIA_cells, c(173L, 45L))
)

cluster_profiles <- scores %>%
  group_by(cluster) %>%
  summarise(
    n_OSTIA_cells = n(),
    mean_n_events = mean(n_events),
    mean_duration_days = mean(mean_duration),
    mean_event_intensity_C = mean(mean_intensity),
    mean_max_intensity_C = mean(mean_max_intensity),
    mean_cumulative_intensity_Cdays = mean(mean_cumulative_intensity),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Cluster-selection diagnostics
# ------------------------------------------------------------------------------

within_cluster_ss <- function(k) {

  groups <- if (k == 1) {
    rep(1L, nrow(cluster_data))
  } else {
    cutree(cell_tree, k = k)
  }

  sum(vapply(
    split(seq_len(nrow(cluster_data)), groups),
    function(i) {
      x <- as.matrix(cluster_data[i, , drop = FALSE])
      centre <- colMeans(x)
      sum((x - rep(centre, each = nrow(x)))^2)
    },
    numeric(1)
  ))
}

diagnostic_k <- 1:8

diagnostics <- data.frame(
  k = diagnostic_k,
  WSS = vapply(diagnostic_k, within_cluster_ss, numeric(1)),
  silhouette_width = NA_real_
)

for (k in 2:8) {
  cl <- cutree(cell_tree, k = k)
  diagnostics$silhouette_width[diagnostics$k == k] <-
    mean(cluster::silhouette(cl, cell_distance)[, "sil_width"])
}

selected_silhouette <- diagnostics$silhouette_width[diagnostics$k == 2]

# ------------------------------------------------------------------------------
# Attach cluster assignments to events
# ------------------------------------------------------------------------------

cell_cluster_lookup <- scores %>%
  select(ostia_cell, cluster)

cell_events_clustered <- cell_events %>%
  left_join(
    cell_cluster_lookup,
    by = "ostia_cell",
    relationship = "many-to-one"
  ) %>%
  arrange(ostia_cell, date_start)

stopifnot(!anyNA(cell_events_clustered$cluster))

# ------------------------------------------------------------------------------
# PCA space and mapped farm locations
# ------------------------------------------------------------------------------

cluster_colours <- c("Cluster 1" = "#2C7BB6", "Cluster 2" = "#D7191C")

pc1_label <- sprintf("PC1 (%.1f%%)", 100 * variance[1])
pc2_label <- sprintf("PC2 (%.1f%%)", 100 * variance[2])

arrow_scale <- 0.35 * min(
  diff(range(scores$PC1)),
  diff(range(scores$PC2))
)

loading_plot <- loadings %>%
  filter(metric %in% pca_variables) %>%
  mutate(
    xend = PC1 * arrow_scale,
    yend = PC2 * arrow_scale,
    label = c(
      "Events", "Duration", "Mean intensity",
      "Maximum intensity", "Cumulative intensity"
    )
  )

p_pca <- ggplot(scores, aes(PC1, PC2, colour = cluster)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey75") +
  geom_point(size = 2, alpha = 0.8) +
  geom_segment(
    data = loading_plot,
    aes(x = 0, y = 0, xend = xend, yend = yend),
    inherit.aes = FALSE,
    arrow = grid::arrow(length = grid::unit(0.15, "cm")),
    linewidth = 0.5
  ) +
  geom_text(
    data = loading_plot,
    aes(x = xend, y = yend, label = label),
    inherit.aes = FALSE,
    size = 3,
    nudge_y = 0.12
  ) +
  scale_colour_manual(values = cluster_colours) +
  labs(x = pc1_label, y = pc2_label, colour = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top")

farm_clusters <- farm_lookup %>%
  distinct(site_id, longitude, latitude, ostia_cell) %>%
  left_join(
    cell_cluster_lookup,
    by = "ostia_cell",
    relationship = "many-to-one"
  )

stopifnot(
  n_distinct(farm_clusters$site_id) == 361,
  !anyNA(farm_clusters$cluster)
)

farm_sf <- st_as_sf(
  farm_clusters,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

uk <- rnaturalearth::ne_countries(
  country = "United Kingdom",
  scale = "medium",
  returnclass = "sf"
)

p_map <- ggplot() +
  geom_sf(data = uk, fill = "grey96", colour = "grey45", linewidth = 0.3) +
  geom_sf(
    data = farm_sf,
    aes(colour = cluster),
    size = 1.2,
    alpha = 0.85
  ) +
  scale_colour_manual(values = cluster_colours) +
  coord_sf(xlim = c(-8.6, -0.5), ylim = c(55.1, 60.9), expand = FALSE) +
  labs(x = NULL, y = NULL, colour = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top")

main_figure <- p_pca + p_map +
  patchwork::plot_annotation(tag_levels = "A")

ggsave(
  main_figure_file,
  main_figure,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

p_wss <- ggplot(diagnostics, aes(k, WSS)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = 2, linetype = "dashed") +
  scale_x_continuous(breaks = diagnostic_k) +
  labs(x = "Number of clusters (k)", y = "Within-cluster sum of squares") +
  theme_classic(base_size = 11)

p_sil <- diagnostics %>%
  filter(k >= 2) %>%
  ggplot(aes(k, silhouette_width)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = 2, linetype = "dashed") +
  scale_x_continuous(breaks = 2:8) +
  labs(x = "Number of clusters (k)", y = "Average silhouette width") +
  theme_classic(base_size = 11)

diagnostic_figure <- p_wss + p_sil +
  patchwork::plot_annotation(tag_levels = "A")

ggsave(
  diagnostic_figure_file,
  diagnostic_figure,
  width = 8,
  height = 3.8,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------

write.csv(cell_metrics, metrics_file, row.names = FALSE)
write.csv(scores, scores_file, row.names = FALSE)
write.csv(loadings, loadings_file, row.names = FALSE)
write.csv(variance_table, variance_file, row.names = FALSE)
write.csv(cluster_profiles, profiles_file, row.names = FALSE)
write.csv(cluster_sizes, sizes_file, row.names = FALSE)
write.csv(diagnostics, diagnostics_file, row.names = FALSE)

saveRDS(cell_pca, pca_model_file)
saveRDS(cell_tree, clustering_file)
saveRDS(cell_events_clustered, clustered_events_file)

print(variance_table[1:2, ])
print(cluster_sizes)
print(cluster_profiles)
cat("Mean silhouette width (k = 2):", selected_silhouette, "\n")
