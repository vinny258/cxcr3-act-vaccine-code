
# Location-independent paths and fonts. See helpers/paths.R
source(file.path(rprojroot_find <- {
  d <- tryCatch(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])),
                error = function(e) getwd())
  if (!length(d) || !nzchar(d)) d <- getwd()
  while (!file.exists(file.path(d, "helpers", "paths.R")) && dirname(d) != d) d <- dirname(d)
  d
}, "helpers", "paths.R"))

# Write outputs beside this script, or wherever FIG_OUTDIR points.
setwd(out_dir())

# ============================================================
# Cluster analysis of Tet+ CD8+ T cells — Exp333 (ACT in LN)
# Loads pre-computed UMAP + expression data from:
#   UMAP_Exp333_Tet_combined_df.rds
# Adds k-means clustering (k = 5) and generates:
#   Figure 3H  — UMAP coloured by cluster (overall and by group)
#   Figure 3I  — Density plots of arcsinh marker expression per cluster
#   Figure S5A — Heatmap: median marker expression per cluster (two
#                orientations, both with z-score colour bar label):
#                  UMAP_Exp333_Tet_cluster_heatmap_vertical.pdf
#                    markers as rows, clusters as columns (portrait 4×7")
#                  UMAP_Exp333_Tet_cluster_heatmap_horizontal.pdf
#                    clusters as rows, markers as columns (landscape 7×4")
#                Colour scale: column z-score of arcsinh-transformed
#                median expression. Label added via grid::grid.text.
#                Saved with pdf() + grid::grid.draw() (not ggsave).
#
# Last updated: 2026-04-27
# ============================================================

# ---- 0. Packages -------------------------------------------
pkgs_cran <- c("ggplot2", "dplyr", "tidyr", "patchwork",
               "RColorBrewer", "viridis", "scales", "pheatmap")
lapply(pkgs_cran, library, character.only = TRUE)

# ---- 1. Paths ----------------------------------------------
analysis_dir <- out_dir()
output_dir <- out_dir()

# ---- 2. Load data ------------------------------------------
cat("Loading combined_df...\n")
combined_df <- readRDS(file.path(data_dir("UMAP_Exp333_Tet_combined_df.rds"), "UMAP_Exp333_Tet_combined_df.rds"))

cat("Dimensions:", nrow(combined_df), "cells x", ncol(combined_df), "columns\n")
cat("Groups:", paste(unique(combined_df$Group), collapse = ", "), "\n")
cat("LN types:", paste(unique(combined_df$LN_type), collapse = ", "), "\n")
cat("Columns:\n"); print(colnames(combined_df))

# ---- 3. Identify arcsinh marker columns --------------------
asin_cols <- grep("_asin$", colnames(combined_df), value = TRUE)
marker_names <- sub("_asin$", "", asin_cols)

cat("\nMarkers for clustering:", paste(marker_names, collapse = ", "), "\n")

# ---- 4. K-means clustering ---------------------------------
# Try k = 4, 5, 6 and pick by visual inspection; default = 5
K <- 5
set.seed(42)

expr_mat <- as.matrix(combined_df[, asin_cols])

cat("Running k-means with k =", K, "...\n")
km <- kmeans(expr_mat, centers = K, nstart = 25, iter.max = 100)

combined_df$cluster <- factor(paste0("C", km$cluster))
cat("Cluster sizes:\n"); print(table(combined_df$cluster))

# ---- 5. Cluster labels -------------------------------------
cluster_labels <- c(
  C1 = "Effector (CD62Lmod)",
  C2 = "Naive-like",
  C3 = "Effector (CD62Llo)",
  C4 = "Central Memory",
  C5 = "Central Memory (CXCR3hi)"
)

cluster_order <- c("Naive-like", "Central Memory", "Central Memory (CXCR3hi)",
                   "Effector (CD62Lmod)", "Effector (CD62Llo)")

combined_df$cluster_name <- factor(
  cluster_labels[as.character(combined_df$cluster)],
  levels = cluster_order
)

# ---- 6. Colour palette (Okabe-Ito, colorblind-friendly) ----
# Okabe-Ito palette: safe for deuteranopia, protanopia, tritanopia
okabe_ito_5 <- c(
  "Naive-like"                  = "#56B4E9",   # sky blue
  "Central Memory"              = "#009E73",   # bluish green
  "Central Memory (CXCR3hi)"   = "#CC79A7",   # reddish purple
  "Effector (CD62Lmod)"        = "#D55E00",   # vermillion
  "Effector (CD62Llo)"         = "#E69F00"    # orange
)

cluster_colours <- okabe_ito_5

group_colours <- setNames(
  brewer.pal(max(3, length(unique(na.omit(combined_df$Group)))), "Set1")[
    seq_along(unique(na.omit(combined_df$Group)))],
  sort(unique(na.omit(combined_df$Group)))
)

# ---- 7. Shared UMAP theme ----------------------------------
theme_umap <- theme_classic(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    axis.line       = element_line(colour = "grey70"),
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 10),
    aspect.ratio    = 1
  )

# ---- 8. UMAP: coloured by cluster --------------------------
p_cluster <- ggplot(combined_df,
                    aes(x = UMAP1, y = UMAP2, colour = cluster_name)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_colour_manual(values = cluster_colours, name = "Cluster") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = paste0("Tet+ CD8+ T cells — k-means clusters (k=", K, ")"),
       x = "UMAP 1", y = "UMAP 2") +
  theme_umap
print(p_cluster)

# ---- 9. UMAP: faceted by Group, coloured by cluster --------
p_cluster_by_group <- ggplot(
    combined_df %>% filter(!is.na(Group)),
    aes(x = UMAP1, y = UMAP2, colour = cluster_name)) +
  geom_point(size = 0.35, alpha = 0.6) +
  scale_colour_manual(values = cluster_colours, name = "Cluster") +
  facet_wrap(~ Group) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "Clusters by group",
       x = "UMAP 1", y = "UMAP 2") +
  theme_umap +
  theme(strip.text       = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey90", colour = NA))
print(p_cluster_by_group)

# ---- 10. UMAP: faceted by Group × LN type ------------------
p_cluster_by_group_ln <- ggplot(
    combined_df %>% filter(!is.na(Group)),
    aes(x = UMAP1, y = UMAP2, colour = cluster_name)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_colour_manual(values = cluster_colours, name = "Cluster") +
  facet_grid(LN_type ~ Group) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "Clusters by group and LN type",
       x = "UMAP 1", y = "UMAP 2") +
  theme_umap +
  theme(strip.text       = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey90", colour = NA),
        aspect.ratio     = 0.9)
print(p_cluster_by_group_ln)

# ---- 11. Stacked bar: cluster proportions per group --------
# Calculate per-sample proportions first, then average per group
# (avoids larger samples dominating)
prop_by_sample <- combined_df %>%
  filter(!is.na(Group)) %>%
  count(file_name, Group, LN_type, cluster_name) %>%
  group_by(file_name) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Mean proportion per group (across samples)
prop_by_group <- prop_by_sample %>%
  group_by(Group, cluster_name) %>%
  summarise(mean_prop = mean(prop), se_prop = sd(prop) / sqrt(n()), .groups = "drop")

p_bar_group <- ggplot(prop_by_group,
                      aes(x = Group, y = mean_prop, fill = cluster_name)) +
  geom_col(colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = cluster_colours, name = "Cluster") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Cluster proportions per group",
       x = NULL, y = "Mean proportion per sample") +
  theme_classic(base_size = 13) +
  theme(
    plot.title  = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.title = element_text(size = 11)
  )
print(p_bar_group)

# ---- 12. Stacked bar: cluster proportions per group × LN type --
prop_by_group_ln <- prop_by_sample %>%
  group_by(Group, LN_type, cluster_name) %>%
  summarise(mean_prop = mean(prop), .groups = "drop")

p_bar_group_ln <- ggplot(prop_by_group_ln,
                          aes(x = Group, y = mean_prop, fill = cluster_name)) +
  geom_col(colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = cluster_colours, name = "Cluster") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  facet_wrap(~ LN_type) +
  labs(title = "Cluster proportions per group × LN type",
       x = NULL, y = "Mean proportion per sample") +
  theme_classic(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 30, hjust = 1),
    strip.text       = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90", colour = NA),
    legend.title     = element_text(size = 11)
  )
print(p_bar_group_ln)

# ---- 13. Dot plot: cluster × group proportions -------------
p_dot_group <- ggplot(prop_by_group,
                       aes(x = Group, y = cluster_name,
                           size = mean_prop, colour = cluster_name)) +
  geom_point() +
  scale_size_continuous(range = c(2, 12), labels = percent_format(accuracy = 1),
                        name = "Mean %") +
  scale_colour_manual(values = cluster_colours, guide = "none") +
  labs(title = "Cluster frequency per group",
       x = NULL, y = "Cluster") +
  theme_classic(base_size = 13) +
  theme(
    plot.title  = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.y = element_line(colour = "grey90")
  )
print(p_dot_group)

# ---- 14. Heatmap: median marker expression per cluster -----
cluster_medians <- combined_df %>%
  group_by(cluster_name) %>%
  summarise(across(all_of(asin_cols), median, .names = "{.col}"),
            .groups = "drop")

heatmap_mat <- as.matrix(cluster_medians[, asin_cols])
rownames(heatmap_mat) <- cluster_medians$cluster_name
colnames(heatmap_mat) <- marker_names

# Scale per marker (column z-score) for easier pattern reading
heatmap_scaled <- scale(heatmap_mat)

ann_row <- data.frame(Cluster = rownames(heatmap_mat),
                      row.names = rownames(heatmap_mat))
ann_colours <- list(Cluster = cluster_colours)

# Vertical heatmap: markers as rows, clusters as columns
p_heatmap_v <- pheatmap(t(heatmap_scaled),
         cluster_rows      = FALSE,
         cluster_cols      = TRUE,
         color             = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
         border_color      = NA,
         cellwidth         = 35,
         cellheight        = 50,
         fontsize          = 12,
         annotation_col    = ann_row,
         annotation_colors = ann_colours,
         main              = paste0("Median marker expression per cluster (z-score, k=", K, ")"),
         silent            = TRUE)

# Horizontal heatmap: clusters as rows, markers as columns
p_heatmap_h <- pheatmap(heatmap_scaled,
         cluster_rows      = TRUE,
         cluster_cols      = FALSE,
         color             = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
         border_color      = NA,
         cellwidth         = 50,
         cellheight        = 35,
         fontsize          = 12,
         annotation_row    = ann_row,
         annotation_colors = ann_colours,
         main              = paste0("Median marker expression per cluster (z-score, k=", K, ")"),
         silent            = TRUE)

cat("Saving heatmaps...\n")
dir.create(file.path(output_dir, "FigS5A"), showWarnings = FALSE)
pdf(file.path(output_dir, "FigS5A", "UMAP_Exp333_Tet_cluster_heatmap_vertical.pdf"), width = 4, height = 7)
grid::grid.newpage()
grid::grid.draw(p_heatmap_v$gtable)
grid::grid.text("z-score", x = 0.97, y = 0.52, rot = 270,
                gp = grid::gpar(fontsize = 9))
dev.off()
cat("  Saved FigS5A/UMAP_Exp333_Tet_cluster_heatmap_vertical.pdf\n")

pdf(file.path(output_dir, "FigS5A", "UMAP_Exp333_Tet_cluster_heatmap_horizontal.pdf"), width = 7, height = 4)
grid::grid.newpage()
grid::grid.draw(p_heatmap_h$gtable)
grid::grid.text("z-score", x = 0.97, y = 0.52, rot = 270,
                gp = grid::gpar(fontsize = 9))
dev.off()
cat("  Saved FigS5A/UMAP_Exp333_Tet_cluster_heatmap_horizontal.pdf\n")

# ---- 16. Save PDFs -----------------------------------------
dir.create(file.path(output_dir, "Fig3H"), showWarnings = FALSE)
ggsave(file.path(output_dir, "Fig3H", "UMAP_Exp333_Tet_clusters.pdf"),
       plot = p_cluster, width = 7, height = 6)

ggsave(file.path(output_dir, "Fig3H", "UMAP_Exp333_Tet_clusters_by_group.pdf"),
       plot = p_cluster_by_group,
       width = 4 * length(unique(na.omit(combined_df$Group))),
       height = 5)

ggsave(file.path(output_dir, "Fig3H", "UMAP_Exp333_Tet_clusters_by_group_LN.pdf"),
       plot = p_cluster_by_group_ln,
       width = 4 * length(unique(na.omit(combined_df$Group))),
       height = 9)

ggsave("UMAP_Exp333_Tet_cluster_bar_group.pdf",
       plot = p_bar_group, width = 6, height = 5)

ggsave("UMAP_Exp333_Tet_cluster_bar_group_LN.pdf",
       plot = p_bar_group_ln, width = 10, height = 5)

ggsave("UMAP_Exp333_Tet_cluster_dot_group.pdf",
       plot = p_dot_group, width = 6, height = 5)

cat("\nAll plots saved:\n")
cat("  UMAP_Exp333_Tet_clusters.pdf\n")
cat("  UMAP_Exp333_Tet_clusters_by_group.pdf\n")
cat("  UMAP_Exp333_Tet_clusters_by_group_LN.pdf\n")
cat("  UMAP_Exp333_Tet_cluster_bar_group.pdf\n")
cat("  UMAP_Exp333_Tet_cluster_bar_group_LN.pdf\n")
cat("  UMAP_Exp333_Tet_cluster_dot_group.pdf\n")
cat("  UMAP_Exp333_Tet_cluster_heatmap.pdf\n")

# ---- 17. Marker density UMAPs (hexbin smoothed) ------------
# Each plot bins the UMAP into a hex grid and colours by median
# arcsinh-transformed expression — gives a smooth density map
# per marker without individual dot overplotting.

theme_density <- theme_void(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 12,
                                    hjust = 0.5, margin = margin(b = 4)),
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 9),
    legend.key.width = unit(0.5, "cm"),
    strip.text       = element_text(face = "bold"),
    aspect.ratio     = 1
  )

density_plots <- lapply(seq_along(asin_cols), function(i) {
  mc <- asin_cols[i]
  mn <- marker_names[i]

  ggplot(combined_df, aes(x = UMAP1, y = UMAP2, z = .data[[mc]])) +
    stat_summary_hex(fun = median, bins = 60) +
    scale_fill_gradientn(
      colours = c("#0D0221", "#3B0F70", "#8C2981",
                  "#DE4968", "#FE9F6D", "#FCFDBF"),
      name    = "Median\n(arcsinh)",
      limits  = c(0, max(combined_df[[mc]], na.rm = TRUE))
    ) +
    labs(title = mn) +
    theme_density
})

p_density_panel <- wrap_plots(density_plots, ncol = 3) +
  plot_annotation(
    title = "Tet+ CD8+ T cells — marker expression density",
    theme = theme(plot.title = element_text(face = "bold", size = 14,
                                             hjust = 0.5))
  )

print(p_density_panel)

dir.create(file.path(output_dir, "Fig3I"), showWarnings = FALSE)
ggsave(file.path(output_dir, "Fig3I", "UMAP_Exp333_Tet_marker_density.pdf"),
       plot = p_density_panel, width = 14, height = 9)
cat("Saved Fig3I/UMAP_Exp333_Tet_marker_density.pdf\n")

# ---- 18. Overall cell density UMAP -------------------------
# Shows where cells are most concentrated regardless of marker
p_cell_density <- ggplot(combined_df, aes(x = UMAP1, y = UMAP2)) +
  stat_density_2d_filled(contour_var = "ndensity",
                         bins = 12, alpha = 0.85) +
  scale_fill_viridis_d(option = "magma", name = "Density") +
  geom_point(size = 0.15, alpha = 0.15, colour = "white") +
  labs(title = "Overall cell density", x = "UMAP 1", y = "UMAP 2") +
  theme_umap

print(p_cell_density)

ggsave("UMAP_Exp333_Tet_cell_density.pdf",
       plot = p_cell_density, width = 7, height = 6)
cat("Saved UMAP_Exp333_Tet_cell_density.pdf\n")

# ---- 19. Save updated combined_df with cluster labels ------
saveRDS(combined_df, "UMAP_Exp333_Tet_combined_df.rds")
cat("\nUpdated combined_df (with cluster column) saved.\n")

message("Done!")
