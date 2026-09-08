
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
# Statistical comparison of cluster proportions across groups
# Exp333 — Tet+ CD8+ T cells (ACT in LN)
#
# Design:
#   3 groups  × 2 LN types × 5 clusters
#   Unit of replication: individual mouse (file_name)
#   Test: Kruskal-Wallis per cluster × LN type,
#         post-hoc pairwise Wilcoxon (BH correction)
# ============================================================

# ---- 0. Packages -------------------------------------------
pkgs <- c("ggplot2", "dplyr", "tidyr", "rstatix",
          "ggpubr", "scales", "writexl")

for (p in pkgs)

lapply(pkgs, library, character.only = TRUE)

# ---- 1. Paths ----------------------------------------------
analysis_dir <- out_dir()

# ---- 2. Load data ------------------------------------------
cat("Loading combined_df...\n")
combined_df <- readRDS(file.path(data_dir("UMAP_Exp333_Tet_combined_df.rds"), "UMAP_Exp333_Tet_combined_df.rds"))
cat("Dimensions:", nrow(combined_df), "cells x", ncol(combined_df), "columns\n")

# ---- 3. Per-sample proportions per cluster × LN type -------
# Each row = one mouse × LN type × cluster combination
prop_df <- combined_df %>%
  filter(!is.na(Group), !is.na(LN_type)) %>%
  count(file_name, Group, LN_type, cluster_name, .drop = FALSE) %>%
  group_by(file_name, LN_type) %>%
  mutate(prop = n / sum(n) * 100) %>%   # express as %
  ungroup()

# Fix factor ordering
group_order <- c("T cells", "T cells + ChAd-i.a.", "T cells + ChAd-P1A")
prop_df$Group <- factor(prop_df$Group, levels = group_order)

cluster_order <- c("Naive-like", "Central Memory", "Central Memory (CXCR3hi)",
                   "Effector (CD62Lmod)", "Effector (CD62Llo)")
prop_df$cluster_name <- factor(prop_df$cluster_name, levels = cluster_order)

cat("\nSamples per group × LN type:\n")
print(
  prop_df %>%
    distinct(file_name, Group, LN_type) %>%
    count(Group, LN_type)
)

# ---- 4. Kruskal-Wallis + pairwise Wilcoxon -----------------
# Run for each cluster × LN type stratum
kw_results <- prop_df %>%
  group_by(cluster_name, LN_type) %>%
  kruskal_test(prop ~ Group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

cat("\nKruskal-Wallis results (BH-adjusted across all strata):\n")
print(kw_results %>% select(cluster_name, LN_type, statistic, df, p, p.adj, p.adj.signif))

# Pairwise Wilcoxon for all cluster × LN type combinations
pwc <- prop_df %>%
  group_by(cluster_name, LN_type) %>%
  pairwise_wilcox_test(prop ~ Group, p.adjust.method = "BH") %>%
  add_significance("p.adj")

cat("\nPairwise Wilcoxon results:\n")
print(pwc %>% select(cluster_name, LN_type, group1, group2, n1, n2, p, p.adj, p.adj.signif))

# ---- 5. Colour palette (Okabe-Ito, matches existing plots) --
okabe_ito_5 <- c(
  "Naive-like"                 = "#56B4E9",
  "Central Memory"             = "#009E73",
  "Central Memory (CXCR3hi)"  = "#CC79A7",
  "Effector (CD62Lmod)"       = "#D55E00",
  "Effector (CD62Llo)"        = "#E69F00"
)

group_shapes <- c("T cells" = 16, "T cells + ChAd-i.a." = 17, "T cells + ChAd-P1A" = 15)

# ---- 6. Plot: jitter + box per cluster, faceted by LN type --
# Attach pairwise stats and filter to significant comparisons only
pwc_signif <- pwc %>%
  add_xy_position(x = "Group", dodge = 0.8) %>%
  filter(p.adj <= 0.05)

p_stats <- ggplot(prop_df,
                  aes(x = Group, y = prop, colour = cluster_name)) +
  geom_boxplot(aes(group = Group), outlier.shape = NA,
               colour = "grey40", fill = NA, linewidth = 0.5) +
  geom_jitter(aes(shape = Group), width = 0.15, size = 2.5, alpha = 0.8) +
  scale_colour_manual(values = okabe_ito_5, name = "Cluster") +
  scale_shape_manual(values = group_shapes, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0.02, 0.15))) +
  facet_grid(cluster_name ~ LN_type, scales = "free_y") +
  labs(
    title = "Cluster proportions per group and LN type",
    subtitle = "Kruskal-Wallis + pairwise Wilcoxon (BH); brackets = p.adj ≤ 0.05",
    x = NULL, y = "% of Tet+ CD8+ T cells"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    axis.text.x      = element_text(angle = 35, hjust = 1),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey92", colour = NA),
    legend.position  = "none"
  )

# Add significance brackets only if there are significant pairs
if (nrow(pwc_signif) > 0) {
  p_stats <- p_stats +
    stat_pvalue_manual(
      pwc_signif,
      label        = "p.adj.signif",
      hide.ns      = TRUE,
      tip.length   = 0.01,
      size         = 3.5
    )
}

print(p_stats)

ggsave("UMAP_Exp333_Tet_cluster_stats_jitter.pdf",
       plot = p_stats,
       width = 3 * length(unique(prop_df$LN_type)),
       height = 2.8 * length(cluster_order))

cat("Saved UMAP_Exp333_Tet_cluster_stats_jitter.pdf\n")

# ---- 7. Alternative: one panel per LN type, clusters on x axis ------
p_overview <- ggplot(prop_df,
                     aes(x = cluster_name, y = prop,
                         colour = Group, group = Group)) +
  geom_boxplot(aes(fill = Group), alpha = 0.15,
               outlier.shape = NA, position = position_dodge(0.7),
               linewidth = 0.5) +
  geom_point(aes(shape = Group),
             position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width  = 0.7),
             size = 2.5, alpha = 0.9) +
  scale_colour_brewer(palette = "Set1", name = "Group") +
  scale_fill_brewer(palette = "Set1", name = "Group") +
  scale_shape_manual(values = group_shapes, name = "Group") +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0.02, 0.15))) +
  facet_wrap(~ LN_type, ncol = 1) +
  labs(
    title    = "Cluster proportions — overview",
    subtitle = "Kruskal-Wallis + pairwise Wilcoxon (BH)",
    x = "Cluster", y = "% of Tet+ CD8+ T cells"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    axis.text.x      = element_text(angle = 35, hjust = 1),
    strip.text       = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92", colour = NA)
  )

# Add significance brackets if any
pwc_overview <- pwc %>%
  add_xy_position(x = "cluster_name", dodge = 0.7) %>%
  filter(p.adj <= 0.05)

if (nrow(pwc_overview) > 0) {
  p_overview <- p_overview +
    stat_pvalue_manual(
      pwc_overview,
      label        = "p.adj.signif",
      hide.ns      = TRUE,
      tip.length   = 0.01,
      size         = 3
    )
}

print(p_overview)

ggsave("UMAP_Exp333_Tet_cluster_stats_overview.pdf",
       plot = p_overview, width = 10, height = 9)
cat("Saved UMAP_Exp333_Tet_cluster_stats_overview.pdf\n")

# ---- 8. Export summary tables to Excel ----------------------
write_xlsx(
  list(
    "Kruskal-Wallis"       = as.data.frame(kw_results),
    "Pairwise Wilcoxon"    = as.data.frame(pwc),
    "Per-sample proportions" = as.data.frame(prop_df)
  ),
  path = "UMAP_Exp333_Tet_cluster_stats.xlsx"
)

cat("Saved UMAP_Exp333_Tet_cluster_stats.xlsx\n")

# ---- 9. Print summary to console ----------------------------
cat("\n============================================================\n")
cat("SUMMARY — significant pairwise differences (p.adj ≤ 0.05)\n")
cat("============================================================\n")
sig <- pwc %>%
  filter(p.adj <= 0.05) %>%
  select(cluster_name, LN_type, group1, group2, n1, n2, p, p.adj, p.adj.signif) %>%
  arrange(cluster_name, LN_type)

if (nrow(sig) == 0) {
  cat("No significant differences after BH correction.\n")
} else {
  print(sig, n = Inf)
}

message("\nDone!")
