
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
# Cluster analysis — absolute cell numbers per LN
# Exp333 — Tet+ CD8+ T cells (ACT in LN)
#
# Strategy:
#   per-cluster absolute count = (cluster % within Tet+)
#                                × (total Tet+ E+3 cells/LN)
#
# Tests: Kruskal-Wallis + pairwise Wilcoxon (BH)
#
# Outputs:
#   Figure S5C — UMAP_Exp333_Tet_cluster_abs_jitter.pdf
#     Absolute counts per LN (log10 scale). Clusters as rows,
#     LN type as columns. facet_wrap(scales = "free_y") gives
#     each panel an independent y-axis. Significance brackets
#     drawn with geom_segment + geom_text; y-positions computed
#     per panel in log10 space (max_y × 10^(0.18 × rank)).
#     No title/subtitle (manuscript-ready).
#
#   Figure S5B — UMAP_Exp333_Tet_cluster_prop_jitter.pdf
#     Cluster proportions (% of Tet+ CD8+), linear scale.
#     facet_wrap(scales = "free_y") per panel. Significance
#     brackets via stat_pvalue_manual (ggpubr).
#     No title/subtitle (manuscript-ready).
#
#   UMAP_Exp333_Tet_cluster_abs_stats.xlsx
#     Kruskal-Wallis, pairwise Wilcoxon, and per-sample
#     absolute counts exported to Excel.
#
# Last updated: 2026-04-27
# ============================================================

# ---- 0. Packages -------------------------------------------
pkgs <- c("ggplot2", "dplyr", "tidyr", "rstatix",
          "ggpubr", "scales", "writexl", "readxl")

for (p in pkgs)

lapply(pkgs, library, character.only = TRUE)

# ---- 1. Paths ----------------------------------------------
analysis_dir <- out_dir()
output_dir <- out_dir()

# ---- 2. Load cluster data from UMAP RDS --------------------
cat("Loading combined_df...\n")
combined_df <- readRDS(file.path(data_dir("UMAP_Exp333_Tet_combined_df.rds"), "UMAP_Exp333_Tet_combined_df.rds"))
cat("Dimensions:", nrow(combined_df), "cells x", ncol(combined_df), "cols\n")

# ---- 3. Per-sample cluster proportions ---------------------
cluster_order <- c("Naive-like", "Central Memory", "Central Memory (CXCR3hi)",
                   "Effector (CD62Lmod)", "Effector (CD62Llo)")

# Strip the numeric suffix that FlowJo appends: "dLN_M_1.fcs_000..." → "dLN_M_1.fcs"
combined_df <- combined_df %>%
  mutate(file_name_short = sub("\\.fcs_.*$", ".fcs", file_name))

prop_df <- combined_df %>%
  filter(!is.na(Group), !is.na(LN_type)) %>%
  count(file_name_short, Group, LN_type, cluster_name, .drop = FALSE) %>%
  group_by(file_name_short, LN_type) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(cluster_name = factor(cluster_name, levels = cluster_order))

cat("Unique file names in RDS (stripped):\n")
print(sort(unique(prop_df$file_name_short)))

# ---- 4. Load Tet+ E+3 cells/LN from Excel -----------------
cat("\nLoading Exp333 Excel sheet...\n")
xls_path <- "20260409_Exp333_LN.xls"

raw_xls <- read_excel(xls_path, sheet = "Exp333", col_names = FALSE)

# The sheet has dLN (cols 1-30) and ndLN (cols 31-60) side by side.
# Col 0  = Group
# Col 1  = dLN sample name
# Col 41 = Tet+ CD8+CD3+ E+3 cells/LN  (right-hand ndLN block, 0-indexed = col 41)
# Col 31 = ndLN sample name
#
# In R (1-indexed):
#   col 1  = Group, col 2 = dLN sample, col 32 = ndLN sample
#   col 12 = Tet+ cells dLN (NOT E+3), col 42 = Tet+ E+3 cells ndLN
#   BUT the user says "Tet+ CD8+CD3+ E+3 cells/LN" — confirmed col 42 (1-indexed)

# Identify header row (row 1 in Excel = row 1 in readxl)
headers <- as.character(raw_xls[1, ])

# The Excel sheet has dLN (left) and ndLN (right) blocks side by side.
# dLN Tet+ column: "Tet+ CD8+CD3+ cells/LN"   (no E+3 in name, same units)
# ndLN Tet+ column: "Tet+ CD8+CD3+ E+3 cells/LN"
# Both are in E+3 (×1000 cells/LN) — only label differs.
tet_dln_col  <- which(grepl("Tet\\+.*CD8.*cells/LN", headers) &
                        !grepl("E\\+3", headers))
tet_ndln_col <- which(grepl("Tet\\+.*CD8.*E\\+3.*cells/LN", headers))

cat("dLN Tet+ column (1-indexed):", tet_dln_col,
    "→", headers[tet_dln_col], "\n")
cat("ndLN Tet+ column (1-indexed):", tet_ndln_col,
    "→", headers[tet_ndln_col], "\n")

# Data rows start at row 2
data_xls <- raw_xls[-1, ]

# Build a tidy sample → Tet+ cells/LN lookup
# dLN: col 2 = dLN sample name
dln_df <- data_xls %>%
  select(sample   = 2,
         group    = 1,
         tet_e3   = all_of(tet_dln_col)) %>%
  filter(!is.na(sample), sample != "") %>%
  mutate(LN_type = "dLN",
         tet_e3  = as.numeric(tet_e3))

# ndLN: col 32 = ndLN sample name
ndln_df <- data_xls %>%
  select(sample   = 32,
         group    = 1,
         tet_e3   = all_of(tet_ndln_col)) %>%
  filter(!is.na(sample), sample != "") %>%
  mutate(LN_type = "ndLN",
         tet_e3  = as.numeric(tet_e3))

tet_counts <- bind_rows(dln_df, ndln_df) %>%
  rename(file_name_short = sample) %>%
  select(file_name_short, LN_type, tet_e3)

cat("\nTet+ E+3 cells/LN lookup table:\n")
print(tet_counts, n = Inf)

# ---- 5. Join and compute absolute cluster counts -----------
abs_df <- prop_df %>%
  left_join(tet_counts, by = c("file_name_short", "LN_type")) %>%
  mutate(
    # absolute cells = proportion × total_Tet × 1000 (E+3 → actual cells)
    abs_cells = prop * tet_e3 * 1000,
    # pseudocount of 1 cell so log10(0) is avoided
    abs_cells_log = pmax(abs_cells, 1)
  )

# Sanity check
n_missing <- sum(is.na(abs_df$tet_e3))
if (n_missing > 0) {
  warning(n_missing, " rows could not be matched to Tet+ cell counts.")
  cat("Unmatched file_names:\n")
  print(
    abs_df %>%
      filter(is.na(tet_e3)) %>%
      distinct(file_name_short, LN_type)
  )
}

# Fix factor levels
group_order <- c("T cells", "T cells + ChAd-i.a.", "T cells + ChAd-P1A")
abs_df$Group <- factor(trimws(abs_df$Group), levels = group_order)

cat("\nTotal Tet+ cells/LN by group (raw values × 1000):\n")
print(
  tet_counts %>%
    left_join(
      combined_df %>%
        filter(!is.na(Group)) %>%
        distinct(file_name_short, Group, LN_type),
      by = c("file_name_short", "LN_type")
    ) %>%
    mutate(Group = factor(trimws(Group), levels = group_order)) %>%
    group_by(Group, LN_type) %>%
    summarise(
      n       = n(),
      median  = median(tet_e3 * 1000, na.rm = TRUE),
      min     = min(tet_e3 * 1000, na.rm = TRUE),
      max     = max(tet_e3 * 1000, na.rm = TRUE),
      .groups = "drop"
    )
)

# ---- 6. Statistical tests on absolute counts ---------------
kw_abs <- abs_df %>%
  filter(!is.na(abs_cells)) %>%
  group_by(cluster_name, LN_type) %>%
  kruskal_test(abs_cells ~ Group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

cat("\nKruskal-Wallis on ABSOLUTE COUNTS (BH-adjusted):\n")
print(kw_abs %>% select(cluster_name, LN_type, statistic, df, p, p.adj, p.adj.signif))

pwc_abs <- abs_df %>%
  filter(!is.na(abs_cells)) %>%
  group_by(cluster_name, LN_type) %>%
  pairwise_wilcox_test(abs_cells ~ Group, p.adjust.method = "BH") %>%
  add_significance("p.adj")

cat("\nPairwise Wilcoxon on ABSOLUTE COUNTS:\n")
print(pwc_abs %>% select(cluster_name, LN_type, group1, group2, n1, n2, p, p.adj, p.adj.signif))

# ---- 7. Colour palette -------------------------------------
okabe_ito_5 <- c(
  "Naive-like"                 = "#56B4E9",
  "Central Memory"             = "#009E73",
  "Central Memory (CXCR3hi)"  = "#CC79A7",
  "Effector (CD62Lmod)"       = "#D55E00",
  "Effector (CD62Llo)"        = "#E69F00"
)

group_shapes <- c(
  "T cells"              = 16,
  "T cells + ChAd-i.a." = 17,
  "T cells + ChAd-P1A"  = 15
)

# ---- 8. Jitter plot — absolute counts ----------------------
# Per-panel max → bracket y positions computed per panel in log space
max_per_panel <- abs_df %>%
  filter(!is.na(abs_cells_log)) %>%
  group_by(cluster_name, LN_type) %>%
  summarise(max_y = max(abs_cells_log, na.rm = TRUE), .groups = "drop")

# Build bracket annotation table with native ggplot2 coordinates:
#   x1 / x2  = discrete x positions (1 = T cells, 2 = ChAd-i.a., 3 = ChAd-P1A)
#   y.position = bracket top in linear cell-count space (log-transformed by scale)
#   y_tip      = bracket tip (0.06 log10 units below bracket top)
#   x_mid      = midpoint for significance label
pwc_abs_signif <- pwc_abs %>%
  filter(p.adj <= 0.05) %>%
  left_join(max_per_panel, by = c("cluster_name", "LN_type")) %>%
  group_by(cluster_name, LN_type) %>%
  arrange(cluster_name, LN_type, group1, group2) %>%
  mutate(
    y.position = max_y * 10^(0.18 * row_number()),
    y_tip      = y.position * 10^(-0.06),
    x1         = match(group1, group_order),
    x2         = match(group2, group_order),
    x_mid      = (x1 + x2) / 2
  ) %>%
  ungroup()

p_abs_jitter <- ggplot(
    abs_df %>% filter(!is.na(abs_cells_log)),
    aes(x = Group, y = abs_cells_log, colour = cluster_name)) +
  geom_boxplot(aes(group = Group), outlier.shape = NA, na.rm = TRUE,
               colour = "grey40", fill = NA, linewidth = 0.5) +
  geom_jitter(aes(shape = Group), width = 0.15, size = 2.5, alpha = 0.8) +
  scale_colour_manual(values = okabe_ito_5) +
  scale_shape_manual(values = group_shapes, guide = "none") +
  scale_y_log10(labels = comma,
                breaks = c(1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000),
                expand = expansion(add = c(0.05, 0.6))) +
  facet_wrap(cluster_name ~ LN_type, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "Cells per LN") +
  theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    axis.text.x      = element_text(angle = 35, hjust = 1),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey92", colour = NA),
    legend.position  = "none"
  )

if (nrow(pwc_abs_signif) > 0) {
  p_abs_jitter <- p_abs_jitter +
    geom_segment(data = pwc_abs_signif, inherit.aes = FALSE,
                 aes(x = x1, xend = x2,
                     y = y.position, yend = y.position),
                 colour = "black", linewidth = 0.4) +
    geom_segment(data = pwc_abs_signif, inherit.aes = FALSE,
                 aes(x = x1, xend = x1,
                     y = y_tip, yend = y.position),
                 colour = "black", linewidth = 0.4) +
    geom_segment(data = pwc_abs_signif, inherit.aes = FALSE,
                 aes(x = x2, xend = x2,
                     y = y_tip, yend = y.position),
                 colour = "black", linewidth = 0.4) +
    geom_text(data = pwc_abs_signif, inherit.aes = FALSE,
              aes(x = x_mid, y = y.position * 10^0.05,
                  label = p.adj.signif),
              size = 3.5, vjust = 0, colour = "black")
}

print(p_abs_jitter)

dir.create(file.path(output_dir, "FigS5C"), showWarnings = FALSE)
ggsave(file.path(output_dir, "FigS5C", "UMAP_Exp333_Tet_cluster_abs_jitter.pdf"),
       plot   = p_abs_jitter,
       width  = 3 * length(unique(abs_df$LN_type)),
       height = 2.8 * length(cluster_order))
save_source_data(p_abs_jitter, file.path(output_dir, "FigS5C", "UMAP_Exp333_Tet_cluster_abs_jitter.pdf"))
cat("Saved FigS5C/UMAP_Exp333_Tet_cluster_abs_jitter.pdf\n")

# ---- 9. Jitter plot — cluster proportions (%) --------------
# Pairwise Wilcoxon on proportions (both LN types)
pwc_prop <- prop_df %>%
  filter(!is.na(Group)) %>%
  mutate(Group = factor(trimws(Group), levels = group_order)) %>%
  group_by(cluster_name, LN_type) %>%
  pairwise_wilcox_test(prop ~ Group, p.adjust.method = "BH") %>%
  add_significance("p.adj")

pwc_prop_signif <- pwc_prop %>%
  add_xy_position(x = "Group", dodge = 0.8) %>%
  filter(p.adj <= 0.05) %>%
  mutate(y.position = y.position * 100)   # convert to % scale

prop_plot_df <- prop_df %>%
  filter(!is.na(Group)) %>%
  mutate(Group    = factor(trimws(Group), levels = group_order),
         prop_pct = prop * 100)

p_prop_jitter <- ggplot(prop_plot_df,
                         aes(x = Group, y = prop_pct,
                             colour = cluster_name)) +
  geom_boxplot(aes(group = Group), outlier.shape = NA,
               colour = "grey40", fill = NA, linewidth = 0.5) +
  geom_jitter(aes(shape = Group), width = 0.15, size = 2.5, alpha = 0.8) +
  scale_colour_manual(values = okabe_ito_5) +
  scale_shape_manual(values = group_shapes, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.18))) +
  facet_wrap(cluster_name ~ LN_type, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "% of Tet+ CD8+ T cells") +
  theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    axis.text.x      = element_text(angle = 35, hjust = 1),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey92", colour = NA),
    legend.position  = "none"
  )

if (nrow(pwc_prop_signif) > 0) {
  p_prop_jitter <- p_prop_jitter +
    stat_pvalue_manual(pwc_prop_signif, label = "p.adj.signif",
                       hide.ns = TRUE, tip.length = 0.01, size = 3.5)
}

print(p_prop_jitter)

dir.create(file.path(output_dir, "FigS5B"), showWarnings = FALSE)
ggsave(file.path(output_dir, "FigS5B", "UMAP_Exp333_Tet_cluster_prop_jitter.pdf"),
       plot   = p_prop_jitter,
       width  = 3 * length(unique(prop_df$LN_type)),
       height = 2.8 * length(cluster_order))
save_source_data(p_prop_jitter, file.path(output_dir, "FigS5B", "UMAP_Exp333_Tet_cluster_prop_jitter.pdf"))
cat("Saved FigS5B/UMAP_Exp333_Tet_cluster_prop_jitter.pdf\n")

# ---- 10. Export to Excel -----------------------------------
write_xlsx(
  list(
    "KW_absolute"         = as.data.frame(kw_abs),
    "Pairwise_absolute"   = as.data.frame(pwc_abs),
    "Per_sample_absolute" = as.data.frame(abs_df %>% select(
      file_name_short, Group, LN_type, cluster_name, n, prop, tet_e3, abs_cells
    ))
  ),
  path = "UMAP_Exp333_Tet_cluster_abs_stats.xlsx"
)
cat("Saved UMAP_Exp333_Tet_cluster_abs_stats.xlsx\n")

# ---- 11. Summary -------------------------------------------
cat("\n============================================================\n")
cat("SUMMARY — significant pairwise differences, ABSOLUTE COUNTS\n")
cat("============================================================\n")
sig_abs <- pwc_abs %>%
  filter(p.adj <= 0.05) %>%
  select(cluster_name, LN_type, group1, group2, n1, n2, p, p.adj, p.adj.signif) %>%
  arrange(cluster_name, LN_type)

if (nrow(sig_abs) == 0) {
  cat("No significant differences after BH correction.\n")
} else {
  print(sig_abs, n = Inf)
}

message("\nDone!")
