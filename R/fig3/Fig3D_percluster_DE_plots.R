# ============================================================
# Figure:  Figure 3D (alternative visualisations)
# Title:   Per-cluster pseudobulk DESeq2 — Cxcl10 & Cxcl9 (D1 vs D0)
#          WT draining lymph node cells, scRNA-seq
# Method:  DESeq2 run separately per cell cluster; paired design
#          (~ mouse + day); sub-clusters merged into broad cell types.
#          Expression quantified from raw RNA counts (RNA assay);
#          log-normalisation (log1p, scale factor 10,000) computed
#          on-the-fly (NormalizeData was not run; object uses SCTransform).
#          Expressing cells defined as raw count >= 1 UMI.
#
# Input:   ../shared data/wt_D0_D1.rds  (Seurat object, SCTransform)
#          Germinal centre B cells (clusters 17_GCB, 23_GCB) excluded
#          due to insufficient cell numbers (n = 19 D0, n = 10 D1).
#
# Output:  Fig3D_optionA_dotplot_significance.pdf
#            — dot plot: dot size = % expressing, colour = log2FC (D1/D0),
#              significance markers from per-cluster DESeq2 (* p<0.05,
#              ** p<0.01, *** p<0.001, Benjamini-Hochberg padj)
#          Fig3D_optionB_lollipop_FC_CI.pdf
#            — lollipop: log2FC ± 95% CI per cluster, coloured by
#              significance (padj < 0.05)
#          Fig3D_optionC_cell_counts.pdf
#            — horizontal bar chart: absolute number of expressing cells
#              per cluster x timepoint (Day 0 grey, Day 1 red)
#          Fig3D_optionD_dotplot_expr.pdf
#            — dot plot: dot size = n expressing cells,
#              colour = mean log-normalised expression among expressors;
#              DESeq2 significance markers overlaid on Day 1 dots
#          Fig3D_optionE_scatter_dominance.pdf
#            — scatter plot: x = n expressing cells (breadth),
#              y = mean expression (intensity), bubble size = total
#              transcriptional output (n x mean expr); coloured by
#              broad lineage (Myeloid / Innate lymphoid / T cells /
#              B cells); Day 1 only; faceted by gene
#          FigS_optionF_scatter_dominance_D0.pdf
#            — same as Option E but for Day 0 (baseline); supplementary
#
# Intermediate files saved to disk (for fast re-plotting):
#          pct_long.rds   — per-cluster x day x gene expression summary
#          dot_df.rds     — merged expression + count + DESeq2 sig data
#          de_merged.rds  — per-cluster DESeq2 results (merged cell types)
# ============================================================

# ---- 0. Packages -------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2",      quietly = TRUE)) BiocManager::install("DESeq2")
for (p in c("ggplot2", "dplyr", "tidyr"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

library(Seurat)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("C:/Users/valmeida/Ludwig Institute for Cancer Research Dropbox/Vinnycius Pereira Almeida/BVDE Lab/Lab members/Vinny/DPhil Clinical Medicine/ATC + vaccines paper - VPA/Science Immunology/R scripts/Fig3D_percluster_DE_plots")

GENES <- c("Cxcl10", "Cxcl9")

# ---- 1. Load & subset Seurat (WT singlets only) ------------
cat("Loading Seurat object...\n")
seu <- readRDS("../shared data/wt_D0_D1.rds")

seu_wt <- subset(seu,
                 group   == "WT" &
                 hash.ID %in% grep("^HTO", unique(seu$hash.ID), value = TRUE) &
                 !grepl("GCB", clusterid))
DefaultAssay(seu_wt) <- "RNA"
cat("WT cells (GCB excluded):", ncol(seu_wt), "\n")

# ---- 2. Per-cluster pseudobulk DESeq2 ----------------------
# Run on ALL genes per cluster so DESeq2 can estimate dispersion
# then extract Cxcl10/Cxcl9 results
cat("\nRunning per-cluster pseudobulk DESeq2...\n")

clusters <- sort(unique(seu_wt$clusterid))
de_list  <- list()

for (cl in clusters) {
  cat(" Cluster:", cl, "... ")

  cells_cl  <- colnames(seu_wt)[seu_wt$clusterid == cl]
  meta_cl   <- seu_wt@meta.data[cells_cl, c("hash.ID", "day")]
  days_present <- names(table(meta_cl$day)[table(meta_cl$day) > 0])

  # Skip if both days not represented
  if (!all(c("D0", "D1") %in% days_present)) {
    cat("skipped (missing day)\n"); next
  }

  seu_cl <- subset(seu_wt, cells = cells_cl)

  pb <- tryCatch(
    AggregateExpression(seu_cl, assays = "RNA",
                        return.seurat = FALSE,
                        group.by = c("hash.ID", "day"))$RNA,
    error = function(e) NULL
  )
  if (is.null(pb)) { cat("skipped (aggregation failed)\n"); next }

  sample_names <- colnames(pb)
  coldata <- data.frame(
    row.names = sample_names,
    mouse     = sub("_D[01]$", "", sample_names),
    day       = sub(".*_(D[01])$", "\\1", sample_names),
    stringsAsFactors = TRUE
  )
  coldata$day   <- relevel(factor(coldata$day),   ref = "D0")
  coldata$mouse <- factor(coldata$mouse)

  # Need ≥ 2 samples per day
  if (any(table(coldata$day) < 2)) {
    cat("skipped (< 2 samples per day)\n"); next
  }

  # Filter low-count genes before DESeq2
  keep <- rowSums(pb >= 10) >= 2
  pb   <- pb[keep, , drop = FALSE]

  if (nrow(pb) < 10) {
    cat("skipped (too few genes after filtering)\n"); next
  }

  # Check target genes are present
  genes_present <- intersect(GENES, rownames(pb))
  if (length(genes_present) == 0) {
    cat("skipped (target genes absent)\n"); next
  }

  tryCatch({
    dds <- DESeqDataSetFromMatrix(pb, colData = coldata, design = ~ mouse + day)
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("day", "D1", "D0"), alpha = 0.05)

    for (g in genes_present) {
      r <- as.data.frame(res[g, ])
      de_list[[paste0(cl, "__", g)]] <- data.frame(
        clusterid     = cl,
        gene          = g,
        log2FC        = r$log2FoldChange,
        lfcSE         = r$lfcSE,
        padj          = r$padj,
        stringsAsFactors = FALSE
      )
    }
    cat("OK\n")
  }, error = function(e) cat("error:", conditionMessage(e), "\n"))
}

de_df <- bind_rows(de_list)
cat("\nClusters with results:", n_distinct(de_df$clusterid), "\n")
cat("Genes recovered:", paste(unique(de_df$gene), collapse = ", "), "\n")

# Significance labels
de_df <- de_df %>%
  mutate(
    sig = case_when(
      is.na(padj)    ~ "",
      padj < 0.001   ~ "***",
      padj < 0.01    ~ "**",
      padj < 0.05    ~ "*",
      TRUE           ~ ""
    ),
    ci_lo = log2FC - 1.96 * lfcSE,
    ci_hi = log2FC + 1.96 * lfcSE
  )

cat("\nSignificant clusters (padj < 0.05):\n")
print(de_df %>% filter(sig != "") %>%
        select(clusterid, gene, log2FC, padj, sig) %>%
        arrange(gene, padj))

# ---- 3. Merge sub-clusters into broad cell types -----------
# Strip numeric prefix: "0_CD8Mem" -> "CD8Mem"
# log2FC  : mean across sub-clusters
# padj    : minimum (most significant sub-cluster represents the type)
# lfcSE   : from the sub-cluster with the minimum padj
de_df <- de_df %>%
  mutate(cell_type = gsub("^[0-9]+_", "", clusterid))

de_merged <- de_df %>%
  group_by(cell_type, gene) %>%
  summarise(
    log2FC = mean(log2FC, na.rm = TRUE),
    padj   = ifelse(all(is.na(padj)), NA, min(padj, na.rm = TRUE)),
    lfcSE  = lfcSE[which.min(replace(padj, is.na(padj), Inf))],
    .groups = "drop"
  ) %>%
  mutate(
    sig = case_when(
      is.na(padj)  ~ "",
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ ""
    ),
    ci_lo = log2FC - 1.96 * lfcSE,
    ci_hi = log2FC + 1.96 * lfcSE
  )

cat("\nSignificant cell types after merging (padj < 0.05):\n")
print(de_merged %>% filter(sig != "") %>%
        select(cell_type, gene, log2FC, padj, sig) %>%
        arrange(gene, padj))

saveRDS(de_merged, "de_merged.rds")

# ---- 4. Compute pct expressing per cell type × day ---------
meta       <- seu_wt@meta.data[, c("clusterid", "day")]
counts_mat <- GetAssayData(seu_wt, assay = "RNA", layer = "counts")

pct_list <- lapply(GENES, function(g) {
  expressed <- as.numeric(counts_mat[g, ]) > 0
  meta %>%
    mutate(expressed  = expressed,
           cell_type  = gsub("^[0-9]+_", "", clusterid)) %>%
    group_by(cell_type, day) %>%
    summarise(
      n_expressing = sum(expressed),
      n_total      = n(),
      pct          = mean(expressed),
      .groups = "drop"
    ) %>%
    mutate(gene = g)
})

pct_long <- bind_rows(pct_list)
saveRDS(pct_long, "pct_long.rds")

pct_df <- pct_long %>%
  pivot_wider(
    id_cols     = c(cell_type, gene),
    names_from  = day,
    values_from = c(pct, n_expressing, n_total),
    names_glue  = "{.value}_{day}"
  )

# Merge pct with merged DE results
plot_df <- pct_df %>%
  left_join(de_merged, by = c("cell_type", "gene")) %>%
  mutate(gene = factor(gene, levels = GENES))

# ---- 5. OPTION A: Dotplot with significance markers --------
cat("\nGenerating Option A — dotplot with significance markers...\n")

# Remove cell types where DESeq2 could not be run (NA log2FC)
valid_types <- de_merged %>%
  group_by(cell_type) %>%
  filter(!any(is.na(log2FC))) %>%
  pull(cell_type) %>%
  unique()

plot_df_A <- plot_df %>% filter(cell_type %in% valid_types)
cat("Cell types excluded (insufficient data for DESeq2):",
    paste(setdiff(unique(plot_df$cell_type), valid_types), collapse = ", "), "\n")

# Rename cell types for display
label_map <- c(
  "Mono"        = "Monocytes",
  "Tgd"         = "\u03b3\u03b4 T cells",
  "Neu"         = "Neutrophils",
  "DCs"         = "Dendritic cells",
  "Tregs"       = "Regulatory T cells",
  "CD8Effector" = "Effector CD8\u207a T cells",
  "Ifn_CD8"     = "Type I IFN-stimulated CD8\u207a T cells",
  "NK"          = "Natural Killer cells",
  "CD8Mem"      = "Memory CD8\u207a T cells",
  "CD4Naive"    = "Na\u00efve CD4\u207a T cells",
  "NaiveB"      = "Na\u00efve B cells",
  "GCB"         = "Germinal centre B cells"
)

plot_df_A <- plot_df_A %>%
  mutate(cell_type = unname(ifelse(cell_type %in% names(label_map),
                                   label_map[cell_type], cell_type)))

sig_label_map <- label_map

# Order rows by Cxcl10 pct_D1 descending
row_order <- plot_df_A %>%
  filter(gene == "Cxcl10") %>%
  arrange(desc(pct_D1)) %>%
  pull(cell_type)
all_labels <- unique(plot_df_A$cell_type)
row_order  <- c(row_order, setdiff(all_labels, row_order))

plot_df_A$cell_type <- factor(plot_df_A$cell_type, levels = rev(row_order))

# Long format for dots
df_long <- plot_df_A %>%
  pivot_longer(cols = c(pct_D0, pct_D1),
               names_to = "timepoint", values_to = "pct") %>%
  mutate(timepoint = ifelse(timepoint == "pct_D1", "Day 1", "Day 0"),
         timepoint = factor(timepoint, levels = c("Day 0", "Day 1")))

# Significance label positions (Day 1 only)
sig_df <- plot_df_A %>%
  filter(!is.na(sig), sig != "") %>%
  select(cell_type, gene, pct_D1, log2FC, sig) %>%
  mutate(timepoint = factor("Day 1", levels = c("Day 0", "Day 1")),
         cell_type = factor(cell_type, levels = levels(plot_df_A$cell_type)))

pA <- ggplot(df_long, aes(x = timepoint, y = cell_type, size = pct)) +

  geom_point(data = subset(df_long, timepoint == "Day 0"),
             colour = "grey60", alpha = 0.8) +

  geom_point(data = subset(df_long, timepoint == "Day 1"),
             aes(colour = log2FC), alpha = 0.9) +

  # Significance markers above Day 1 dots
  geom_text(data = sig_df,
            aes(label = sig, x = timepoint, y = cell_type),
            size = 4, vjust = -0.8, colour = "black", inherit.aes = FALSE) +

  scale_colour_gradient2(
    low      = "#3498DB",
    mid      = "#F5CBA7",
    high     = "#C0392B",
    midpoint = 0,
    name     = "log2FC\n(D1 vs D0)",
    na.value = "grey80"
  ) +

  scale_size_continuous(
    name   = "% expressing",
    range  = c(1, 10),
    breaks = c(0.1, 0.3, 0.5, 0.7, 0.9),
    labels = c("10%", "30%", "50%", "70%", "90%")
  ) +

  facet_grid(. ~ gene) +

  labs(
    title    = "Cxcl10 and Cxcl9 expression across cell clusters",
    subtitle = "WT cells | pseudobulk DESeq2 | * p<0.05  ** p<0.01  *** p<0.001",
    x = NULL, y = NULL
  ) +

  theme_classic(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    axis.text.x      = element_text(size = 12),
    axis.text.y      = element_text(size = 10),
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    strip.text       = element_text(size = 12, face = "bold.italic"),
    strip.background = element_rect(fill = "grey92", colour = NA),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
    panel.spacing    = unit(1.5, "lines")
  )

print(pA)
ggsave("Fig3D_optionA_dotplot_significance.pdf",
       plot = pA, width = 7, height = 10, device = "pdf")
cat("Option A saved to Fig3D_optionA_dotplot_significance.pdf\n")

# ---- 6. OPTION B: Lollipop plot log2FC ± 95% CI ------------
cat("\nGenerating Option B — lollipop plot log2FC ± 95% CI...\n")

lollipop_df <- de_merged %>%
  filter(!is.na(log2FC)) %>%
  mutate(
    gene      = factor(gene, levels = GENES),
    sig_col   = ifelse(sig != "", "Significant (padj < 0.05)", "Not significant"),
    cell_type = reorder(cell_type, log2FC)
  )

pB <- ggplot(lollipop_df,
             aes(x = log2FC, y = cell_type, colour = sig_col)) +

  # CI segment
  geom_segment(aes(x = ci_lo, xend = ci_hi,
                   y = cell_type, yend = cell_type),
               linewidth = 0.8, alpha = 0.7) +

  # Point
  geom_point(size = 3.5) +

  # Zero line
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +

  scale_colour_manual(
    values = c("Significant (padj < 0.05)" = "#C0392B",
               "Not significant"            = "grey60"),
    name = NULL
  ) +

  facet_grid(. ~ gene, scales = "free_x") +

  labs(
    title    = "Cxcl10 and Cxcl9: log2FC per cluster (D1 vs D0)",
    subtitle = "Pseudobulk DESeq2 | Error bars = 95% CI | Ordered by log2FC",
    x        = "log2 Fold Change (D1 / D0)",
    y        = NULL
  ) +

  theme_classic(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    axis.text.y      = element_text(size = 10),
    axis.text.x      = element_text(size = 11),
    axis.title.x     = element_text(size = 12),
    legend.position  = "bottom",
    legend.text      = element_text(size = 11),
    strip.text       = element_text(size = 12, face = "bold.italic"),
    strip.background = element_rect(fill = "grey92", colour = NA),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.4),
    panel.spacing    = unit(2, "lines")
  )

print(pB)
ggsave("Fig3D_optionB_lollipop_FC_CI.pdf",
       plot = pB, width = 10, height = 8, device = "pdf")
cat("Option B saved to Fig3D_optionB_lollipop_FC_CI.pdf\n")

# ---- 7. OPTION C: Expressing cell count bar chart ----------
cat("\nGenerating Option C — expressing cell counts per cluster...\n")

count_df <- pct_long %>%
  mutate(
    cell_type_label = unname(ifelse(cell_type %in% names(label_map),
                                    label_map[cell_type], cell_type)),
    timepoint       = factor(day, levels = c("D0", "D1"),
                             labels = c("Day 0", "Day 1")),
    gene            = factor(gene, levels = GENES)
  )

# Order by total Cxcl10 n_expressing (D0+D1), highest at top
count_order <- count_df %>%
  filter(gene == "Cxcl10") %>%
  group_by(cell_type_label) %>%
  summarise(total = sum(n_expressing), .groups = "drop") %>%
  arrange(total) %>%
  pull(cell_type_label)
count_df$cell_type_label <- factor(count_df$cell_type_label, levels = count_order)

pC <- ggplot(count_df,
             aes(x = n_expressing, y = cell_type_label, fill = timepoint)) +

  geom_col(position = position_dodge(width = 0.7), width = 0.6) +

  geom_text(aes(label = n_expressing),
            position = position_dodge(width = 0.7),
            hjust = -0.15, size = 3.2, colour = "grey20") +

  scale_fill_manual(
    values = c("Day 0" = "grey70", "Day 1" = "#C0392B"),
    name   = NULL
  ) +

  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +

  facet_grid(. ~ gene, scales = "free_x") +

  labs(
    title    = "Cells expressing Cxcl10 and Cxcl9 per cluster",
    subtitle = "WT cells | count > 0 | ordered by Cxcl10 Day 1",
    x        = "Number of expressing cells",
    y        = NULL
  ) +

  theme_classic(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 9, colour = "grey40"),
    axis.text.y        = element_text(size = 10),
    axis.text.x        = element_text(size = 11),
    axis.title.x       = element_text(size = 12),
    legend.position    = "bottom",
    legend.text        = element_text(size = 11),
    strip.text         = element_text(size = 12, face = "bold.italic"),
    strip.background   = element_rect(fill = "grey92", colour = NA),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.4),
    panel.spacing      = unit(2, "lines")
  )

print(pC)
ggsave("Fig3D_optionC_cell_counts.pdf",
       plot = pC, width = 10, height = 8, device = cairo_pdf)
cat("Option C saved to Fig3D_optionC_cell_counts.pdf\n")

# ---- 8. OPTION D: Dot plot — n expressing + mean expression ----
cat("\nGenerating Option D — dot plot: n expressing + mean expression...\n")

# RNA data layer contains raw counts (NormalizeData was never run — SCTransform only)
# Normalize on the fly: log1p(counts / library_size * 10000)
counts_raw <- GetAssayData(seu_wt, assay = "RNA", layer = "counts")
lib_sizes  <- colSums(counts_raw)
norm_mat   <- log1p(sweep(counts_raw, 2, lib_sizes, "/") * 10000)

expr_list <- lapply(GENES, function(g) {
  expr_vals <- as.numeric(norm_mat[g, ])
  meta %>%
    mutate(
      expr      = expr_vals,
      cell_type = gsub("^[0-9]+_", "", clusterid)
    ) %>%
    group_by(cell_type, day) %>%
    summarise(
      mean_expr = { v <- expr[expr > 0]; if (length(v) > 0) mean(v) else NA_real_ },
      .groups = "drop"
    ) %>%
    mutate(gene = g)
})

expr_df <- bind_rows(expr_list)

dot_df <- pct_long %>%
  left_join(expr_df, by = c("cell_type", "day", "gene")) %>%
  left_join(de_merged %>% select(cell_type, gene, sig, padj),
            by = c("cell_type", "gene")) %>%
  mutate(
    cell_type_label = unname(ifelse(cell_type %in% names(label_map),
                                    label_map[cell_type], cell_type)),
    timepoint       = factor(day, levels = c("D0", "D1"),
                             labels = c("Day 0", "Day 1")),
    gene            = factor(gene, levels = GENES)
  )

# Order by mean log-normalised Cxcl10 expression at Day 1 (highest at top)
dot_order <- dot_df %>%
  filter(gene == "Cxcl10", day == "D1") %>%
  arrange(mean_expr) %>%
  pull(cell_type_label)
dot_df$cell_type_label <- factor(dot_df$cell_type_label, levels = dot_order)

saveRDS(dot_df, "dot_df.rds")

sig_d1 <- dot_df %>%
  filter(day == "D1", !is.na(sig), sig != "")

pD <- ggplot(dot_df, aes(x = timepoint, y = cell_type_label)) +

  geom_point(aes(size = n_expressing, colour = mean_expr), alpha = 0.9) +

  geom_text(data = sig_d1,
            aes(label = sig, x = timepoint, y = cell_type_label),
            size = 3.5, nudge_y = 0.35, colour = "black",
            inherit.aes = FALSE) +

  scale_colour_gradient(
    low      = "#F5CBA7",
    high     = "#C0392B",
    name     = "Mean expression\n(log-normalised)",
    na.value = "grey85"
  ) +

  scale_size_continuous(
    name   = "Expressing cells (n)",
    range  = c(1, 10)
  ) +

  facet_grid(. ~ gene) +

  labs(
    title    = "Cxcl10 and Cxcl9 expression per cluster",
    subtitle = "Dot size = no. expressing cells  |  Colour = mean log-normalised expression in expressing cells",
    x = NULL, y = NULL
  ) +

  theme_classic(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    axis.text.x      = element_text(size = 12),
    axis.text.y      = element_text(size = 10),
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 9),
    strip.text       = element_text(size = 12, face = "bold.italic"),
    strip.background = element_rect(fill = "grey92", colour = NA),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
    panel.spacing    = unit(1.5, "lines")
  )

print(pD)
ggsave("Fig3D_optionD_dotplot_expr.pdf",
       plot = pD, width = 8, height = 10, device = cairo_pdf)
cat("Option D saved to Fig3D_optionD_dotplot_expr.pdf\n")

# ---- 9. OPTION E: Scatter plot — transcriptional dominance ----
cat("\nGenerating Option E — scatter plot: transcriptional dominance...\n")

if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

scatter_df <- dot_df %>%
  filter(day == "D1", !is.na(mean_expr)) %>%
  mutate(
    total_expr = n_expressing * mean_expr,
    gene       = factor(gene, levels = GENES),
    lineage    = case_when(
      cell_type_label %in% c("Monocytes", "Neutrophils", "Dendritic cells") ~ "Myeloid",
      grepl("Natural Killer|γδ", cell_type_label)                  ~ "Innate lymphoid",
      grepl("B cells", cell_type_label)                                       ~ "B cells",
      TRUE                                                                     ~ "T cells"
    )
  )

lineage_colors <- c(
  "Myeloid"         = "#E74C3C",
  "Innate lymphoid" = "#F39C12",
  "T cells"         = "#3498DB",
  "B cells"         = "#95A5A6"
)

pE <- ggplot(scatter_df,
             aes(x = n_expressing, y = mean_expr,
                 size = total_expr, colour = lineage,
                 label = cell_type_label)) +

  geom_point(alpha = 0.8) +

  geom_text_repel(size = 3.2, show.legend = FALSE,
                  max.overlaps = 20, seed = 42) +

  scale_colour_manual(values = lineage_colors, name = "Lineage") +

  scale_size_continuous(
    name  = "Total output\n(n × mean expr)",
    range = c(2, 12)
  ) +

  facet_wrap(~ gene, scales = "free_x") +

  labs(
    title    = "Cxcl10 and Cxcl9 — transcriptional dominance per cluster",
    subtitle = "Day 1 | X = no. expressing cells (breadth)  |  Y = mean expression (intensity)  |  Size = total output",
    x        = "Number of expressing cells",
    y        = "Mean log-normalised expression\n(expressing cells only)"
  ) +

  theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    strip.text       = element_text(size = 12, face = "bold.italic"),
    strip.background = element_rect(fill = "grey92", colour = NA),
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 9),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.4)
  )

print(pE)
ggsave("Fig3D_optionE_scatter_dominance.pdf",
       plot = pE, width = 12, height = 6, device = cairo_pdf)
cat("Option E saved to Fig3D_optionE_scatter_dominance.pdf\n")

# ---- 10. OPTION F: Scatter — transcriptional dominance, Day 0 (baseline) ----
cat("\nGenerating Option F — scatter plot: transcriptional dominance at Day 0...\n")

scatter_df_D0 <- dot_df %>%
  filter(day == "D0", !is.na(mean_expr)) %>%
  mutate(
    total_expr = n_expressing * mean_expr,
    gene       = factor(gene, levels = GENES),
    lineage    = case_when(
      cell_type_label %in% c("Monocytes", "Neutrophils", "Dendritic cells") ~ "Myeloid",
      grepl("Natural Killer|γδ", cell_type_label)                  ~ "Innate lymphoid",
      grepl("B cells", cell_type_label)                                       ~ "B cells",
      TRUE                                                                     ~ "T cells"
    )
  )

pF <- ggplot(scatter_df_D0,
             aes(x = n_expressing, y = mean_expr,
                 size = total_expr, colour = lineage,
                 label = cell_type_label)) +

  geom_point(alpha = 0.8) +

  geom_text_repel(size = 3.2, show.legend = FALSE,
                  max.overlaps = 20, seed = 42) +

  scale_colour_manual(values = lineage_colors, name = "Lineage") +

  scale_size_continuous(
    name  = "Total output\n(n × mean expr)",
    range = c(2, 12)
  ) +

  facet_wrap(~ gene, scales = "free_x") +

  labs(
    title    = "Cxcl10 and Cxcl9 — transcriptional dominance per cluster (baseline)",
    subtitle = "Day 0 | X = no. expressing cells (breadth)  |  Y = mean expression (intensity)  |  Size = total output",
    x        = "Number of expressing cells",
    y        = "Mean log-normalised expression\n(expressing cells only)"
  ) +

  theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    strip.text       = element_text(size = 12, face = "bold.italic"),
    strip.background = element_rect(fill = "grey92", colour = NA),
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 9),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.4)
  )

print(pF)
ggsave("FigS_optionF_scatter_dominance_D0.pdf",
       plot = pF, width = 12, height = 6, device = cairo_pdf)
cat("Option F saved to FigS_optionF_scatter_dominance_D0.pdf\n")

message("Done!")
