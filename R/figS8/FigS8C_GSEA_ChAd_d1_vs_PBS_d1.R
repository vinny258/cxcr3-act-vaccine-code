# ============================================================
# Figure:  Figure S8C
# Title:   GSEA — upregulated pathways in ChAd vs PBS at Day 1
#          (bulk RNA-seq of B16-OVA tumours, mouse)
# Method:  Gene ranking by DESeq2 log2FC (ChAd/PBS direction);
#          fgsea via clusterProfiler; top 10 per database shown
#          Gene sets: MSigDB Hallmark + Reactome (Mus musculus)
# Note:    log2FC column in source is PBS/ChAd — negated in script
#          to give ChAd/PBS direction; positive NES = up in ChAd
# Input:   FigS8C_GSEA_ChAd_d1_vs_PBS_d1.xlsx (sheet: biomarker_raw)
# Output:  FigS8C_GSEA_ChAd_d1_vs_PBS_d1.pdf  (bubble plot)
#          GSEA_ChAd_d1_vs_PBS_d1_Hallmark_Reactome.xlsx (results)
# ============================================================

# ---- 0. Install packages -----------------------------------

bioc_packages <- c("clusterProfiler", "enrichplot")
bioc_installed <- bioc_packages %in% rownames(installed.packages())

cran_packages <- c("msigdbr", "ggplot2", "readxl", "openxlsx", "dplyr", "stringr")
cran_installed <- cran_packages %in% rownames(installed.packages())

lapply(c(bioc_packages, cran_packages), library, character.only = TRUE)

# Location-independent paths and fonts. See helpers/paths.R
source(file.path({
  a <- commandArgs(FALSE); ff <- sub("^--file=", "", a[grep("^--file=", a)])
  dd <- if (length(ff)) dirname(ff) else getwd()
  while (!file.exists(file.path(dd, "helpers", "paths.R")) && dirname(dd) != dd) dd <- dirname(dd)
  dd
}, "helpers", "paths.R"))
setwd(out_dir())


# ---- 1. Load data ------------------------------------------
# working directory handled by helpers/paths.R
df <- read_excel(file.path(data_dir("biomarker.xlsx"), "biomarker.xlsx"), sheet = "biomarker_raw")
df <- as.data.frame(df)
cat("Loaded:", nrow(df), "genes\n")

# ---- 2. Build ranked gene list (ChAd Day 1 vs PBS Day 1) ---
fc_col <- "PBSday1_vs_ChAdday1_DESeq2_log2FC"

df[[fc_col]] <- suppressWarnings(as.numeric(df[[fc_col]]))

df <- df[
  !is.na(df[[fc_col]]) &
  !is.na(df$gene_name)  &
  df$gene_name != ""    &
  df$gene_name != "--",
]

# Positive log2FC in this column = upregulated in ChAd vs PBS
df$rank_metric <- df[[fc_col]]

df <- df[order(-abs(df$rank_metric)), ]
df <- df[!duplicated(df$gene_name), ]

gene_list <- setNames(df$rank_metric, df$gene_name)
gene_list <- sort(gene_list, decreasing = TRUE)
cat("Total genes in ranked list:", length(gene_list), "\n")

# ---- 3. Download gene sets for mouse -----------------------
cat("\nFetching MSigDB gene sets for Mus musculus...\n")

hallmark <- msigdbr(species = "Mus musculus", category = "H")
reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")

term2gene_hallmark <- hallmark %>% dplyr::select(gs_name, gene_symbol)
term2gene_reactome <- reactome %>% dplyr::select(gs_name, gene_symbol)

cat("Hallmark gene sets:", n_distinct(hallmark$gs_name), "\n")
cat("Reactome gene sets:", n_distinct(reactome$gs_name), "\n")

# ---- Helper ------------------------------------------------
has_results <- function(x) !is.null(x) && nrow(as.data.frame(x)) > 0

# ---- 4. GSEA — Hallmark ------------------------------------
cat("\nRunning GSEA - Hallmark...\n")
gsea_hallmark <- tryCatch(
  GSEA(
    geneList     = gene_list,
    TERM2GENE    = term2gene_hallmark,
    minGSSize    = 15,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    verbose      = FALSE,
    eps          = 0
  ),
  error = function(e) { cat("Hallmark GSEA error:", conditionMessage(e), "\n"); NULL }
)
cat("Significant Hallmark gene sets:",
    if (has_results(gsea_hallmark)) nrow(as.data.frame(gsea_hallmark)) else 0, "\n")

# ---- 5. GSEA — Reactome ------------------------------------
cat("\nRunning GSEA - Reactome...\n")
gsea_reactome <- tryCatch(
  GSEA(
    geneList     = gene_list,
    TERM2GENE    = term2gene_reactome,
    minGSSize    = 15,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    verbose      = FALSE,
    eps          = 0
  ),
  error = function(e) { cat("Reactome GSEA error:", conditionMessage(e), "\n"); NULL }
)
cat("Significant Reactome gene sets:",
    if (has_results(gsea_reactome)) nrow(as.data.frame(gsea_reactome)) else 0, "\n")

# ---- 6. Build combined data frame for plotting -------------
extract_top <- function(gsea_obj, db_name, prefix, top_n = 15) {
  if (!has_results(gsea_obj)) return(NULL)
  df <- as.data.frame(gsea_obj)
  df <- df[df$NES > 0, ]
  if (nrow(df) == 0) return(NULL)
  df <- df[order(-df$NES), ]
  df <- head(df, top_n)
  df$database <- db_name
  df$label    <- str_wrap(gsub("_", " ", tolower(sub(prefix, "", df$ID))), width = 45)
  df
}

hall_df  <- extract_top(gsea_hallmark, "Hallmark", "HALLMARK_", top_n = 10)
react_df <- extract_top(gsea_reactome, "Reactome", "REACTOME_", top_n = 10)

plot_df <- bind_rows(hall_df, react_df)

# ---- 7. Combined dotplot -----------------------------------
p_combined <- NULL

if (!is.null(plot_df) && nrow(plot_df) > 0) {

  plot_df <- plot_df %>%
    arrange(database, NES) %>%
    mutate(
      label    = factor(label, levels = unique(label)),
      database = factor(database, levels = c("Hallmark", "Reactome"))
    )

  p_combined <- ggplot(plot_df,
                       aes(x = NES, y = label, size = setSize, colour = p.adjust)) +
    geom_point() +
    scale_colour_gradient(
      low   = "#B22222",
      high  = "#6495ED",
      name  = "p.adjust",
      guide = guide_colorbar(reverse = TRUE)
    ) +
    scale_size_continuous(name = "Gene set size", range = c(3, 9)) +
    facet_grid(database ~ ., scales = "free_y", space = "free_y") +
    coord_cartesian(xlim = c(1.8, NA)) +
    labs(
      title    = "GSEA — Hallmark & Reactome (ChAd Day 1 vs PBS Day 1)",
      subtitle = "Upregulated in ChAd Day 1 (NES > 0), top 10 per database",
      x        = "Normalised Enrichment Score (NES)",
      y        = NULL
    ) +
    theme_classic(base_size = 16) +
    theme(
      plot.title       = element_text(face = "bold", size = 16),
      plot.subtitle    = element_text(size = 16),
      axis.title.x     = element_text(size = 16),
      axis.text.x      = element_text(size = 16),
      axis.text.y      = element_text(size = 11, lineheight = 0.85),
      legend.title     = element_text(size = 16),
      legend.text      = element_text(size = 16),
      strip.text       = element_text(size = 16, face = "bold"),
      strip.background = element_rect(fill = "grey90", colour = NA),
      panel.spacing    = unit(1, "lines")
    )

  tryCatch(print(p_combined), error = function(e) {
    grid::grid.newpage()
    grid::grid.draw(ggplotGrob(p_combined))
  })

} else {
  cat("No upregulated pathways found in either database\n")
}

# ---- 8. Save results to Excel ------------------------------
wb <- createWorkbook()

if (has_results(gsea_hallmark)) {
  hall_out <- as.data.frame(gsea_hallmark)[, c("ID", "setSize", "enrichmentScore",
                                                "NES", "pvalue", "p.adjust",
                                                "core_enrichment")]
  addWorksheet(wb, "Hallmark")
  writeData(wb, "Hallmark", hall_out)
}

if (has_results(gsea_reactome)) {
  react_out <- as.data.frame(gsea_reactome)[, c("ID", "setSize", "enrichmentScore",
                                                  "NES", "pvalue", "p.adjust",
                                                  "core_enrichment")]
  addWorksheet(wb, "Reactome")
  writeData(wb, "Reactome", react_out)
}

saveWorkbook(wb, "GSEA_ChAd_d1_vs_PBS_d1_Hallmark_Reactome.xlsx", overwrite = TRUE)
cat("Results saved to GSEA_ChAd_d1_vs_PBS_d1_Hallmark_Reactome.xlsx\n")

# ---- 9. Save PDF ------------------------------------------
tryCatch({
  if (!is.null(p_combined))
    ggsave("GSEA_ChAd_d1_vs_PBS_d1_Hallmark_Reactome_dotplot.pdf",
           plot = p_combined, width = 14, height = 12, device = "pdf")
    save_source_data(p_combined, "GSEA_ChAd_d1_vs_PBS_d1_Hallmark_Reactome_dotplot.pdf")
  cat("Plot saved to GSEA_ChAd_d1_vs_PBS_d1_Hallmark_Reactome_dotplot.pdf\n")
}, error = function(e) cat("PDF saving skipped:", conditionMessage(e), "\n"))

message("Done!")
