# ============================================================
# Figure:  Figure 4J
# Title:   Cxcl9 expression boxplot — ChAd vs PBS, Days 1, 2 & 3
#          (bulk RNA-seq of B16-OVA tumours, mouse)
# Method:  log2(FPKM + 1); DESeq2 adjusted p-value annotated
#          Individual replicates shown as points
# Input:   Fig4J_cxcl9_boxplot.xlsx  (sheet: biomarker_raw)
# Output:  Fig4J_cxcl9_boxplot.pdf   (boxplot)
# ============================================================

# ---- 0. Packages -------------------------------------------
packages <- c("ggplot2", "readxl", "dplyr", "tidyr")
installed <- packages %in% rownames(installed.packages())

lapply(packages, library, character.only = TRUE)

# Location-independent paths and fonts. See helpers/paths.R
source(file.path({
  a <- commandArgs(FALSE); ff <- sub("^--file=", "", a[grep("^--file=", a)])
  dd <- if (length(ff)) dirname(ff) else getwd()
  while (!file.exists(file.path(dd, "helpers", "paths.R")) && dirname(dd) != dd) dd <- dirname(dd)
  dd
}, "helpers", "paths.R"))
setwd(out_dir())


# ---- 1. Load data ------------------------------------------
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
# working directory handled by helpers/paths.R
}

df <- read_excel(file.path(data_dir("biomarker.xlsx"), "biomarker.xlsx"), sheet = "biomarker_raw")
df <- as.data.frame(df)

# ---- 2. Extract CXCL9 row ----------------------------------
cxcl9 <- df[!is.na(df$gene_name) & tolower(df$gene_name) == "cxcl9", ]

if (nrow(cxcl9) == 0) stop("Cxcl9 not found in gene_name column.")

if (nrow(cxcl9) > 1) {
  cat("Multiple Cxcl9 rows found — using the one with highest mean FPKM.\n")
  fpkm_cols_all <- grep("_FPKM$", colnames(cxcl9), value = TRUE)
  means <- rowMeans(sapply(cxcl9[, fpkm_cols_all], as.numeric), na.rm = TRUE)
  cxcl9 <- cxcl9[which.max(means), ]
}

cat("Found Cxcl9 row — gene ID:", cxcl9$ID[1], "\n")

# ---- 3. Pivot FPKM columns to long format ------------------
fpkm_map <- list(
  list(day = "d-1", treatment = "ChAd", cols = c("ChAdd1a_FPKM", "ChAdd1b_FPKM", "ChAdd1c_FPKM")),
  list(day = "d-1", treatment = "PBS",  cols = c("PBSd1a_FPKM",  "PBSd1b_FPKM",  "PBSd1c_FPKM")),
  list(day = "d-2", treatment = "ChAd", cols = c("ChAdd2a_FPKM", "ChAdd2b_FPKM", "ChAdd2c_FPKM")),
  list(day = "d-2", treatment = "PBS",  cols = c("PBSd2a_FPKM",  "PBSd2b_FPKM",  "PBSd2c_FPKM")),
  list(day = "d-3", treatment = "ChAd", cols = c("ChAdd3a_FPKM", "ChAdd3b_FPKM", "ChAdd3c_FPKM")),
  list(day = "d-3", treatment = "PBS",  cols = c("PBSd3a_FPKM",  "PBSd3b_FPKM",  "PBSd3c_FPKM"))
)

long_list <- lapply(fpkm_map, function(entry) {
  vals <- suppressWarnings(as.numeric(cxcl9[, entry$cols]))
  data.frame(
    day       = entry$day,
    treatment = entry$treatment,
    fpkm      = vals,
    replicate = c("a", "b", "c"),
    stringsAsFactors = FALSE
  )
})

plot_df <- do.call(rbind, long_list)
plot_df$log2fpkm  <- log2(plot_df$fpkm + 1)
plot_df$day       <- factor(plot_df$day,       levels = c("d-1", "d-2", "d-3"))
plot_df$treatment <- factor(plot_df$treatment, levels = c("ChAd", "PBS"))

cat("log2(FPKM+1) values:\n")
print(plot_df[, c("day", "treatment", "replicate", "fpkm", "log2fpkm")])

# ---- 4. Extract DESeq2 p-values ----------------------------
pval_cols <- c(
  "d-1" = "PBSday1_vs_ChAdday1_DESeq2_Pvalue",
  "d-2" = "PBSday2_vs_ChAdday2_DESeq2_Pvalue",
  "d-3" = "PBSday3_vs_ChAdday3_DESeq2_Pvalue"
)

pval_df <- data.frame(
  day   = names(pval_cols),
  padj  = sapply(pval_cols, function(col) {
    suppressWarnings(as.numeric(cxcl9[[col]]))
  }),
  stringsAsFactors = FALSE
)
pval_df$day <- factor(pval_df$day, levels = c("d-1", "d-2", "d-3"))

# Format p-value label
pval_df$label <- ifelse(
  is.na(pval_df$padj), "p = NA",
  ifelse(pval_df$padj < 0.001, "p < 0.001",
  ifelse(pval_df$padj < 0.01,  paste0("p = ", formatC(pval_df$padj, format = "f", digits = 3)),
  ifelse(pval_df$padj < 0.05,  paste0("p = ", formatC(pval_df$padj, format = "f", digits = 3)),
                                paste0("p = ", formatC(pval_df$padj, format = "f", digits = 2)))))
)

# Position the label just above the max value in each day
y_max <- tapply(plot_df$log2fpkm, plot_df$day, max, na.rm = TRUE)
pval_df$y_pos <- as.numeric(y_max[as.character(pval_df$day)]) + 0.3

cat("\nDESeq2 adjusted p-values:\n")
print(pval_df[, c("day", "padj", "label")])

# ---- 5. Plot -----------------------------------------------
# Arial: windowsFonts() exists only on Windows. Register the family in a way that
# works on macOS and Linux too, and fall back to the default family if Arial is
# unavailable so the script still runs headless.
if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
  FONT_FAMILY <- "Arial"
} else {
  FONT_FAMILY <- if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else ""
}
colour_map <- c("ChAd" = "#E74C3C", "PBS" = "grey60")

p <- ggplot(plot_df, aes(x = treatment, y = log2fpkm, fill = treatment)) +

  geom_boxplot(
    outlier.shape = NA,
    alpha         = 0.6,
    width         = 0.55,
    colour        = "grey30",
    linewidth     = 0.5
  ) +

  geom_jitter(
    aes(colour = treatment),
    width = 0.12,
    size  = 2.5,
    alpha = 0.9
  ) +

  # DESeq2 p-value annotation
  geom_text(
    data        = pval_df,
    aes(x = 1.5, y = y_pos, label = label),
    inherit.aes = FALSE,
    size        = 3,
    colour      = "grey30"
  ) +

  facet_wrap(~ day, nrow = 1) +

  scale_fill_manual(values   = colour_map) +
  scale_colour_manual(values = colour_map) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +

  labs(
    title = NULL,
    x     = NULL,
    y     = "Cxcl9 log2 (FPKM + 1)"
  ) +

  theme_classic(base_size = 9, base_family = FONT_FAMILY) +
  theme(
    aspect.ratio     = 2,
    axis.title.y     = element_text(size = 9),
    axis.text        = element_text(size = 9),
    strip.text       = element_text(size = 9, face = "bold"),
    strip.background = element_rect(fill = "grey92", colour = NA),
    legend.position  = "none",
    panel.spacing    = unit(1, "lines")
  )

# ---- 6. Preview --------------------------------------------
# Only preview interactively: the default headless device has no Arial
# metrics, which aborts rendering before the file is written.
if (interactive()) {
  tryCatch(
    print(p),
    error = function(e) { grid::grid.newpage(); grid::grid.draw(ggplotGrob(p)) }
  )
}

ggsave("Fig4J_cxcl9_boxplot.pdf", plot = p, width = 9, height = 7, device = cairo_pdf)
cat("Plot saved to Fig4J_cxcl9_boxplot.pdf\n")

message("Done!")
