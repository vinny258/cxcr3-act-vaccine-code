# ============================================================
# Figure:  Figure 5J
# Title:   Volcano plots — ChAd vs PBS at Day 1, Day 2 and Day 3
#          (bulk RNA-seq of B16-OVA tumours, mouse)
# Method:  Pre-computed DESeq2 log2FC and adjusted p-values
#          Highlights innate immune and chemokine genes
# Input:   biomarker.xlsx  (sheet: default — DESeq2 results table)
# Output:  biomarker_volcano_Day1.pdf   (Day 1 volcano)
#          biomarker_volcano_Day2.pdf   (Day 2 volcano)
#          biomarker_volcano_Day3.pdf   (Day 3 volcano)
#          biomarker_volcano_combined.pdf      (3-panel combined)
#          biomarker_volcano_combined_TEXT.pdf (text layer, for
#            overlay with biomarker_volcano_combined_DATA.png)
# ============================================================

# Load packages
packages <- c("ggplot2", "ggrepel", "patchwork", "readxl", "openxlsx", "dplyr")
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
# working directory handled by helpers/paths.R
df <- read_excel(file.path(data_dir("biomarker.xlsx"), "biomarker.xlsx"))
df <- as.data.frame(df)
colnames(df)[1] <- "ID"

# ---- 2. Day-specific genes to highlight --------------------

# Day 1: innate immune activation after adenovirus vaccination
genes_d1 <- c(
  # Innate DNA/RNA sensing
  "Cgas", "Tmem173", "Tlr9", "Ddx58", "Ifih1", "Irf3",
  # Type I Interferon response
  "Ifnb1", "Irf7", "Stat1", "Stat2", "Mx1", "Mx2", "Isg15",
  # Chemokines
  "Cxcl10", "Cxcl9", "Cxcl11", "Cxcl2", "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl7",
  # Pro-inflammatory cytokines
  "Il6", "Il1b", "Il12b", "Tnf", "Il18", "Il15", "Ifng",
  # MHC-I antigen presentation
  "H2-K1", "H2-K2", "H2-D1", "B2m", "Nlrc5", "Cd80", "Cd86"
)

# Day 2: innate-to-adaptive transition
genes_d2 <- c(
  "Igha",     # IgA antibody - early B cell response
  "C4a",      # Complement component 4a
  "Il12rb1",  # IL-12 receptor - Th1 polarisation
  "Apol6",    # Antiviral ISG
  "H19"       # Immune regulatory lncRNA
)

# Day 3: resolution phase with residual antiviral signal
genes_d3 <- c(
  "Rnasel",   # OAS-RNaseL antiviral pathway
  "Clec1a",   # C-type lectin innate receptor
  "Fmod"      # TGF-beta/inflammation resolution
)

# ---- 3. Prepare stats from pre-computed DESeq2 columns -----
# Note: positive log2FC = upregulated in ChAd
# Note: DESeq2_Pvalue is used directly as padj - no re-adjustment
prepare_stats <- function(data, pval_col, log2fc_col) {
  padj   <- as.numeric(data[[pval_col]])
  log2FC <- as.numeric(data[[log2fc_col]])

  data.frame(
    ID        = data$ID,
    Symbol    = data$gene_name,
    log2FC    = log2FC,
    padj      = padj,
    neglog10p = -log10(padj)
  )
}

# ---- 4. Thresholds -----------------------------------------
FC_THRESH <- log2(1.5)   # ~0.585, equivalent to 1.5-fold change
P_THRESH  <- 0.05

# ---- 5. Volcano plot function ------------------------------
# compact = FALSE → full-size individual plots
# compact = TRUE  → sized for 3/4 A4 landscape combined figure
make_volcano <- function(stats, title, genes_to_label, compact = FALSE) {

  bs       <- if (compact) 7   else 16    # base font size
  pt_size  <- if (compact) 0.8 else 2.0   # dot size
  lab_size <- if (compact) 2.0 else 5.7   # gene label size

  stats$group <- "NS"
  stats$group[stats$padj < P_THRESH & stats$log2FC >  FC_THRESH] <- "Up in ChAd"
  stats$group[stats$padj < P_THRESH & stats$log2FC < -FC_THRESH] <- "Up in PBS"

  counts  <- table(stats$group)
  up_chad <- ifelse("Up in ChAd" %in% names(counts), counts["Up in ChAd"], 0)
  up_pbs  <- ifelse("Up in PBS"  %in% names(counts), counts["Up in PBS"],  0)
  ns      <- ifelse("NS"         %in% names(counts), counts["NS"],         0)

  # Labels only on upregulated in ChAd genes
  stats$highlight <- tolower(stats$Symbol) %in% tolower(genes_to_label)
  label_genes     <- stats[stats$highlight & stats$group == "Up in ChAd", ]

  colour_map <- c(
    "Up in ChAd" = "#E74C3C",
    "Up in PBS"  = "#3498DB",
    "NS"         = "grey70"
  )

  ggplot(stats, aes(x = log2FC, y = neglog10p, colour = group)) +

    geom_point(alpha = 0.5, size = pt_size) +

    scale_colour_manual(
      values = colour_map,
      labels = c(
        "Up in ChAd" = paste0("Up in ChAd (", up_chad, ")"),
        "Up in PBS"  = paste0("Up in PBS (",  up_pbs,  ")"),
        "NS"         = paste0("NS (",         ns,      ")")
      )
    ) +

    geom_vline(xintercept = c(-FC_THRESH, FC_THRESH),
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = -log10(P_THRESH),
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +

    geom_text_repel(
      data               = label_genes,
      aes(label          = Symbol),
      size               = lab_size,
      fontface           = "plain",
      colour             = "black",
      box.padding        = 0.3,
      point.padding      = 0.2,
      force              = 2,
      max.overlaps       = Inf,
      min.segment.length = 0,
      segment.colour     = "black",
      segment.size       = 0.25,
      seed               = 42
    ) +

    labs(
      title    = title,
      subtitle = if (compact) NULL else
        paste0("Thresholds: |log2FC| >= ", round(FC_THRESH, 3),
               " (1.5-fold)  |  DESeq2 padj < ", P_THRESH),
      x        = "log2 Fold Change (ChAd / PBS)",
      y        = "-log10 (adjusted p-value)",
      colour   = "Expression"
    ) +
    theme_classic(base_size = bs) +
    theme(
      aspect.ratio    = 1,
      plot.title      = element_text(face = "bold", size = bs),
      axis.title      = element_text(size = bs),
      axis.text       = element_text(size = bs),
      legend.position = "none"
    )
}

# ---- 6. Compute stats for each day -------------------------
stats_d1 <- prepare_stats(df,
                          pval_col   = "PBSday1_vs_ChAdday1_DESeq2_Pvalue",
                          log2fc_col = "PBSday1_vs_ChAdday1_DESeq2_log2FC")

stats_d2 <- prepare_stats(df,
                          pval_col   = "PBSday2_vs_ChAdday2_DESeq2_Pvalue",
                          log2fc_col = "PBSday2_vs_ChAdday2_DESeq2_log2FC")

stats_d3 <- prepare_stats(df,
                          pval_col   = "PBSday3_vs_ChAdday3_DESeq2_Pvalue",
                          log2fc_col = "PBSday3_vs_ChAdday3_DESeq2_log2FC")

# ---- 7. Generate and print plots ---------------------------
p1 <- make_volcano(stats_d1, "ChAd vs PBS - Day 1", genes_d1)
p2 <- make_volcano(stats_d2, "ChAd vs PBS - Day 2", genes_d2)
p3 <- make_volcano(stats_d3, "ChAd vs PBS - Day 3", genes_d3)

print(p1)
print(p2)
print(p3)

# ---- 8. Top/Bottom 30 DEGs per timepoint -------------------
get_top_bottom <- function(stats, day_label, p_thresh = 0.05, fc_thresh = 0) {

  stats <- stats[!is.na(stats$padj), ]

  up   <- stats[stats$padj < p_thresh & stats$log2FC >  fc_thresh, ]
  down <- stats[stats$padj < p_thresh & stats$log2FC < -fc_thresh, ]

  if (nrow(up) == 0) {
    cat("\nNo upregulated genes found for", day_label, "\n")
    top30 <- data.frame()
  } else {
    up    <- up[order(-up$log2FC), ]
    top30 <- head(up, 30)
    top30$direction <- "Upregulated"
  }

  if (nrow(down) == 0) {
    cat("\nNo downregulated genes found for", day_label, "\n")
    bottom30 <- data.frame()
  } else {
    down     <- down[order(down$log2FC), ]
    bottom30 <- head(down, 30)
    bottom30$direction <- "Downregulated"
  }

  result <- rbind(top30, bottom30)
  if (nrow(result) == 0) return(NULL)

  result$day    <- day_label
  result        <- result[, c("day", "direction", "Symbol", "log2FC", "padj")]
  result$log2FC <- round(result$log2FC, 3)
  result$padj   <- signif(result$padj, 3)
  return(result)
}

top_bottom_d1 <- get_top_bottom(stats_d1, "Day 1")
top_bottom_d2 <- get_top_bottom(stats_d2, "Day 2")
top_bottom_d3 <- get_top_bottom(stats_d3, "Day 3")

if (!is.null(top_bottom_d1)) { cat("\n===== DAY 1 =====\n"); print(top_bottom_d1) }
if (!is.null(top_bottom_d2)) { cat("\n===== DAY 2 =====\n"); print(top_bottom_d2) }
if (!is.null(top_bottom_d3)) { cat("\n===== DAY 3 =====\n"); print(top_bottom_d3) }

# ---- 9. Save DEG table to Excel ----------------------------
wb <- createWorkbook()
if (!is.null(top_bottom_d1)) { addWorksheet(wb, "Day 1"); writeData(wb, "Day 1", top_bottom_d1) }
if (!is.null(top_bottom_d2)) { addWorksheet(wb, "Day 2"); writeData(wb, "Day 2", top_bottom_d2) }
if (!is.null(top_bottom_d3)) { addWorksheet(wb, "Day 3"); writeData(wb, "Day 3", top_bottom_d3) }
saveWorkbook(wb, "biomarker_top_bottom_DEGs.xlsx", overwrite = TRUE)
cat("DEG table saved to biomarker_top_bottom_DEGs.xlsx\n")

# ---- 10. Save plots as PDF ---------------------------------

# Individual full-size plots
ggsave("biomarker_volcano_Day1.pdf", plot = p1, width = 8, height = 8)
save_source_data(p1, "biomarker_volcano_Day1.pdf")
ggsave("biomarker_volcano_Day2.pdf", plot = p2, width = 8, height = 8)
save_source_data(p2, "biomarker_volcano_Day2.pdf")
ggsave("biomarker_volcano_Day3.pdf", plot = p3, width = 8, height = 8)
save_source_data(p3, "biomarker_volcano_Day3.pdf")

# Compact combined figure — 3/4 of A4 portrait width (158 x 110 mm)
p1c <- make_volcano(stats_d1, "ChAd vs PBS — Day 1", genes_d1, compact = TRUE)
p2c <- make_volcano(stats_d2, "ChAd vs PBS — Day 2", genes_d2, compact = TRUE)
p3c <- make_volcano(stats_d3, "ChAd vs PBS — Day 3", genes_d3, compact = TRUE)

p_combined <- p1c | p2c | p3c

# Preview: print each plot individually — use the arrows in the Plots pane to navigate
tryCatch(print(p1c), error = function(e) { grid::grid.newpage(); grid::grid.draw(ggplotGrob(p1c)) })
tryCatch(print(p2c), error = function(e) { grid::grid.newpage(); grid::grid.draw(ggplotGrob(p2c)) })
tryCatch(print(p3c), error = function(e) { grid::grid.newpage(); grid::grid.draw(ggplotGrob(p3c)) })

message("--- Preview shown (3 plots). Run the line below to save the combined PDF ---")
# ggsave("biomarker_volcano_combined.pdf", plot = p_combined, width = 158, height = 110, units = "mm", device = "pdf")
 save_source_data(p_combined, "biomarker_volcano_combined.pdf")

message("Done!")

# ---- 11. Illustrator-friendly export ---------------------------------------
# PNG  → data only  (dots + threshold lines, no text)   — raster layer
# PDF  → text only  (all labels + gene names, transparent background) — vector layer
# Both use identical canvas dimensions → place as layers in Illustrator

W_MM <- 158; H_MM <- 110   # matches biomarker_volcano_combined.pdf

# Data layer: keep dots and lines, hide all text
make_volcano_data_layer <- function(stats, title, compact = TRUE) {
  bs      <- if (compact) 7   else 16
  pt_size <- if (compact) 0.8 else 2.0

  stats$group <- "NS"
  stats$group[stats$padj < P_THRESH & stats$log2FC >  FC_THRESH] <- "Up in ChAd"
  stats$group[stats$padj < P_THRESH & stats$log2FC < -FC_THRESH] <- "Up in PBS"

  colour_map <- c("Up in ChAd" = "#E74C3C", "Up in PBS" = "#3498DB", "NS" = "grey70")

  ggplot(stats, aes(x = log2FC, y = neglog10p, colour = group)) +
    geom_point(alpha = 0.5, size = pt_size) +
    scale_colour_manual(values = colour_map) +
    geom_vline(xintercept = c(-FC_THRESH, FC_THRESH),
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = -log10(P_THRESH),
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    labs(title = title, x = "log2 Fold Change (ChAd / PBS)",
         y = "-log10 (adjusted p-value)") +
    theme_classic(base_size = bs) +
    theme(
      aspect.ratio    = 1,
      legend.position = "none",
      plot.title      = element_text(colour = "transparent"),
      axis.title      = element_text(colour = "transparent"),
      axis.text       = element_text(colour = "transparent"),
      axis.ticks      = element_line(colour = "transparent")
    )
}

# Text layer: invisible points (preserve coord space), hide lines, show all text
make_volcano_text_layer <- function(stats, title, genes_to_label, compact = TRUE) {
  bs       <- if (compact) 7   else 16
  pt_size  <- if (compact) 0.8 else 2.0
  lab_size <- if (compact) 2.0 else 5.7

  stats$group <- "NS"
  stats$group[stats$padj < P_THRESH & stats$log2FC >  FC_THRESH] <- "Up in ChAd"
  stats$group[stats$padj < P_THRESH & stats$log2FC < -FC_THRESH] <- "Up in PBS"

  stats$highlight <- tolower(stats$Symbol) %in% tolower(genes_to_label)
  label_genes     <- stats[stats$highlight & stats$group == "Up in ChAd", ]

  ggplot(stats, aes(x = log2FC, y = neglog10p, colour = group)) +
    geom_point(alpha = 0, size = pt_size) +   # invisible — keeps coord space identical
    scale_colour_manual(
      values = c("Up in ChAd" = "transparent", "Up in PBS" = "transparent", "NS" = "transparent")
    ) +
    geom_vline(xintercept = c(-FC_THRESH, FC_THRESH),
               linetype = "dashed", colour = "transparent", linewidth = 0.5) +
    geom_hline(yintercept = -log10(P_THRESH),
               linetype = "dashed", colour = "transparent", linewidth = 0.5) +
    geom_text_repel(
      data               = label_genes,
      aes(label          = Symbol),
      size               = lab_size,
      fontface           = "plain",
      colour             = "black",
      box.padding        = 0.3,
      point.padding      = 0.2,
      force              = 2,
      max.overlaps       = Inf,
      min.segment.length = 0,
      segment.colour     = "black",
      segment.size       = 0.25,
      seed               = 42
    ) +
    labs(title = title, x = "log2 Fold Change (ChAd / PBS)",
         y = "-log10 (adjusted p-value)") +
    theme_classic(base_size = bs) +
    theme(
      aspect.ratio     = 1,
      plot.title       = element_text(face = "bold", size = bs, colour = "black"),
      axis.title       = element_text(size = bs, colour = "black"),
      axis.text        = element_text(size = bs, colour = "black"),
      axis.ticks       = element_line(colour = "black"),
      axis.line        = element_line(colour = "transparent"),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background  = element_rect(fill = "transparent", colour = NA),
      legend.position  = "none"
    )
}

# Build and save PNG (data layer)
p1_data <- make_volcano_data_layer(stats_d1, "ChAd vs PBS \u2014 Day 1")
p2_data <- make_volcano_data_layer(stats_d2, "ChAd vs PBS \u2014 Day 2")
p3_data <- make_volcano_data_layer(stats_d3, "ChAd vs PBS \u2014 Day 3")

ggsave("biomarker_volcano_combined_DATA.png",
       plot  = p1_data | p2_data | p3_data,
       width = W_MM, height = H_MM, units = "mm",
       dpi   = 300, bg = "white")
save_source_data(p1_data | p2_data | p3_data, "biomarker_volcano_combined_DATA.png")

# Build and save PDF (text layer)
p1_txt <- make_volcano_text_layer(stats_d1, "ChAd vs PBS \u2014 Day 1", genes_d1)
p2_txt <- make_volcano_text_layer(stats_d2, "ChAd vs PBS \u2014 Day 2", genes_d2)
p3_txt <- make_volcano_text_layer(stats_d3, "ChAd vs PBS \u2014 Day 3", genes_d3)

p_txt_combined <- (p1_txt | p2_txt | p3_txt) &
  theme(plot.background  = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA))

ggsave("biomarker_volcano_combined_TEXT.pdf",
       plot   = p_txt_combined,
       width  = W_MM, height = H_MM, units = "mm",
       device = cairo_pdf, bg = "transparent")
save_source_data(p_txt_combined, "biomarker_volcano_combined_TEXT.pdf")

message("Saved: biomarker_volcano_combined_DATA.png + biomarker_volcano_combined_TEXT.pdf")
