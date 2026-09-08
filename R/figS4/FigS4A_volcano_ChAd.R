
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
# Figure:  Figure S4A
# Title:   Volcano plot — ChAdOx1 prime vs baseline
#          (human peripheral blood bulk RNA-seq, edgeR)
# Method:  Pre-computed edgeR log2FC and FDR; highlights
#          innate immune and chemokine genes (human gene symbols)
# Input:   FigS4A_volcano_ChAd.csv
#          (edgeR DE results: logFC, FDR, gene symbols)
# Output:  FigS4A_volcano_ChAd.pdf  (volcano plot)
# Highlighted genes: CXCL10, CXCL11, IRF7, DDX58, IFIH1, CCL2, CXCL9,
#                    CGAS, IL1B, IL6, CXCL2, IFNG, CCL7, IL18
# ============================================================

# ---- 0. Packages -------------------------------------------
packages <- c("ggplot2", "ggrepel", "dplyr")
installed <- packages %in% rownames(installed.packages())

lapply(packages, library, character.only = TRUE)

# ---- 1. Load data ------------------------------------------
# working directory handled by helpers/paths.R

df <- read.csv(file.path(data_dir("FigS4A_volcano_ChAd.csv"), "FigS4A_volcano_ChAd.csv"), stringsAsFactors = FALSE)
cat("Loaded:", nrow(df), "genes\n")

# ---- 2. Genes to highlight (uppercase — human convention) --
genes_highlight <- toupper(c(
  "Cxcl10", "Cxcl11", "Irf7", "Ddx58", "Ifih1", "Ccl2", "Cxcl9",
  "Cgas", "Il1b", "Il6", "Cxcl2", "Ifng", "Ccl7", "Il18"
))

# ---- 3. Prepare plotting data ------------------------------
FC_THRESH <- 0.5
P_THRESH  <- 0.05

# Replace zero FDR with smallest non-zero value
min_fdr <- min(df$FDR[df$FDR > 0], na.rm = TRUE)
cat("Replacing", sum(df$FDR == 0, na.rm = TRUE), "zero FDR values with", min_fdr, "\n")
df$FDR[df$FDR == 0] <- min_fdr

df <- df %>%
  filter(!is.na(logFC), !is.na(FDR)) %>%
  mutate(
    neglog10p    = -log10(FDR),
    group = case_when(
      FDR < P_THRESH & logFC >  FC_THRESH ~ "Up (ChAd)",
      FDR < P_THRESH & logFC < -FC_THRESH ~ "Down (ChAd)",
      TRUE                                ~ "NS"
    ),
    is_highlight = toupper(gene_name) %in% genes_highlight & group != "NS",
    label        = ifelse(is_highlight, gene_name, "")
  )

counts   <- table(df$group)
up_chad  <- ifelse("Up (ChAd)"   %in% names(counts), counts["Up (ChAd)"],   0)
dn_chad  <- ifelse("Down (ChAd)" %in% names(counts), counts["Down (ChAd)"], 0)
ns       <- ifelse("NS"          %in% names(counts), counts["NS"],          0)

cat("Up (ChAd):", up_chad, "| Down (ChAd):", dn_chad, "| NS:", ns, "\n")
cat("Highlighted genes found:", sum(df$is_highlight), "\n")
cat("Highlighted genes present:", paste(df$gene_name[df$is_highlight], collapse = ", "), "\n")

# ---- 4. Plot -----------------------------------------------
colour_map <- c(
  "Up (ChAd)"   = "#E74C3C",
  "Down (ChAd)" = "#3498DB",
  "NS"          = "grey70"
)

df_bg <- df[!df$is_highlight, ]
df_fg <- df[ df$is_highlight, ]

p <- ggplot() +

  # Background: all non-highlighted genes
  geom_point(
    data  = df_bg,
    aes(x = logFC, y = neglog10p, colour = group),
    size  = 1.2, alpha = 0.4
  ) +

  # Foreground: highlighted genes (on top, larger)
  geom_point(
    data   = df_fg,
    aes(x  = logFC, y = neglog10p),
    shape  = 21, fill = "#E74C3C", colour = "black",
    size   = 3.0, stroke = 0.6, alpha = 0.95
  ) +

  # Threshold lines
  geom_vline(xintercept = c(-FC_THRESH, FC_THRESH),
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(P_THRESH),
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +

  # Labels for highlighted genes
  geom_text_repel(
    data               = df_fg,
    aes(x              = logFC, y = neglog10p, label = gene_name),
    size               = 4.5,
    fontface           = "plain",
    colour             = "black",
    box.padding        = 0.4,
    point.padding      = 0.3,
    force              = 3,
    max.overlaps       = Inf,
    min.segment.length = 0,
    segment.colour     = "black",
    segment.size       = 0.3,
    seed               = 42
  ) +

  scale_colour_manual(
    values = colour_map,
    labels = c(
      "Up (ChAd)"   = paste0("Up in ChAd (", up_chad, ")"),
      "Down (ChAd)" = paste0("Down in ChAd (", dn_chad, ")"),
      "NS"          = paste0("NS (", ns, ")")
    )
  ) +

  labs(
    title    = "ChAdOx1 prime — human blood bulk RNA-seq",
    subtitle = paste0("Thresholds: |logFC| > ", FC_THRESH, "  |  FDR < ", P_THRESH),
    x        = "log2 Fold Change (ChAdOx1 / baseline)",
    y        = "-log10 (FDR)",
    colour   = "Expression"
  ) +

  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", size = 14),
    axis.title   = element_text(size = 13),
    axis.text    = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11),
    aspect.ratio = 1
  )

# ---- 5. Preview in VS Code ---------------------------------
tryCatch(
  print(p),
  error = function(e) { grid::grid.newpage(); grid::grid.draw(ggplotGrob(p)) }
)

ggsave("FigS4A_volcano_ChAd.pdf", plot = p, width = 8, height = 8, device = "pdf")
save_source_data(p, "FigS4A_volcano_ChAd.pdf")
cat("Plot saved to FigS4A_volcano_ChAd.pdf\n")

message("Done!")
