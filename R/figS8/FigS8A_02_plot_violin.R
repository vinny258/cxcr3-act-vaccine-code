# ============================================================
# Figure S8A, step 2 of 2: draw the panel
#
# Reads the annotated Seurat object from step 1 and draws the published panel.
# Deterministic: no clustering, so the same figure on any Seurat version.
#
# The panel is a horizontal, jittered, faceted violin. Seurat's VlnPlot() does
# not produce this layout, so it is built directly in ggplot2. The ordering,
# the "Proliferating T/NK cells" label and the choice of two features were all
# applied by hand in Illustrator for the paper; they are set explicitly here so
# the script alone reproduces the published panel.
#
# Input :  B16F10_GSE307143_annotated.rds
# Output:  FigS8A_GSE307143_violin.pdf
# ============================================================

suppressMessages({ library(Seurat); library(ggplot2); library(dplyr); library(tidyr); library(scales) })

script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(normalizePath(dirname(f)))
  if (!is.null(sys.frames()[[1]]$ofile)) return(normalizePath(dirname(sys.frames()[[1]]$ofile)))
  getwd()
})
OUTDIR <- Sys.getenv("FIG_OUTDIR", unset = script_dir)
RDS    <- Sys.getenv("FIGS8A_RDS",
                     unset = file.path(script_dir, "B16F10_GSE307143_annotated.rds"))

if (!file.exists(RDS))
  stop("Annotated object not found: ", RDS, "\nRun FigS8A_01_build_annotate.R first.")

seu <- readRDS(RDS)
GENES <- c("Cxcl10", "Cxcl9")           # Il21r was in the original call, not in the panel

PANEL_ORDER <- c(
  "Inflammatory Macrophages", "Macrophages", "Monocytes",
  "Regulatory T cells", "Dendritic cells", "B cells",
  "Activated DCs", "Neutrophils", "Tissue-resident Macrophages",
  "Immature lymphocytes", "CD8+ T cells", "CD4+ T cells",
  "T cells", "NK cells", "NKT cells", "Proliferating T/NK cells"
)

ct <- as.character(seu$cell_type)
ct[ct == "Proliferating cells"] <- "Proliferating T/NK cells"   # relabelled for the figure

missing <- setdiff(PANEL_ORDER, unique(ct))
if (length(missing))
  stop("Cell types in the published panel are absent from this object: ",
       paste(missing, collapse = ", "),
       "\nThe object was probably built with a different Seurat version. ",
       "See 00_README_FigS8A.md.")

# Colours follow the ORIGINAL annotation order, not the panel order. In the paper
# the rows were reordered by hand after plotting, so each cell type kept the colour
# ggplot had assigned it from the annotation order. Reproduce that mapping.
orig_levels <- levels(factor(seu$cell_type))
orig_levels[orig_levels == "Proliferating cells"] <- "Proliferating T/NK cells"
PALETTE <- setNames(scales::hue_pal()(length(orig_levels)), orig_levels)

expr <- FetchData(seu, vars = GENES, layer = "data")
df <- data.frame(cell_type = ct, expr, check.names = FALSE) |>
  pivot_longer(all_of(GENES), names_to = "gene", values_to = "expression") |>
  mutate(cell_type = factor(cell_type, levels = rev(PANEL_ORDER)),
         gene      = factor(gene, levels = GENES))

p <- ggplot(df, aes(x = expression, y = cell_type, fill = cell_type)) +
  geom_violin(scale = "width", linewidth = 0.3, colour = "grey20") +
  geom_jitter(height = 0.25, size = 0.05, alpha = 0.35, colour = "grey20") +
  scale_fill_manual(values = PALETTE) +
  facet_wrap(~ gene, nrow = 1) +
  labs(x = "Normalized expression", y = NULL,
       title = "B16F10 tumors from GSE307143, Wang et al. 2025") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none",
        panel.grid      = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black"),
        plot.title      = element_text(size = 8, hjust = 0.5))

ggsave(file.path(OUTDIR, "FigS8A_GSE307143_violin.pdf"), p, width = 9, height = 5)
save_source_data(p, file.path(OUTDIR, "FigS8A_GSE307143_violin.pdf"))
cat("Wrote", file.path(OUTDIR, "FigS8A_GSE307143_violin.pdf"), "\n")
cat("Cell types drawn:", nlevels(df$cell_type), "| genes:", paste(GENES, collapse = ", "), "\n")
