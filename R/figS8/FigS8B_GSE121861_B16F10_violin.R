
# Location-independent paths and fonts. See helpers/paths.R
source(file.path({
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(f) else getwd()
  while (!file.exists(file.path(d, "helpers", "paths.R")) && dirname(d) != d) d <- dirname(d)
  d
}, "helpers", "paths.R"))

# Write outputs beside this script, or wherever FIG_OUTDIR points.
setwd(out_dir())

# ============================================================
# Figure:  Figure S8B
# Title:   Cxcl10 & Cxcl9 expression per cell type in B16F10 tumours
#          (scRNA-seq, GSE121861, Kumar et al. 2018, Cell Reports)
# Method:  Subset B16F10 cells from pan-tumour Seurat object;
#          re-normalise using B16F10 cells only; violin plots
#          per annotated cell type (cell types with < 5 cells excluded)
# Input:   AllTumors_seurat.rds  (pan-tumour Seurat object, GSE121861)
#          Located in: .../Bulk RNA seq/Analysis/Decon/
# Output:  FigS8B_GSE121861_B16F10_violin.pdf  (violin plot — Fig S8B)
# ============================================================

# ---- 0. Packages -------------------------------------------
pkgs <- c("Seurat", "ggplot2", "patchwork")
installed <- pkgs %in% rownames(installed.packages())

lapply(pkgs, library, character.only = TRUE)
# working directory handled by helpers/paths.R
# ---- 1. Load pan-tumour Seurat object ----------------------
seu_all <- readRDS(file.path(data_dir("AllTumors_seurat.rds"), "AllTumors_seurat.rds"))
cat("Loaded AllTumors object:", ncol(seu_all), "cells\n")
cat("Cell types in B16F10:\n")
print(table(seu_all$cell_type[seu_all$tumor_model == "B16F10"]))

# ---- 2. Subset to B16F10 and re-normalise ------------------
seu_b16 <- subset(seu_all, subset = tumor_model == "B16F10")
seu_b16 <- JoinLayers(seu_b16)
seu_b16 <- NormalizeData(seu_b16)
cat("B16F10 cells after subset:", ncol(seu_b16), "\n")

# ---- 3. Filter rare cell types -----------------------------
keep_types <- c("Macrophages", "Inflammatory Macrophages", "Monocytes",
                "NK/T cells", "Dendritic cells", "B cells", "Tumour cells")
seu_b16 <- subset(seu_b16, subset = cell_type %in% keep_types)
cat("Cells after filtering rare types:", ncol(seu_b16), "\n")
print(table(seu_b16$cell_type))

# ---- 4. Violin plots ---------------------------------------
p <- VlnPlot(seu_b16, features = c("Cxcl10", "Cxcl9"),
             group.by = "cell_type", pt.size = 0.5)

print(p)

ggsave("FigS8B_GSE121861_B16F10_violin.pdf",
       plot = p, width = 14, height = 5)
save_source_data(p, "FigS8B_GSE121861_B16F10_violin.pdf")
cat("Plot saved to FigS8B_GSE121861_B16F10_violin.pdf\n")

message("Done!")
