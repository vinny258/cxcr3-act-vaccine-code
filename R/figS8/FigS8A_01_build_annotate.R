# ============================================================
# Figure S8A, step 1 of 2: build and annotate
#
# Loads the three GSE307143 B16F10 samples, does QC, clustering and manual
# annotation, and saves the annotated Seurat object. The panel itself is drawn
# by FigS8A_02_plot_violin.R from that object.
#
# WHY THIS IS SPLIT
# Annotation maps *cluster numbers* to cell type names. Cluster numbering is not
# stable across Seurat versions: re-running under Seurat 5.5.0 on 2026-09-07
# produced 16 clusters where the original run produced 17, so the label for
# cluster "16" (Neutrophils) was silently never applied and that population
# vanished from the figure with no error. The check below now makes that a hard
# failure rather than a silent one.
#
# Pin Seurat to the version recorded in the renv lockfile before running this.
#
# Input :  GSM9217284/5/6 matrix, features and barcodes files (this folder)
# Output:  B16F10_GSE307143_annotated.rds  plus QC and UMAP figures
# ============================================================


# ---- 0. Packages -------------------------------------------
pkgs <- c("Seurat", "ggplot2", "dplyr", "ggrepel", "patchwork", "scales")
installed <- pkgs %in% rownames(installed.packages())
if (any(!installed)) install.packages(pkgs[!installed])
lapply(pkgs, library, character.only = TRUE)

script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(normalizePath(dirname(f)))
  if (!is.null(sys.frames()[[1]]$ofile)) return(normalizePath(dirname(sys.frames()[[1]]$ofile)))
  getwd()
})
setwd(script_dir)
output_dir <- Sys.getenv("FIG_OUTDIR", unset = script_dir)
# ---- 1. Load each sample with ReadMtx ----------------------
# Sample 1 barcodes file has an unusual name (_. instead of _barcodes.)
# ReadMtx() lets us specify exact paths regardless of filename

load_sample <- function(prefix, barcodes_file, features_file, matrix_file, sample_name) {
  cat("Loading", sample_name, "...\n")
  mat <- ReadMtx(
    mtx      = matrix_file,
    cells    = barcodes_file,
    features = features_file
  )
  seu <- CreateSeuratObject(counts = mat, project = sample_name,
                             min.cells = 3, min.features = 100)
  seu$sample <- sample_name
  cat("  Cells loaded:", ncol(seu), "\n")
  seu
}

s1 <- load_sample(
  sample_name    = "B16F10_1",
  barcodes_file  = "GSM9217284_SS_20_0527_.tsv.gz",
  features_file  = "GSM9217284_SS_20_0527_features.tsv.gz",
  matrix_file    = "GSM9217284_SS_20_0527_matrix.mtx.gz"
)

s2 <- load_sample(
  sample_name    = "B16F10_2",
  barcodes_file  = "GSM9217285_SS_20_0528_barcodes.tsv.gz",
  features_file  = "GSM9217285_SS_20_0528_features.tsv.gz",
  matrix_file    = "GSM9217285_SS_20_0528_matrix.mtx.gz"
)

s3 <- load_sample(
  sample_name    = "B16F10_3",
  barcodes_file  = "GSM9217286_SS_20_0529_barcodes.tsv.gz",
  features_file  = "GSM9217286_SS_20_0529_features.tsv.gz",
  matrix_file    = "GSM9217286_SS_20_0529_matrix.mtx.gz"
)

# ---- 2. Merge and QC ---------------------------------------
seu <- merge(s1, y = list(s2, s3),
             add.cell.ids = c("B16F10_1", "B16F10_2", "B16F10_3"),
             project = "B16F10_GSE307143")
rm(s1, s2, s3); gc()

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

cat("\nCells per sample before QC:\n")
print(table(seu$sample))

seu <- subset(seu,
              subset = nFeature_RNA > 200 &
                       nFeature_RNA < 7000 &
                       percent.mt  < 25)

cat("\nCells per sample after QC:\n")
print(table(seu$sample))
cat("Total cells:", ncol(seu), "\n")

# ---- 3. Normalise and cluster ------------------------------
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.4)
seu <- RunUMAP(seu, dims = 1:20, seed.use = 42)

cat("\nClusters found:", length(levels(seu$seurat_clusters)), "\n")
print(table(seu$seurat_clusters))

# ---- 4. Find markers and print for annotation --------------
seu <- JoinLayers(seu)

cat("\nFinding cluster markers (for annotation)...\n")
cluster_markers <- FindAllMarkers(seu, only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.5,
                                  verbose = FALSE)
top5 <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  select(cluster, gene, avg_log2FC, pct.1, p_val_adj)

cat("\nTop 5 markers per cluster:\n")
print(as.data.frame(top5))

cat("\nCxcl10 mean expression per cluster:\n")
cxcl10_expr <- AverageExpression(seu, features = "Cxcl10",
                                  group.by = "seurat_clusters")
print(cxcl10_expr$RNA)

# ---- 5. UMAP and dot plot for annotation -------------------
markers_known <- list(
  "Tumour cells"    = c("Dct", "Tyrp1", "Mlana", "Mitf", "Sox10", "Pmel"),
  "Macrophage"      = c("Csf1r", "Adgre1", "Cd68", "C1qa", "C1qb"),
  "Inf. Macrophage" = c("Csf1r", "Adgre1", "Cd68", "Nos2", "Il1b", "Cxcl10"),
  "DC"              = c("Itgax", "H2-Eb1", "Xcr1", "Siglech", "Cd209a"),
  "T/NK cell"       = c("Cd3e", "Cd3d", "Cd8a", "Cd4", "Nkg7", "Klrd1"),
  "B cell"          = c("Cd79a", "Ms4a1", "Cd19"),
  "Monocyte"        = c("Ly6c2", "Ccr2", "S100a8", "S100a9"),
  "Endothelial"     = c("Pecam1", "Cdh5", "Kdr"),
  "Pericyte"        = c("Pdgfrb", "Rgs5", "Acta2"),
  "Fibroblast"      = c("Col1a1", "Col1a2", "Pdgfra")
)

all_markers <- unique(unlist(markers_known))
all_markers <- all_markers[all_markers %in% rownames(seu)]

p_umap <- DimPlot(seu, reduction = "umap", label = TRUE,
                  label.size = 4, pt.size = 0.3) +
  labs(title = "B16F10 — Seurat clusters (GSE307143)",
       subtitle = paste0(ncol(seu), " cells")) +
  theme(plot.title = element_text(face = "bold"))

p_umap_sample <- DimPlot(seu, reduction = "umap",
                          group.by = "sample", pt.size = 0.3) +
  labs(title = "Coloured by replicate") +
  theme(plot.title = element_text(face = "bold"))

p_dot <- DotPlot(seu, features = all_markers) +
  RotatedAxis() +
  labs(title = "Canonical markers per cluster", x = NULL, y = "Cluster") +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 8))

p_cxcl10 <- FeaturePlot(seu, features = "Cxcl10", reduction = "umap",
                         pt.size = 0.3, cols = c("lightgrey", "#C0392B")) +
  labs(title = "Cxcl10 expression") +
  theme(plot.title = element_text(face = "bold"))

print(p_umap)
print(p_umap_sample)
print(p_dot)
print(p_cxcl10)

ggsave("B16F10_UMAP_clusters.pdf",    plot = p_umap,        width = 7, height = 6)
ggsave("B16F10_UMAP_sample.pdf",      plot = p_umap_sample, width = 7, height = 6)
ggsave("B16F10_dotplot_markers.pdf",  plot = p_dot,         width = 14, height = 6)
ggsave("B16F10_Cxcl10_UMAP.pdf",     plot = p_cxcl10,      width = 7, height = 6)

message("--- Check the plots and top markers above, then annotate clusters in SECTION 6 ---")

# ---- 6. ANNOTATE CLUSTERS AND REMOVE LOW-QUALITY CELLS ----

# Remove cluster 8 (low-quality: high mt%, low nFeature, Gm42418-driven)
seu <- subset(seu, idents = "8", invert = TRUE)

new_labels <- c(
  "0"  = "Macrophages",
  "1"  = "Monocytes",
  "2"  = "NKT cells",
  "3"  = "Dendritic cells",
  "4"  = "NK cells",
  "5"  = "Proliferating cells",
  "6"  = "Inflammatory Macrophages",
  "7"  = "T cells",
  "9"  = "CD4+ T cells",
  "10" = "B cells",
  "11" = "Activated DCs",
  "12" = "CD8+ T cells",
  "13" = "Regulatory T cells",
  "14" = "Immature lymphocytes",
  "15" = "Tissue-resident Macrophages",
  "16" = "Neutrophils"
)


# ---- guard: the cluster to label map must apply completely -------------
# If Seurat returns a different number of clusters, some labels silently go
# unused and whole populations disappear from the figure. Fail loudly instead.
found_clusters <- levels(seu$seurat_clusters)
expected <- names(new_labels)
unmatched <- setdiff(expected, found_clusters)
if (length(unmatched)) {
  stop("Cluster/label mismatch. Clusters expected by the label map but not found: ",
       paste(unmatched, collapse = ", "),
       "\nFound ", length(found_clusters), " clusters: ",
       paste(found_clusters, collapse = ", "),
       "\nThis means the Seurat version differs from the one used for the paper. ",
       "Pin Seurat via renv and re-run, or use the archived ",
       "B16F10_GSE307143_annotated.rds instead.")
}
cat("Cluster/label map applied cleanly:", length(expected), "labels\n")

seu <- RenameIdents(seu, new_labels)
seu$cell_type <- Idents(seu)

cat("\nCell counts per type:\n")
print(sort(table(seu$cell_type), decreasing = TRUE))

# Annotated UMAP
p_umap_annotated <- DimPlot(seu, reduction = "umap", label = TRUE,
                             label.size = 3.5, pt.size = 0.3, repel = TRUE) +
  labs(title = "B16F10 — Annotated cell types (GSE307143)",
       subtitle = paste0(ncol(seu), " cells")) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none")

print(p_umap_annotated)
ggsave("B16F10_UMAP_annotated.pdf", plot = p_umap_annotated, width = 8, height = 7)

# ---- 7. Cxcl10 dot plot per cell type ----------------------

cxcl10_df <- FetchData(seu, vars = c("Cxcl10", "cell_type")) %>%
  group_by(cell_type) %>%
  summarise(
    pct_expressing  = mean(Cxcl10 > 0) * 100,
    mean_expression = mean(Cxcl10[Cxcl10 > 0], na.rm = TRUE)
  ) %>%
  arrange(pct_expressing)

cat("\nCxcl10 expression per cell type:\n")
print(as.data.frame(cxcl10_df))

p_cxcl10_celltype <- ggplot(cxcl10_df,
  aes(x = pct_expressing, y = reorder(cell_type, pct_expressing),
      size = pct_expressing, colour = mean_expression)) +
  geom_point() +
  scale_colour_gradient(low = "floralwhite", high = "#C0392B",
                        name = "Mean expression\n(expressors only)") +
  scale_size_continuous(name = "% expressing",
                        range = c(1, 10),
                        breaks = c(10, 30, 50),
                        labels = c("10%", "30%", "50%")) +
  labs(title = "Cxcl10 expression in B16F10",
       subtitle = "GSE307143 — 3 biological replicates",
       x = "% cells expressing Cxcl10", y = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

print(p_cxcl10_celltype)
dir.create(file.path(output_dir, "FigS8A"), showWarnings = FALSE)
ggsave(file.path(output_dir, "FigS8A", "Cxcl10_B16F10_celltype_GSE307143.pdf"),
       plot = p_cxcl10_celltype, width = 8, height = 5)

# ---- 8. Cxcl9 dot plot per cell type (Fig S8C) -------------

cxcl9_df <- FetchData(seu, vars = c("Cxcl9", "cell_type")) %>%
  group_by(cell_type) %>%
  summarise(
    pct_expressing  = mean(Cxcl9 > 0) * 100,
    mean_expression = mean(Cxcl9[Cxcl9 > 0], na.rm = TRUE)
  ) %>%
  arrange(pct_expressing)

cat("\nCxcl9 expression per cell type:\n")
print(as.data.frame(cxcl9_df))

p_cxcl9_celltype <- ggplot(cxcl9_df,
  aes(x = pct_expressing, y = reorder(cell_type, pct_expressing),
      size = pct_expressing, colour = mean_expression)) +
  geom_point() +
  scale_colour_gradient(low = "floralwhite", high = "#C0392B",
                        name = "Mean expression\n(expressors only)") +
  scale_size_continuous(name = "% expressing",
                        range = c(1, 10),
                        breaks = c(10, 30, 50),
                        labels = c("10%", "30%", "50%")) +
  labs(title = "Cxcl9 expression in B16F10",
       subtitle = "GSE307143 — 3 biological replicates",
       x = "% cells expressing Cxcl9", y = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

print(p_cxcl9_celltype)
ggsave(file.path(output_dir, "FigS8A", "Cxcl9_B16F10_celltype_GSE307143.pdf"),
       plot = p_cxcl9_celltype, width = 8, height = 5)

# ---- 9. Il21r dot plot per cell type -----------------------

il21r_df <- FetchData(seu, vars = c("Il21r", "cell_type")) %>%
  group_by(cell_type) %>%
  summarise(
    pct_expressing  = mean(Il21r > 0) * 100,
    mean_expression = mean(Il21r[Il21r > 0], na.rm = TRUE)
  ) %>%
  arrange(pct_expressing)

cat("\nIl21r expression per cell type:\n")
print(as.data.frame(il21r_df))

p_il21r_celltype <- ggplot(il21r_df,
  aes(x = pct_expressing, y = reorder(cell_type, pct_expressing),
      size = pct_expressing, colour = mean_expression)) +
  geom_point() +
  scale_colour_gradient(low = "floralwhite", high = "#C0392B",
                        name = "Mean expression\n(expressors only)") +
  scale_size_continuous(name = "% expressing",
                        range = c(1, 10),
                        breaks = c(10, 30, 50),
                        labels = c("10%", "30%", "50%")) +
  labs(title = "Il21r expression in B16F10",
       subtitle = "GSE307143 — 3 biological replicates",
       x = "% cells expressing Il21r", y = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

print(p_il21r_celltype)
ggsave(file.path(output_dir, "FigS8A", "Il21r_B16F10_celltype_GSE307143.pdf"),
       plot = p_il21r_celltype, width = 8, height = 5)

# ---- 10. Violin plots — Cxcl10, Cxcl9 & Il21r per cell type

# ---- save the annotated object; step 2 draws the panel ----------------
saveRDS(seu, file.path(output_dir, "B16F10_GSE307143_annotated.rds"))
cat("Saved B16F10_GSE307143_annotated.rds\n")
cat("Now run FigS8A_02_plot_violin.R to draw the panel.\n")
