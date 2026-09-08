# ============================================================
# Figure:  Figure 3C
# Title:   Pseudobulk differential expression - Day 1 vs Day 0
#          (WT draining lymph node cells, scRNA-seq)
# Method:  Pseudobulk (raw RNA counts) per mouse x timepoint -> DESeq2
#          Design: ~ day  (UNPAIRED: day 0 and day 1 are independent
#          cohorts of mice, n = 5 per timepoint; lymph-node harvest is
#          terminal, so no mouse is shared across timepoints).
#          Mouse identity = cellhashR consensus call ("consensuscall";
#          negatives/doublets already removed).
#          Gene filter: expressed in >= 5% of cells in either timepoint
#          AND >= 10 counts in >= 2 pseudobulk samples (removes genes
#          driven by very few cells, e.g. Lcn2, Cxcl2).
# Input:   shared data/wt_D0_D1.rds  (Seurat object)
# Output:  Fig3C_pseudobulk_volcano_d1_vs_d0.pdf  (volcano plot)
#          pseudobulk_de_d1_vs_d0.csv             (DE results table)
#          FigS4F_pseudobulk_GSEA.csv             (ranked input for fig S4F)
# Refs:    Seurat (Hao 2021); DESeq2 (Love 2014); cellhashR (McGinnis 2019)
# ============================================================

library(Seurat); library(DESeq2); library(Matrix)
library(ggplot2); library(ggrepel); library(dplyr)

# Location-independent paths and fonts. See helpers/paths.R
source(file.path(rprojroot_find <- {
  d <- tryCatch(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])),
                error = function(e) getwd())
  if (!length(d) || !nzchar(d)) d <- getwd()
  while (!file.exists(file.path(d, "helpers", "paths.R")) && dirname(d) != d) d <- dirname(d)
  d
}, "helpers", "paths.R"))


OUTDIR <- out_dir()
RDS <- file.path(data_dir("wt_D0_D1.rds"), "wt_D0_D1.rds")
PCT_THRESHOLD <- 0.05     # min fraction of cells expressing a gene (either timepoint)
FC_THRESH <- 0.5; P_THRESH <- 0.05

# ---- 1. Load WT singlets; mouse = consensuscall ------------
seu <- readRDS(RDS)
seu_wt <- subset(seu, group == "WT" & !is.na(consensuscall))
DefaultAssay(seu_wt) <- "RNA"
md <- seu_wt@meta.data
md$consensuscall <- as.character(md$consensuscall); md$day <- as.character(md$day)
cat("WT singlet cells:", ncol(seu_wt), "\n")
cat("Cells per mouse x day:\n"); print(table(md$consensuscall, md$day))

cnts <- GetAssayData(seu_wt, assay = "RNA", layer = "counts")

# ---- 2. Pseudobulk per mouse x day (raw counts) ------------
grp <- paste(md$consensuscall, md$day, sep = "_")
samples <- sort(unique(grp))
pb <- sapply(samples, function(s) Matrix::rowSums(cnts[, grp == s, drop = FALSE]))
rownames(pb) <- rownames(cnts)
cat("\nPseudobulk:", nrow(pb), "genes x", ncol(pb), "samples\n")

coldata <- data.frame(
  row.names = colnames(pb),
  mouse = sub("_D[01]$", "", colnames(pb)),
  day   = relevel(factor(sub(".*_(D[01])$", "\\1", colnames(pb))), ref = "D0")
)

# ---- 3. Gene filters ---------------------------------------
expr   <- cnts > 0
pct_D0 <- Matrix::rowMeans(expr[, md$day == "D0", drop = FALSE])
pct_D1 <- Matrix::rowMeans(expr[, md$day == "D1", drop = FALSE])
keep_pct <- (pct_D0 >= PCT_THRESHOLD) | (pct_D1 >= PCT_THRESHOLD)
keep_cnt <- rowSums(pb >= 10) >= 2
keep <- keep_pct & keep_cnt
cat(sprintf("Genes kept: %d (>= %.0f%% expressing in either group AND >=10 counts in >=2 samples)\n",
            sum(keep), 100 * PCT_THRESHOLD))

# ---- 4. DESeq2 (design ~ day) ------------------------------
dds <- DESeqDataSetFromMatrix(round(pb[keep, ]), coldata, design = ~ day)
dds <- DESeq(dds)
res <- results(dds, contrast = c("day", "D1", "D0"), alpha = P_THRESH)
cat("\nDE summary (D1 vs D0):\n"); summary(res)

df <- as.data.frame(res) %>% tibble::rownames_to_column("gene") %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>% arrange(padj)
write.csv(df, file.path(OUTDIR, "pseudobulk_de_d1_vs_d0.csv"), row.names = FALSE)
# ranked input for fig S4F GSEA (needs gene, log2FoldChange, stat, padj)
write.csv(df[, c("gene","log2FoldChange","stat","padj")],
          file.path(OUTDIR, "FigS4F_pseudobulk_GSEA.csv"), row.names = FALSE)

cat("\nUp (log2FC>0.5, padj<0.05):", sum(df$log2FoldChange > FC_THRESH & df$padj < P_THRESH),
    "| Down:", sum(df$log2FoldChange < -FC_THRESH & df$padj < P_THRESH), "\n")

# ---- 5. Volcano --------------------------------------------
genes_highlight <- c("Cxcl10","Irf7","Ddx58","Ifih1","Ccl2","Cxcl9",
                     "Cgas","Il1b","Il6","Cxcl2","Ifng","Ccl7","Il18")
min_p <- min(df$padj[df$padj > 0], na.rm = TRUE)
df <- df %>% mutate(
  neglog10p = -log10(pmax(padj, min_p)),
  group = case_when(padj < P_THRESH & log2FoldChange >  FC_THRESH ~ "Up in D1",
                    padj < P_THRESH & log2FoldChange < -FC_THRESH ~ "Up in D0",
                    TRUE ~ "NS"),
  is_highlight = gene %in% genes_highlight)
present <- intersect(genes_highlight, df$gene)
cat("Highlighted genes present after filtering:", paste(present, collapse=", "), "\n")
cat("Highlighted genes removed by filter:", paste(setdiff(genes_highlight, present), collapse=", "), "\n")

nU <- sum(df$group=="Up in D1"); nD <- sum(df$group=="Up in D0"); nN <- sum(df$group=="NS")
cmap <- c("Up in D1"="#E74C3C","Up in D0"="#3498DB","NS"="grey70")
p <- ggplot() +
  geom_point(data=df[!df$is_highlight,], aes(log2FoldChange, neglog10p, colour=group), size=1.2, alpha=0.4) +
  geom_point(data=df[df$is_highlight,], aes(log2FoldChange, neglog10p), shape=21, fill="#E74C3C",
             colour="black", size=3, stroke=0.6) +
  geom_vline(xintercept=c(-FC_THRESH,FC_THRESH), linetype="dashed", colour="grey40") +
  geom_hline(yintercept=-log10(P_THRESH), linetype="dashed", colour="grey40") +
  geom_text_repel(data=df[df$is_highlight,], aes(log2FoldChange, neglog10p, label=gene),
                  size=4.5, colour="black", box.padding=0.4, force=3, max.overlaps=Inf,
                  min.segment.length=0, segment.size=0.3, seed=42, fontface="italic") +
  scale_colour_manual(values=cmap, labels=c("Up in D1"=paste0("Up in D1 (",nU,")"),
                      "Up in D0"=paste0("Up in D0 (",nD,")"), "NS"=paste0("NS (",nN,")"))) +
  labs(title="Day 1 vs Day 0 - pseudobulk DESeq2 (~ day, n=5 mice/timepoint)",
       subtitle=paste0("Filter: >=",100*PCT_THRESHOLD,"% expressing | |log2FC|>",FC_THRESH," | padj<",P_THRESH),
       x="log2 Fold Change (D1 / D0)", y="-log10 (adjusted p-value)", colour="Expression") +
  theme_classic(base_size=14) + theme(plot.title=element_text(face="bold", size=13), aspect.ratio=1)

ggsave(file.path(OUTDIR,"Fig3C_pseudobulk_volcano_d1_vs_d0.pdf"), p, width=8, height=8)
cat("Saved volcano + DE tables to", OUTDIR, "\n")
cat("FIG3C_DONE\n")
