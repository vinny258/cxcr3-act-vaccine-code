# ============================================================
# Figure:  Figure S4F
# Title:   GSEA - upregulated pathways at Day 1 vs Day 0
#          (WT draining lymph node cells, scRNA-seq pseudobulk)
# Method:  Genes ranked by DESeq2 Wald statistic; fgsea.
#          Gene sets: MSigDB Hallmark + Reactome (Mus musculus).
#          DOWNSTREAM of Fig3C: uses the ~ day (unpaired) pseudobulk
#          DESeq2 ranking. Re-run Fig3C first if the design changes.
# Input:   ../Fig3C/FigS4F_pseudobulk_GSEA.csv (gene, log2FoldChange, stat, padj)
# Output:  FigS4F_pseudobulk_GSEA.pdf ; pseudobulk_gsea_results_upregulated.csv
# Refs:    fgsea (Korotkevich 2021); msigdbr; DESeq2 (Love 2014)
# ============================================================
library(fgsea); library(msigdbr); library(dplyr); library(ggplot2); library(tibble)

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


BASE   <- data_dir("the upstream analysis outputs")
INPUT  <- file.path(BASE, "Fig3C", "FigS4F_pseudobulk_GSEA.csv")  # produced by Fig3C (~ day)
OUTDIR <- file.path(BASE, "FigS4F")

de <- read.csv(INPUT, stringsAsFactors = FALSE)
stopifnot(all(c("gene","log2FoldChange","stat","padj") %in% colnames(de)))
cat("Loaded ranked DE genes:", nrow(de), "from", INPUT, "\n")

de <- de %>% filter(!is.na(stat), !is.na(gene), gene != "") %>%
  group_by(gene) %>% slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>% ungroup()
ranks <- sort(deframe(select(de, gene, stat)), decreasing = TRUE)
cat("Genes in ranked list:", length(ranks), "| range", round(min(ranks),2), "to", round(max(ranks),2), "\n")

h_sets <- msigdbr(species="Mus musculus", collection="H") %>%
  select(gs_name, gene_symbol) %>% split(f=.$gs_name) %>% lapply(function(x) x$gene_symbol)
reactome_sets <- msigdbr(species="Mus musculus", collection="C2", subcollection="CP:REACTOME") %>%
  select(gs_name, gene_symbol) %>% split(f=.$gs_name) %>% lapply(function(x) x$gene_symbol)
all_sets <- c(h_sets, reactome_sets)
cat("Gene sets - Hallmark:", length(h_sets), "| Reactome:", length(reactome_sets), "\n")

set.seed(42)
fgsea_res <- fgsea(pathways=all_sets, stats=ranks, minSize=15, maxSize=500, nPermSimple=10000) %>%
  arrange(desc(NES)) %>%
  mutate(source = ifelse(grepl("^HALLMARK_", pathway), "Hallmark", "Reactome"),
         pathway_clean = pathway %>% gsub("^HALLMARK_|^REACTOME_","",.) %>% gsub("_"," ",.) %>%
           tolower() %>% tools::toTitleCase())
cat("\nSignificant upregulated pathways (padj<0.05, NES>0):",
    sum(fgsea_res$padj<0.05 & fgsea_res$NES>0, na.rm=TRUE), "\n")

hm  <- fgsea_res %>% filter(padj<0.05, NES>0, source=="Hallmark") %>% slice_max(NES, n=10) %>% pull(pathway)
rea <- fgsea_res %>% filter(padj<0.05, NES>0, source=="Reactome") %>% slice_max(NES, n=10) %>% pull(pathway)
selected <- c(hm, rea)

top_paths <- fgsea_res %>% filter(pathway %in% selected) %>%
  mutate(source = factor(source, levels=c("Hallmark","Reactome")),
         pathway_nice = pathway_clean %>%
           gsub("Interferon Gamma Response","IFNγ Response",.) %>%
           gsub("Interferon Alpha Response","IFNα Response",.) %>%
           gsub("Interferon Alpha Beta Signaling","IFNα/β Signaling",.) %>%
           gsub("Tnfa Signaling Via Nfkb","TNFα Signaling via NFκB",.,ignore.case=TRUE) %>%
           gsub("Dna","DNA",.) %>% gsub("Rna","RNA",.) %>%
           gsub("Mhc Class I([^I])","MHC Class I\\1",.) %>% gsub("Mhc Class Ii","MHC Class II",.) %>%
           gsub("Isg15","ISG15",.) %>% gsub("Irf([0-9])","IRF\\1",.)) %>%
  arrange(source, NES) %>% mutate(pathway_nice = factor(pathway_nice, levels=unique(pathway_nice)))
cat("Pathways in plot - Hallmark:", length(hm), "| Reactome:", length(rea), "\n")

fgsea_res %>% filter(padj<0.05, NES>0) %>% arrange(desc(NES)) %>%
  mutate(in_plot = pathway %in% selected, leadingEdge = sapply(leadingEdge, paste, collapse=";")) %>%
  write.csv("pseudobulk_gsea_results_upregulated.csv", row.names=FALSE)

p <- ggplot(top_paths, aes(NES, pathway_nice, size=size, colour=-log10(padj))) +
  geom_point(alpha=0.85) +
  scale_colour_gradient(low="#F0997B", high="#C85A30", name="-log10(padj)") +
  scale_size_continuous(name="Gene set size", range=c(3,10)) +
  facet_grid(source ~ ., scales="free_y", space="free_y") +
  labs(title="Upregulated pathways: D1 vaccinated vs D0 unvaccinated",
       subtitle="Pseudobulk DESeq2 (~ day) | Hallmark + Reactome | padj < 0.05",
       x="Normalized Enrichment Score (NES)", y=NULL) +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(size=14, face="bold"), axis.text.y=element_text(size=13),
        strip.text=element_text(size=12, face="bold"),
        strip.background=element_rect(fill="grey92", colour=NA), panel.spacing=unit(1,"lines"))
ggsave("FigS4F_pseudobulk_GSEA.pdf", p, width=11, height=max(6, nrow(top_paths)*0.4+2), device=cairo_pdf)
save_source_data(p, "FigS4F_pseudobulk_GSEA.pdf")
cat("Saved FigS4F outputs to", OUTDIR, "\nFIGS4F_DONE\n")
