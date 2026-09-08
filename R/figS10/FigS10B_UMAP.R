## =====================================================================
## Exp 315 - STEP 11: UMAP of CD3+CD8+ T cells
## Embedding on phenotypic + functional markers (CD44, CD62L, CD127, KLRG1,
## IFNg, TNF, IL2, CD107a). Coloured by donor/host (CD45.1), genotype/group,
## and per-marker expression. Cells from OVA(SIINFEKL)-stimulated tubes of the
## four vaccinated transfer groups.
## Outputs (R analysis/figures/):
##   Fig8a_UMAP_donor_by_group.png
##   Fig8b_UMAP_markers.png
##   Fig8c_UMAP_genotype.png
## =====================================================================
suppressMessages({library(flowCore); library(dplyr); library(tidyr)
                  library(ggplot2); library(uwot); library(viridis)})

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

set.seed(315)
base <- data_dir("the Exp333 analysis outputs")
setwd(base); fdir <- "OMIQ analysis/From OMIQ"
caldir <- "R analysis/cal"; figdir <- "R analysis/figures"

ch <- readRDS(file.path(caldir,"channels.rds")); lg <- readRDS(file.path(caldir,"logicle.rds"))
th <- readRDS(file.path(caldir,"thresholds.rds")); meta <- readRDS(file.path(caldir,"samplesheet.rds"))
fr0 <- read.FCS(file.path(fdir, meta$Filename[1]), transformation=FALSE, truncate_max_range=FALSE)
spm <- spillover(fr0)[["$SPILLOVER"]]

groups <- c("wt T cells","cxcr3 ko T cell",
            "wt T cells + vaccines","cxcr3 ko T cells + vaccines")
glab <- c("wt T cells"="Cxcr3+/+ T cells","cxcr3 ko T cell"="Cxcr3-/- T cells",
          "wt T cells + vaccines"="Cxcr3+/+ T cells + vaccines",
          "cxcr3 ko T cells + vaccines"="Cxcr3-/- T cells + vaccines")
geno <- c("wt T cells"="Cxcr3+/+","cxcr3 ko T cell"="Cxcr3-/-",
          "wt T cells + vaccines"="Cxcr3+/+","cxcr3 ko T cells + vaccines"="Cxcr3-/-")

embed_markers <- c("CD44","CD62L","CD127","KLRG1","IFNg","TNF","IL2","CD107a")
keep_markers  <- c(embed_markers, "CD45.1","CD45.2")
NPER <- 7000

## gate CD3+CD8+ of one OVA_stim file; return transformed marker matrix + donor flag
get_cd8 <- function(file){
  fr <- read.FCS(file.path(fdir,file), transformation=FALSE, truncate_max_range=FALSE)
  fr <- compensate(fr, spm); e <- exprs(fr); n <- nrow(e)
  sc <- e[,c("FSC-A","SSC-A")]; ssub <- sc[sample(n,min(n,20000)),]
  rc <- tryCatch(MASS::cov.rob(ssub,method="mcd",quantile.used=floor(.6*nrow(ssub))),
                 error=function(er) list(center=colMeans(ssub),cov=cov(ssub)))
  m_lymph <- mahalanobis(sc, rc$center, rc$cov) < qchisq(.95,2)
  r <- e[,"FSC-H"]/pmax(e[,"FSC-A"],1); m_sing <- r>(median(r)-3*mad(r)) & r<(median(r)+3*mad(r))
  T <- exprs(transform(fr, lg))
  keep <- m_lymph & m_sing & T[,ch["LD"]]<th$LD & T[,ch["CD3"]]>th$CD3 & T[,ch["CD8"]]>th$CD8
  M <- T[keep, ch[keep_markers], drop=FALSE]; colnames(M) <- keep_markers
  data.frame(M, donor = T[keep, ch["CD45.1"]] > th$CD45.1)
}

## collect cells
dat <- list()
for (g in groups) {
  fs <- meta$Filename[meta$Group==g & meta$Stim_condition=="OVA_stim"]
  M  <- bind_rows(lapply(fs, get_cd8))
  if (nrow(M) > NPER) M <- M[sample(nrow(M), NPER), ]
  M$group <- factor(glab[g], levels=unname(glab)); M$geno <- geno[g]
  dat[[g]] <- M
  cat(glab[g], "cells:", nrow(M), " donor:", sum(M$donor), "\n")
}
D <- bind_rows(dat)
D$Population <- factor(ifelse(D$donor, "Transferred (CD45.1+)", "Host (CD45.1-)"),
                       levels=c("Host (CD45.1-)","Transferred (CD45.1+)"))

## ---- UMAP on scaled embedding markers ----
X <- scale(as.matrix(D[, embed_markers]))
um <- umap(X, n_neighbors=15, min_dist=0.25, metric="euclidean", n_threads=4, verbose=TRUE)
D$UMAP1 <- um[,1]; D$UMAP2 <- um[,2]
saveRDS(D, file.path(caldir,"umap_data.rds"))

th_um <- theme_classic(base_size=11) +
  theme(axis.text=element_blank(), axis.ticks=element_blank(),
        strip.background=element_blank(), strip.text=element_text(face="bold"))

## ---- Fig8a: donor vs host, faceted by group ----
ord <- order(D$donor)  # draw transferred cells on top
pa <- ggplot(D[ord,], aes(UMAP1, UMAP2, colour=Population)) +
  geom_point(size=.35, alpha=.7) +
  scale_colour_manual(values=c("Host (CD45.1-)"="grey80","Transferred (CD45.1+)"="#d7301f")) +
  facet_wrap(~group, nrow=1) +
  guides(colour=guide_legend(override.aes=list(size=3, alpha=1))) +
  labs(title="CD8+ T cells (SIINFEKL stim.)",
       x="UMAP1", y="UMAP2") + th_um
ggsave(file.path(figdir,"Fig8a_UMAP_donor_by_group.png"), pa, width=14, height=4.2, dpi=300, bg="white")
ggsave(file.path(figdir,"Fig8a_UMAP_donor_by_group.pdf"), pa, width=14, height=4.2, bg="white")

## ---- Fig8b: per-marker expression (pooled embedding) ----
ML <- D %>% select(UMAP1, UMAP2, all_of(keep_markers)) %>%
  pivot_longer(all_of(keep_markers), names_to="marker", values_to="expr") %>%
  group_by(marker) %>% mutate(expr = pmin(pmax(expr, quantile(expr,.01)), quantile(expr,.99))) %>%
  ungroup()
ML$marker <- factor(ML$marker, levels=keep_markers)
pb <- ggplot(ML, aes(UMAP1, UMAP2, colour=expr)) +
  geom_point(size=.2) + scale_colour_viridis_c(option="turbo", name="expr\n(logicle)") +
  facet_wrap(~marker, nrow=2) +
  labs(title="CD8+ T cells (SIINFEKL stim.)",
       x="UMAP1", y="UMAP2") + th_um
ggsave(file.path(figdir,"Fig8b_UMAP_markers.png"), pb, width=13, height=6, dpi=300, bg="white")
ggsave(file.path(figdir,"Fig8b_UMAP_markers.pdf"), pb, width=13, height=6, bg="white")

## ---- Fig8c: genotype overlap (donor cells only) ----
Dd <- D[D$donor, ]
pc <- ggplot(Dd, aes(UMAP1, UMAP2, colour=geno)) +
  geom_point(size=.4, alpha=.7) +
  scale_colour_manual(values=c("Cxcr3+/+"="#4292c6","Cxcr3-/-"="#e6550d"), name="Genotype") +
  guides(colour=guide_legend(override.aes=list(size=3, alpha=1))) +
  labs(title="Transferred (CD45.1+) cells only: Cxcr3+/+ vs Cxcr3-/- occupy the same UMAP space",
       x="UMAP1", y="UMAP2") + th_um
ggsave(file.path(figdir,"Fig8c_UMAP_genotype.png"), pc, width=6.5, height=5, dpi=300, bg="white")
ggsave(file.path(figdir,"Fig8c_UMAP_genotype.pdf"), pc, width=6.5, height=5, bg="white")

cat("\nSaved Fig8a/b/c UMAP figures.\n")
