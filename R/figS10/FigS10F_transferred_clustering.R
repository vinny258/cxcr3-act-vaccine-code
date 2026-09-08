## =====================================================================
## Exp 315 - STEP 17: clustering WITHIN the transferred (CD45.1+) cells
## Linear-story analysis: having shown comparable expansion (Fig 2) and IFNg
## (Fig 4/5), and the differentiation landscape (Fig 8a/b), we now FOCUS on the
## transferred (CD45.1+) OVA-specific cells, cluster THEM (not all CD8), and ask
## whether the subpopulations differ between WT+Vax and KO+Vax.
##  - Cells: OVA-stim tubes, CD45.1+ only, WT+Vax & KO+Vax (non-MVA)
##  - Markers: CD44, CD62L, CD127, KLRG1, IFNg, TNF, IL2, CD107a (as Fig 8a/b)
##  - k chosen objectively; per-mouse cluster frequencies; WT vs KO (MWU + BH)
## Outputs: figures/Fig10_*.{png,pdf,_AI.pdf} ; results/transferred_cluster_*.csv
## =====================================================================
suppressMessages({library(flowCore); library(dplyr); library(tidyr); library(ggplot2)
                  library(cluster); library(uwot); library(pheatmap); library(ggpubr); library(ggrastr)})
set.seed(315)
base <- "C:/Users/valmeida/Ludwig Institute for Cancer Research Dropbox/Vinnycius Pereira Almeida/BVDE Lab/Mouse Exp/B6/Exp 315"
setwd(base); fdir <- "OMIQ analysis/From OMIQ"
caldir <- "R analysis/cal"; resdir <- "R analysis/results"; figdir <- "R analysis/figures"

markers <- c("CD44","CD62L","CD127","KLRG1","IFNg","TNF","IL2","CD107a")
MIN_TRANSF <- 50
K <- 6   ## set after viewing Fig10_chooseK (gap tends to over-split continuous data)

## ---- self-contained extraction (needs only the cal/ calibration from step 03) ----
ch <- readRDS(file.path(caldir,"channels.rds")); lg <- readRDS(file.path(caldir,"logicle.rds"))
th <- readRDS(file.path(caldir,"thresholds.rds")); meta <- readRDS(file.path(caldir,"samplesheet.rds"))
fr0 <- read.FCS(file.path(fdir, meta$Filename[1]), transformation=FALSE, truncate_max_range=FALSE)
spm <- spillover(fr0)[["$SPILLOVER"]]
grp_geno <- c("wt T cells + vaccines"="Cxcr3+/+","cxcr3 ko T cells + vaccines"="Cxcr3-/-")

get_transferred <- function(file, mouse, group){
  fr <- read.FCS(file.path(fdir,file), transformation=FALSE, truncate_max_range=FALSE)
  fr <- compensate(fr, spm); e <- exprs(fr); n <- nrow(e)
  sc <- e[,c("FSC-A","SSC-A")]; ssub <- sc[sample(n,min(n,20000)),]
  rc <- tryCatch(MASS::cov.rob(ssub,method="mcd",quantile.used=floor(.6*nrow(ssub))),
                 error=function(er) list(center=colMeans(ssub),cov=cov(ssub)))
  m_lymph <- mahalanobis(sc, rc$center, rc$cov) < qchisq(.95,2)
  r <- e[,"FSC-H"]/pmax(e[,"FSC-A"],1); m_sing <- r>(median(r)-3*mad(r)) & r<(median(r)+3*mad(r))
  T <- exprs(transform(fr, lg))
  keep <- m_lymph & m_sing & T[,ch["LD"]]<th$LD & T[,ch["CD3"]]>th$CD3 &
          T[,ch["CD8"]]>th$CD8 & T[,ch["CD45.1"]]>th$CD45.1          # CD45.1+ transferred
  M <- T[keep, ch[markers], drop=FALSE]; colnames(M) <- markers
  if (nrow(M)==0) return(NULL)
  data.frame(Mouse=mouse, geno=unname(grp_geno[group]), M)
}
ss <- meta %>% filter(Group %in% names(grp_geno), Stim_condition=="OVA_stim")
D <- bind_rows(lapply(seq_len(nrow(ss)), function(i)
       get_transferred(ss$Filename[i], ss$Mouse_number[i], ss$Group[i])))
D$geno <- factor(D$geno, levels=c("Cxcr3+/+","Cxcr3-/-"))
cat("Transferred cells: WT+Vax =", sum(D$geno=="Cxcr3+/+"),
    " KO+Vax =", sum(D$geno=="Cxcr3-/-"), "\n")

## ---- re-scale on the transferred cells themselves ----
ctr <- colMeans(D[,markers]); sdv <- apply(D[,markers],2,sd)
Z <- scale(as.matrix(D[,markers]), center=ctr, scale=sdv)

## ---- choose k (elbow / silhouette / gap) ----
ks<-1:12; wss<-sapply(ks,function(k){set.seed(1);kmeans(Z[sample(nrow(Z),min(20000,nrow(Z))),],k,nstart=10,iter.max=50)$tot.withinss})
ss<-Z[sample(nrow(Z),3000),]; dss<-dist(ss); ks2<-2:10
sil<-sapply(ks2,function(k){set.seed(1);mean(silhouette(kmeans(ss,k,nstart=10)$cluster,dss)[,3])})
gp<-clusGap(Z[sample(nrow(Z),2000),],FUN=function(x,k)list(cluster=kmeans(x,k,nstart=10,iter.max=50)$cluster),K.max=10,B=25,verbose=FALSE)
gt<-as.data.frame(gp$Tab); k_gap<-maxSE(gt$gap,gt$SE.sim,"Tibs2001SEmax")
cat(sprintf("Suggested k: silhouette=%d ; gap=%d ; using K=%d\n", ks2[which.max(sil)], k_gap, K))
p1<-ggplot(data.frame(k=ks,wss=wss),aes(k,wss))+geom_line()+geom_point()+scale_x_continuous(breaks=ks)+labs(title="Elbow",y="within SS")+theme_bw(base_size=11)
p2<-ggplot(data.frame(k=ks2,sil=sil),aes(k,sil))+geom_line()+geom_point()+scale_x_continuous(breaks=ks2)+labs(title="Silhouette")+theme_bw(base_size=11)
p3<-ggplot(data.frame(k=seq_len(nrow(gt)),gap=gt$gap,se=gt$SE.sim),aes(k,gap))+geom_line()+geom_point()+geom_errorbar(aes(ymin=gap-se,ymax=gap+se),width=.2)+scale_x_continuous(breaks=1:10)+labs(title="Gap statistic")+theme_bw(base_size=11)
ggsave(file.path(figdir,"Fig10_chooseK.png"),ggarrange(p1,p2,p3,ncol=3),width=13,height=4,dpi=300,bg="white")

## ---- k-means on transferred cells ----
km <- kmeans(Z, centers=K, nstart=30, iter.max=100)
D$cluster <- factor(km$cluster, levels=as.character(sort(unique(km$cluster))))

## ---- cluster phenotype (median logicle) ----
cph <- D %>% group_by(cluster) %>% summarise(across(all_of(markers), median), n=n())
write.csv(cph, file.path(resdir,"transferred_cluster_phenotype.csv"), row.names=FALSE)
cat("\nCluster phenotype (median logicle):\n"); print(as.data.frame(cph), digits=2)

## ---- per-mouse cluster frequencies + WT vs KO (MWU, BH) ----
tf <- D %>% group_by(Mouse) %>% mutate(ntot=n()) %>% ungroup() %>% filter(ntot>=MIN_TRANSF)
freq <- tf %>% group_by(Mouse,geno,cluster,ntot) %>% summarise(nc=n(),.groups="drop") %>%
  mutate(pct=100*nc/ntot) %>% complete(nesting(Mouse,geno,ntot),cluster,fill=list(nc=0,pct=0))
write.csv(freq, file.path(resdir,"transferred_cluster_freq_permouse.csv"), row.names=FALSE)
stat <- freq %>% group_by(cluster) %>%
  summarise(n_WT=sum(geno=="Cxcr3+/+"), n_KO=sum(geno=="Cxcr3-/-"),
            WT_mean=mean(pct[geno=="Cxcr3+/+"]), KO_mean=mean(pct[geno=="Cxcr3-/-"]),
            p=tryCatch(suppressWarnings(wilcox.test(pct~geno)$p.value),error=function(e)NA),.groups="drop") %>%
  mutate(p_BH=p.adjust(p,"BH"))
write.csv(stat, file.path(resdir,"transferred_cluster_WTvsKO_stats.csv"), row.names=FALSE)
cat("\nWT+Vax vs KO+Vax per cluster:\n"); print(as.data.frame(stat), digits=3)

## =========================================================
## FIGURES
## =========================================================
## UMAP of the transferred cells
um <- umap(Z, n_neighbors=15, min_dist=0.25, n_threads=4)
D$UMAP1<-um[,1]; D$UMAP2<-um[,2]
th_um<-theme_classic(base_size=11)+theme(axis.text=element_blank(),axis.ticks=element_blank(),strip.background=element_blank(),strip.text=element_text(face="bold"))
pal_k<-setNames(scales::hue_pal()(K), levels(D$cluster))
## draw smaller clusters on top
sizes<-table(D$cluster); D<-D[order(-as.integer(sizes[as.character(D$cluster)])),]
cent<-D %>% group_by(cluster) %>% summarise(UMAP1=median(UMAP1),UMAP2=median(UMAP2))

## Fig10a: clusters, WT vs KO side by side
p10a<-ggplot(D,aes(UMAP1,UMAP2,colour=cluster))+geom_point(size=.45,alpha=.75)+
  scale_colour_manual(values=pal_k,name="Cluster",drop=TRUE)+
  geom_label(data=cent,aes(x=UMAP1,y=UMAP2,label=paste0("C",cluster)),inherit.aes=FALSE,colour="black",fill="white",alpha=.7,size=2.8,label.size=0)+
  facet_wrap(~geno, labeller=as_labeller(c("Cxcr3+/+"="Cxcr3+/+ T cells + vaccines",
                                           "Cxcr3-/-"="Cxcr3-/- T cells + vaccines")))+
  guides(colour=guide_legend(override.aes=list(size=3,alpha=1)))+
  labs(title="CD45.1+ T cells (SIINFEKL stim.)")+th_um
ggsave(file.path(figdir,"Fig10a_transferred_clusters.png"),p10a,width=11,height=5,dpi=300,bg="white")
ggsave(file.path(figdir,"Fig10a_transferred_clusters.pdf"),p10a,width=11,height=5,bg="white")
ggsave(file.path(figdir,"Fig10a_AI.pdf"),rasterise(p10a,layers="Point",dpi=300),width=11,height=5,bg="white")

## Fig10b: cluster phenotype heatmap
hm<-as.matrix(cph[,markers]); rownames(hm)<-paste0("C",cph$cluster)
png(file.path(figdir,"Fig10b_cluster_heatmap.png"),width=1000,height=850,res=160)
pheatmap(scale(hm),cluster_cols=FALSE,main="CD45.1+ cluster phenotype (z-scored median)",display_numbers=round(hm,1),fontsize=11);dev.off()
pdf(file.path(figdir,"Fig10b_cluster_heatmap.pdf"),width=6.5,height=5.5)
pheatmap(scale(hm),cluster_cols=FALSE,main="CD45.1+ cluster phenotype (z-scored median)",display_numbers=round(hm,1),fontsize=11);dev.off()

## Fig10c: per-cluster frequency, WT+Vax vs KO+Vax, with stats
freq$geno<-factor(freq$geno,levels=c("Cxcr3+/+","Cxcr3-/-"))
p10c<-ggplot(freq,aes(geno,pct,fill=geno))+geom_boxplot(outlier.shape=NA,alpha=.7)+geom_jitter(width=.15,height=0,size=1.4)+
  facet_wrap(~cluster,scales="free_y",nrow=2,labeller=label_both)+
  scale_fill_manual(values=c("Cxcr3+/+"="#4292c6","Cxcr3-/-"="#e6550d"),guide="none")+
  scale_x_discrete(labels=c("Cxcr3+/+"="Cxcr3+/+\n+ vaccines","Cxcr3-/-"="Cxcr3-/-\n+ vaccines"))+
  stat_compare_means(method="wilcox.test",label="p.format",size=3,label.x.npc="center")+
  labs(x=NULL,y="% of transferred cells in cluster")+
  theme_bw(base_size=10)+theme(panel.grid.minor=element_blank())
ggsave(file.path(figdir,"Fig10c_cluster_freq_WTvsKO.png"),p10c,width=11,height=6,dpi=300,bg="white")
ggsave(file.path(figdir,"Fig10c_cluster_freq_WTvsKO.pdf"),p10c,width=11,height=6,bg="white")

cat("\nSaved Fig10a/b/c and transferred_cluster_*.csv\n")
