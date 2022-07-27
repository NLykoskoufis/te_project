#!/usr/bin/env Rscript 
setwd("~/Documents/PROJECTS/te_project/paper/rerun/Figures/SuppFigure21")

library(ggplot2)
library(ggpubr)
COL = RColorBrewer::brewer.pal(9,"Set3")



BN <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/tumor/TUMOR.BNlearn.ALL.extraINFO.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
BN.NORMAL <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/normal/NORMAL.BNlearn.ALL.extraINFO.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))

UNION <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/union/UNION_normal_tumor_with_posteriors.txt.gz",header=T,stringsAsFactors=F,sep="\t"))
UNION$triplet <- paste(UNION$var_id, UNION$gene, UNION$te, sep=";")
UNION$switch <- "switch"
UNION[which(UNION$normal_best == UNION$tumor_best),]$switch <- "no_switch"
UNION[which(UNION$normal_best != UNION$tumor_best & UNION$tumor_best == "m1"),]$switch <- "switch_causal"


TUMOR <- UNION[UNION$triplet %in% BN$triplet,]
NORMAL <- UNION[UNION$triplet %in% BN.NORMAL$triplet,]
CAUSAL <- TUMOR[which(TUMOR$switch == "switch_causal"),]


############################################################################
ANALYSIS_ = "Differential expression of TEs for Switch CAUSAL vs Others"   #
############################################################################
TEs <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/analysis_DGE/DGE_results_DESEq2_all_TE_genes.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
TEs <- TEs[which(TEs$padj <= 0.05),]
TES <- TEs[!grepl("ENSG", TEs$id),]




OTHER <- TUMOR[which(TUMOR$switch != "switch_causal"),]
OTHER_texp <- merge(OTHER, TEs, by.x="te", by.y="id")
CAUSAL_texp <- merge(CAUSAL, TEs, by.x="te", by.y="id")
common <- intersect(CAUSAL_texp$te, OTHER_texp$te)


OTHER_texp <- OTHER_texp[!OTHER_texp$te %in% common,]
CAUSAL_texp <- CAUSAL_texp[!CAUSAL_texp$te %in% common,]
OTHER_texp$type <- "Other"
CAUSAL_texp$type <- "Switch to causal"
TUMOR_texp <- rbind(CAUSAL_texp, OTHER_texp)
median(TUMOR_texp[which(TUMOR_texp$type == "Switch to causal"),]$log2FoldChange)
median(TUMOR_texp[which(TUMOR_texp$type == "Other"),]$log2FoldChange)

ggmedian <- ggplot(data=TUMOR_texp, aes(x= type, y=log2FoldChange))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test",label.x=1.4,label.y=8)+
  #annotate("text", x=1.4, y=6, label=paste("Switch to causal\nTEs median log2FC =", median(TUMOR_texp[which(TUMOR_texp$type == "Switch to causal"),]$ratio),"\nOther TEs\nmedian log2FC =",median(TUMOR_texp[which(TUMOR_texp$type == "Other"),]$ratio), sep=""))+
  labs(x="", y="log2(tumor expression / normal expression)", title="TE expression log2 fold change difference\nbetween triplets switching to causal and other triplets")+
  theme_bw()+
  theme(panel.grid=element_blank(), 
        plot.title = element_text(hjust=0.5, face="bold"))
ggsave(ggmedian, filename="te_median_expression_difference_switch_causal_others.pdf",height=5,width=5.5,device="pdf")






