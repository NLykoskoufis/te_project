#!/usr/bin/env Rscript 

setwd("/Users/nikolaoslykoskoufis/Documents/PROJECTS/te_project/paper/rerun/Figures/SuppFigure11/")
##LIBRARIES##
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)


filterEnrichment <- function(t)
{
  
  t.sig.both <- t[which(t$sig.x == TRUE & t$sig.y == TRUE),]
  t.sig.x.not.y <- t[which(t$sig.x == TRUE & t$sig.y == FALSE),]
  t.sig.y.not.x <- t[which(t$sig.y == TRUE & t$sig.x == FALSE),]
  
  pb1 <- t.sig.y.not.x[t.sig.y.not.x$OddsRatio.x > t.sig.y.not.x$OddsRatio.y,]
  pb2 <- t.sig.x.not.y[t.sig.x.not.y$OddsRatio.y > t.sig.x.not.y$OddsRatio.x,]
  print(c(pb1$mark,pb2$mark))
  t.sig.y.not.x <- t.sig.y.not.x[!t.sig.y.not.x$OddsRatio.x > t.sig.y.not.x$OddsRatio.y,]
  t.sig.x.not.y <- t.sig.x.not.y[!t.sig.x.not.y$OddsRatio.y > t.sig.x.not.y$OddsRatio.x,]
  toKeep <- unique(c(t.sig.both$mark, t.sig.x.not.y$mark, t.sig.y.not.x$mark))
  return(toKeep)
}



shared_te <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/specific/shared/shared_TEeQTL_enrichment_results_maf002_w2500_k100.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
tspe_te <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/specific/tumor_specific/tumor_specific_TEeQTL_enrichment_results_maf002_w2500_k100.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))


tspe_te$sig <- tspe_te$fdr <= 0.05
tspe_te$sig10 <- tspe_te$fdr <= 0.1
shared_te$sig <- shared_te$fdr <= 0.05
shared_te$sig10 <- shared_te$fdr <= 0.1



# combined tumor-spe and shared plot
comb_te <- merge(tspe_te, shared_te, by="mark")
t <- comb_te 

toKeep_te <- filterEnrichment(t)
toKeep <- t[which(t$sig.x | t$sig.y),]$mark




# if tumor specific is not significant but has higher OR than normal and exclude. If shared is not significant and has higher OR than tumor-specific, exclude. 

TMP_tspe_te <- tspe_te[tspe_te$mark %in% toKeep_te,c(1,7,8)]
TMP_tspe_te$type <- "tumor-specific"

TMP_shared_te <- shared_te[shared_te$mark %in% toKeep_te, c(1,7,8)]
TMP_shared_te$type <- "shared"

TMP_te <- rbind(TMP_shared_te,TMP_tspe_te)
TMP_te$OddsRatio[which(TMP_te$OddsRatio > 30)] <- 30
TMP_te$ypos <- 0
for(i in unique(TMP_te$mark)){
  x <- TMP_te[which(TMP_te$mark == i),]
  mm <- max(x$OddsRatio)
  TMP_te[which(TMP_te$mark == i & TMP_te$type == "tumor-specific"),]$ypos <- mm + 1.5
  TMP_te[which(TMP_te$mark == i & TMP_te$type == "shared"),]$ypos <- mm + 4
}


subset_to_order <- subset(TMP_te,type == "tumor-specific")
subset_to_order$mark <- with(subset_to_order, reorder(mark,-OddsRatio))
TMP_te$mark = factor(TMP_te$mark, levels = levels(subset_to_order$mark))


gg_combined_te <- ggplot(data=TMP_te, aes(x=mark, y=OddsRatio, fill=type,label=round(-log10(pvalue),1)))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_hline(yintercept=1, lty=2, colour="grey")+
  geom_text(aes(x=mark, y=ypos),angle=90,size=3)+
  theme_classic(base_size=20)+
  labs(x="",y="Enrichment over the null",title="Functional enrichment for tumor-specific and shared TE-eQTLs")+
  scale_fill_manual(values=c("tumor-specific"="#e41a1c", "shared"="black"))+
  theme(axis.text.x = element_text(size=9,angle=90,hjust=1,vjust=0.5),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.9,0.9), legend.title=element_blank())
ggsave(gg_combined_te, filename="functional_enrichment_combined_tspe_shared_TE_eQTLs_k100_maf002_w2500.pdf",device="pdf",height=7,width=20)
