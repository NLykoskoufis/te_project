#!/usr/bin/env Rscript 

setwd("/Users/nikolaoslykoskoufis/Documents/PROJECTS/te_project/paper/rerun/Figures/SuppFigure13")
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

#************************************#
#    Functional enrichment for TEs   # 
#************************************#

# Individual plots
tspe_te$sig <- tspe_te$fdr <= 0.05
tspe_te$sig10 <- tspe_te$fdr <= 0.1
shared_te$sig <- shared_te$fdr <= 0.05
shared_te$sig10 <- shared_te$fdr <= 0.1


# combined tumor-spe and shared plot
comb_te <- merge(tspe_te, shared_te, by="mark")
t <- comb_te 

toKeep_te <- filterEnrichment(t)
toKeep <- t[which(t$sig.x | t$sig.y),]$mark




# CHECKING WHETHER STRONGER ENRICHMENT TUMOR SPECIFIC QTLs ARE INACTIVE / ACTIVE FOR OTHER FEATURES ----------------------------
vtf <- as.data.frame(data.table::fread("/Users/nikolaoslykoskoufis/Documents/PROJECTS/te_project/paper/Files/forPlotting/SuppFigure11/VARIANT_LoVo_ChIPseq_overlap.bed.gz",stringsAsFactors=FALSE,sep="\t",header=FALSE))
vtf$V9 <- gsub("LoVo_|_peak.*","", vtf$V9)
QTL <- as.data.frame(data.table::fread("/Users/nikolaoslykoskoufis/Documents/PROJECTS/te_project/paper/Files/forPlotting/SuppFigure11/TUMOR.specific_TE_eQTLs.txt.gz",header=TRUE,stringsAsFactors=FALSE,sep=" "))
NQTL <- as.data.frame(data.table::fread("/Users/nikolaoslykoskoufis/Documents/PROJECTS/te_project/paper/Files/forPlotting/SuppFigure11/NORMAL.conditional_All.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep= " "))

comb_te$ratio <- log2(comb_te$OddsRatio.x / comb_te$OddsRatio.y)
tspe_strong_enrich <- gsub("LoVo ","",comb_te[which(comb_te$ratio >0 & comb_te$mark %in% toKeep_te),]$mark)

df <- vtf[vtf$V9 %in% tspe_strong_enrich,]
df <- merge(QTL,df,by.x="var_id",by.y="V4")
df$qtl <- paste(df$var_id, df$phe_id, sep=";")

nn <- data.frame()
for(i in unique(df$var_id))
{
  rrr <- NQTL[which(NQTL$var_id == i),]
  if(nrow(rrr) != 0)
  {
    if(length(unique(rrr$bwd_sig)) > 1)
    {
      nn <- rbind(nn, data.frame("var_id"=i, "bwd_sig"=1))
    }else
      {
        nn <- rbind(nn, data.frame("var_id"=i, "bwd_sig"=unique(rrr$bwd_sig)))
      }
  }else{
    nn <- rbind(nn, data.frame("var_id"=i, "bwd_sig"=0))
  }
}

pdf("functional_enrichment_tumor_spe_eQTL_stronger_eqtl_activity_in_normal.pdf",6,5.5)
bp<-barplot(table(nn$bwd_sig),main="tumor-specific TE-eQTLs with stronger enrichment\nfor TFs than shared TE eQTLs\nactivity in normal for other features",names.arg = c("Inactive\n (genome wide 5% FDR)","Active\n(genome wide 5% FDR)"))
balance <- table(nn$bwd_sig)
text(bp, balance/2, labels = round(balance, digits = 2))
dev.off()
