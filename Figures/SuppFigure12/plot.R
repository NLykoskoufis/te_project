#!/usr/bin/env Rscript 

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


### DGE FILES ### 
tf <- as.data.frame(data.table::fread("~/Documents/PROJECTS/te_project/paper/Files/forPlotting/SuppFigure10/DGE_results_DESEq2_all_TE_genes.txt",header=T,stringsAsFactors = F,sep="\t"))
convert <- read.table("~/Documents/PROJECTS/te_project/paper/Files/forPlotting/SuppFigure10ensemblToGeneName.txt",header=F,sep=" ",stringsAsFactors = F)
colnames(convert) <- c("ensemblID","geneName")
tf <- merge(convert,tf,by.x="ensemblID",by.y="id")

#### changing names of tf genes to match ddte. 
tf[which(tf$ensemblID == "ENSG00000085276"),]$geneName <- "EVI1"
tf[which(tf$ensemblID == "ENSG00000125820"),]$geneName <- "NKX2.2"
tf[which(tf$ensemblID == "ENSG00000157557"),]$geneName <- "Ets-2"
tf <- tf[tf$sig,]
### ADD GENE FOLD CHANGE ###

comb_te$geneName <- sapply(strsplit(comb_te$mark, "\\ "),"[[",2)
comb_te <- comb_te[comb_te$mark %in% toKeep_te,]
comb_te$ratio <- log2(comb_te$OddsRatio.x / comb_te$OddsRatio.y)

df <- merge(comb_te, tf, by="geneName")

ggscat_te <- ggscatter(df[df$mark %in% toKeep_te,], x = "ratio", y = "log2FoldChange", 
          add = "reg.line", conf.int = FALSE, 
          cor.coef = TRUE, cor.method = "pearson",add.params = list(color="red"),
          xlab = "log2(tumor-specific / shared TE-eQTL enrichment)", ylab = "log2(tumor / normal TE expression)")
ggsave(ggscat_te, filename="tf_exp_vs_enrichment_tspe_shared_TE_eQTL.pdf",height=5,width=5.5,device="pdf")


