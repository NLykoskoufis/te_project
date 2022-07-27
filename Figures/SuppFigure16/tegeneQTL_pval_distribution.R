#!/usr/bin/env Rscript 
library(qvalue)

setwd("/Users/nikolaoslykoskoufis/Documents/PROJECTS/te_project/paper/rerun/Figures/SuppFigure16")

pdf("TEgeneQTLs_association_pval_distribution.pdf",8,4.5)
par(mfrow=c(1,2))

df <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/permute/normal/NORMAL_teqtls.chrALL.txt", header=T,stringsAsFactors=F,sep=" ")
QP = qvalue(df$adj_beta_pval)


hist(df$adj_beta_pval, breaks=20, xlab="Adjusted P-values", ylab="Frequency", main="A. SNP-TE-gene associations\nin normal")
legend("topright", legend=c(paste("pi1=", signif(100-QP$pi0*100, 3), "%", sep=""), paste("#hits=", sum(QP$qvalue <= 0.05, na.rm=TRUE), " (5% FDR)", sep="")), bty="n")

df <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/permute/tumor/TUMOR_teqtls.chrALL.txt", header=T,stringsAsFactors=F,sep=" ")
QP = qvalue(df$adj_beta_pval)

hist(df$adj_beta_pval, breaks=20, xlab="Adjusted P-values", ylab="Frequency", main="B. SNP-TE-gene associations\nin tumor")
legend("topright", legend=c(paste("pi1=", signif(100-QP$pi0*100, 3), "%", sep=""), paste("#hits=", sum(QP$qvalue <= 0.05, na.rm=TRUE), " (5% FDR)", sep="")), bty="n")
dev.off()
