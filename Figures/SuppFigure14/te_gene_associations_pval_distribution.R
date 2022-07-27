#!/usr/bin/env Rscript

setwd("/Users/nikolaoslykoskoufis/Documents/PROJECTS/te_project/paper/rerun/Figures/SuppFigure14")

library(qvalue)

#### NUMBER OF GENES ASSOCIATED PER TE #####s

pdf("N_genes_per_TE.pdf",8,4.5)
par(mfrow=c(1,2))
df <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/emods/nominal/normal/NORMAL.TE.GENE.nominal.chrAll.txt.gz", header=F,stringsAsFactors=F)
df <- as.data.frame(table(df$V8))
hist(df$Freq, breaks=20, xlab="Number of genes associated per TE", ylab="Frequency", main="Number of genes associated\nper TE in normal")
legend("topright", legend=c(paste("#mean number of\nassociated genes\nper TE=", round(mean(df$Freq),1), sep="")), bty="n")

df <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/emods/nominal/tumor/TUMOR.TE.GENE.nominal.chrAll.txt.gz", header=F,stringsAsFactors=F)
df <- as.data.frame(table(df$V8))
hist(df$Freq, breaks=20, xlab="Number of genes associated per TE", ylab="Frequency", main="Number of genes associated\nper TE in tumor")
legend("topright", legend=c(paste("#mean number of\nassociated genes\nper TE=", round(mean(df$Freq),1), sep="")), bty="n")

dev.off()
