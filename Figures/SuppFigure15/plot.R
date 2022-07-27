#!/usr/bin/env Rscript

library(qvalue)
setwd("/Users/nikolaoslykoskoufis/Documents/PROJECTS/te_project/paper/rerun/Figures/SuppFigure15")

pdf("TE_gene_associations_effectSize_distribution_ALL.pdf",8,4.5)
par(mfrow=c(1,2))
df <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/emods/nominal/normal/NORMAL.TE.GENE.nominal.chrAll.txt.gz", header=F,stringsAsFactors=F)
hist(df$V14, breaks=100,xlab="Effect size (regression slope)",ylab="Frequency",main="A. TE-Gene association\neffect size in normal")

df <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/emods/nominal/tumor/TUMOR.TE.GENE.nominal.chrAll.txt.gz", header=F,stringsAsFactors=F)
hist(df$V14, breaks=100,xlab="Effect size (regression slope)",ylab="Frequency",main="B. TE-Gene association\neffect size in tumor")
dev.off()
