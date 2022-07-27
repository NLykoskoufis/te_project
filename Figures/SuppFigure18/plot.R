#!/usr/bin/env Rscript 


suppressMessages(library(data.table))
suppressMessages(library(foreach))
suppressMessages(library(doMC))
suppressMessages(library(lmerTest))
suppressMessages(library(GenABEL))
suppressMessages(source("/home/users/l/lykoskou/bin/tools.R"))
suppressMessages(library(ggplot2))
suppressMessages(library(qvalue))
library(dplyr)
library(RColorBrewer)
COL = brewer.pal(9,"Set3")
cpus <- 20
registerDoMC(cpus)

results_normal <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/normal/NORMAL.BNlearn.ALL.extraINFO.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
results_tumor <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/tumor/TUMOR.BNlearn.ALL.extraINFO.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))

results_tumor$max <- pmax(results_tumor$M1, results_tumor$M2,results_tumor$M3)
results_normal$max <- pmax(results_normal$M1, results_normal$M2,results_normal$M3)



pdf("probability_of_most_likely_models_normal_tumor.pdf",8,4.5)
par(mfrow=c(1,2))
hist(results_normal$max, xlab="Probability of the most likely model", breaks=40, main="A. Normal",ylim=c(0,8000))
hist(results_tumor$max, xlab="Probability of the most likely model", breaks=40, main="B. Tumor",ylim=c(0,8000))
dev.off()



