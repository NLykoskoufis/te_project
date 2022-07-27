#!/usr/bin/env Rscript 

args = commandArgs(trailingOnly=TRUE)

library(dplyr)

file <- args[1]
out <- args[2]

pc <- read.table(file,header=T,stringsAsFactors=F)
rownames(pc) <- pc$SampleID 
samples <- colnames(pc)

gen <- read.table("/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/2backup/hongen/SYSCOL/qtltools/PCA/syscol_normal_gene_geno_merged.60pcs",header=T,stringsAsFactors=F)
gen <- gen[c(1:3),]
rownames(gen) <- gen$SampleID
gen <- gen %>% select(samples)

if(identical(colnames(gen), colnames(pc))){
      pcs <- rbind(gen,pc)
  colnames(pcs) <- gsub("X","",colnames(pcs))
    write.table(pcs, file=paste0(out,sep=""),row.names=F,col.names=T,quote=F,sep="\t")
}else{
      print("ERROR")
}

