#!/usr/bin/env Rscript 

setwd("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn")
BN <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/tumor/TUMOR.BNlearn.ALL.extraINFO.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
#BN$triplet <- paste(BN$var_id, BN$gene, BN$te, sep="_")
BN.NORMAL <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/normal/NORMAL.BNlearn.ALL.extraINFO.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
#BN.NORMAL$triplet <- paste(BN.NORMAL$var_id, BN.NORMAL$gene, BN.NORMAL$te, sep="_")

UNION <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/union/UNION_normal_tumor_with_posteriors.txt.gz",header=T,stringsAsFactors=F,sep="\t"))
UNION$triplet <- paste(UNION$var_id, UNION$gene, UNION$te, sep=";")
UNION$switch <- "switch"
UNION[which(UNION$normal_best == UNION$tumor_best),]$switch <- "no_switch"
UNION[which(UNION$normal_best != UNION$tumor_best & UNION$tumor_best == "m1"),]$switch <- "switch_causal"


TUMOR <- UNION[UNION$triplet %in% BN$triplet,]
NORMAL <- UNION[UNION$triplet %in% BN.NORMAL$triplet,]
CAUSAL <- TUMOR[which(TUMOR$switch == "switch_causal"),]

ensembl<- read.table("/Users/srv/beegfs/scratch/groups/data/nikos/annotations/bed/ensemblID_to_genename.txt.gz",header=FALSE,stringsAsFactors=F)
ensembl$V1 <- sapply(strsplit(ensembl$V1,"\\."),"[[",1)

TCGT <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/analysis_TcTs/TCGT_PROCESSED_ALL_withFisher.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))

findTCGT <- function(tcgts, file){
  ndata <- data.frame()
  for(i in 1:nrow(tcgts))
  {
    r <- tcgts[i,]
    bb <- file[which(file$te == r$id),]
    b <- bb[grepl(r$tcgt_gene, bb$gene),]
    if(nrow(b) != 0){
      ndata <- rbind(ndata,cbind(b,r))
    }
  }
  return(ndata)
}

CAUSAL.TCGT <- findTCGT(TCGT[which(TCGT$fisher_pvalue < 0.05),],CAUSAL)
TUMOR.TCGT <- findTCGT(TCGT[which(TCGT$fisher_pvalue < 0.05),],TUMOR)
# 46 triplets 
# 14 unique genes  
# 39 unique TEs
ensembl[ensembl$V1 %in% CAUSAL.TCGT$gene ,]
