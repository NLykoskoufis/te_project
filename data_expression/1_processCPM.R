#!/usr/bin/env Rscript 

setwd("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression")

if (! file.exists(c("SYSCOL_normal_275_gene_TE.raw.cpm.percentage.txt","SYSCOL_tumor_276_gene_TE.raw.cpm.percentage.txt")))
{
  te <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/raw/SYSCOL_allSamples_TEs.raw.cpm.bed.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
  gene <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/raw/SYSCOL_allSamples_genes.raw.cpm.bed.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
    
  # I need to change the names of the samples that have R to match with genotypes ==> Done in previous step when computing CPMs but does not hurt to redo.
  names(te) <- gsub("R","T", names(te))
  names(gene) <- gsub("R","T", names(gene))
  
  n.samples <- as.character(read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/stat/syscol_normal_samples.txt",header=F,sep="\t",stringsAsFactors = FALSE))
  t.samples <- as.character(read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/stat/syscol_tumor_samples.txt",header=F,sep="\t",stringsAsFactors = FALSE))
  
  
  ereinfo <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/info_files/ereinfoTable.txt.gz", header=TRUE,stringsAsFactors = FALSE,sep="\t"))
  ereinfo$teID <- paste(ereinfo$unMerged, ereinfo$start, ereinfo$end, ereinfo$stran, sep="_")
  
  nrow(te)
  te <- merge(te, ereinfo[,c(1,9)], by.x="id",by.y="repmaskID")
  te <- dplyr::select(te, c("id","#chr", "start", "end", "teID", "info","strand", names(te)[7:590]))
  te$info <- paste0(te$info, ";RM=",te$id)
  te <- te[,-1]
  names(te)[4] <- "id"
  
  TMP <- rbind(gene,te)
  
  TMP$"#chr" <- gsub("^chr","", TMP$"#chr")
  TMP <- TMP[! TMP$"#chr" %in% c("M","X","Y"),]
  TMP$"#chr" <- as.numeric(as.character(TMP$"#chr"))
  TMP <- TMP[grepl("protein_coding|lincRNA|transposable_element",TMP$info),]
  
  TMP <- TMP[with(TMP, order(`#chr`,start, end)),]
  data.table::fwrite(TMP, file = "raw/SYSCOL_ALLsamples_normal_tumor_gene_TE.raw.cpm.bed", sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
  
  TMP.n <- dplyr::select(TMP, c("#chr","start","end","id","info", "strand",n.samples))
  TMP.t <- dplyr::select(TMP, c("#chr","start","end","id","info", "strand",t.samples))
  
  data.table::fwrite(TMP.n, file = "raw/SYSCOL_normal_275_gene_TE.raw.cpm.bed", sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)  
  data.table::fwrite(TMP.t, file = "raw/SYSCOL_tumor_276_gene_TE.raw.cpm.bed", sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}else
{
  TMP.n <- as.data.frame(data.table::fread("raw/SYSCOL_normal_275_gene_TE.raw.cpm.bed",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
  TMP.t <- as.data.frame(data.table::fread("raw/SYSCOL_tumor_276_gene_TE.raw.cpm.bed",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
  # load percentages #
  n.perc <- as.data.frame(data.table::fread("raw/SYSCOL_normal_275_gene_TE.raw.cpm.percentage.txt", header=TRUE,stringsAsFactors = FALSE,sep="\t"))
  t.perc <- as.data.frame(data.table::fread("raw/SYSCOL_tumor_276_gene_TE.raw.cpm.percentage.txt", header=TRUE,stringsAsFactors = FALSE,sep="\t"))
}

n.keep <- n.perc[which(n.perc$perc_missing < 0.5),]$gene
t.keep <- t.perc[which(t.perc$perc_missing < 0.5),]$gene

toKeep <- union(n.keep, t.keep)


TMP.n <- TMP.n[TMP.n$id %in% toKeep ,]
TMP.t <- TMP.t[TMP.t$id %in% toKeep ,]

data.table::fwrite(TMP.n, file = "normal/SYSCOL_normal_275_gene_TE.FM05.raw.cpm.bed", sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)  
data.table::fwrite(TMP.t, file = "tumor/SYSCOL_tumor_276_gene_TE.FM05.raw.cpm.bed", sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)


##### COMPARING WITH PREVIOUS ###### 

old <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_expression/normal/raw/CPM_trimmed_Genes_proteinCoding_lincRNA_TEs_ALL_275_normal_F05normal.raw.bed.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))

library(ggvenn)
old.gene <- old[grepl("ENSG",old$gene),]$gene 
old.te <- old[!grepl("ENSG",old$gene),]$gene 
new.gene <- TMP.n[grepl("ENSG",TMP.n$id),]$id 
new.te <- TMP.n[!grepl("ENSG",TMP.n$id),]$id 

ggGene <- ggvenn(data=list("old genes"=old.gene, "new genes"=new.gene))
ggTE <- ggvenn(data=list("old TEs"=old.te, "new TEs"=new.te))
ggsave(ggGene, filename="stat/venn_old_vs_new_genes_afterFiltering.pdf",device="pdf", height=6,width=6.5)
ggsave(ggTE, filename="stat/venn_old_vs_new_TEs_afterFiltering.pdf",device="pdf", height=6,width=6.5)
