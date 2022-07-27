#!/usr/bin/env Rscript 


library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggsci)
library(qvalue)

tumor <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/tumor/SYSCOL_tumor_276_gene_TE.FM05.resid.noRNT.inR.cpm.bed.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
MSS.samples <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/revisions/info/MSS.samples",header=FALSE,stringsAsFactors=FALSE,sep="\n")$V1

head(tumor)
rownames(tumor) <- tumor$gene
tumor <- tumor[,-c(1:6)]

TMP <- tumor[!grepl("ENSG",rownames(tumor)),]

SS <- data.frame(samples = names(tumor))
SS$status <- "MSI" 
SS$status[SS$samples %in% MSS.samples] <- "MSS"

dd <- data.frame()
pdf("mss_msi_tePlots_Log2Exp.pdf",8,8)
par(mfrow=c(2,2))
for (rrr in 1:nrow(TMP))
{
  cat("Processed",rrr,"-",nrow(TMP), "\n")
  te <- rownames(TMP[rrr,])
  x <- as.data.frame(t(TMP[rrr,]),stringsAsFactors = FALSE)
  colnames(x) <- c("exp")
  x$samples <- rownames(x)
  
  x <- merge(x,SS, by="samples")
  
  MSS <- x[which(x$status == "MSS"),]$exp  
  MSI <- x[which(x$status == "MSI"),]$exp  
  tt <- wilcox.test(MSS,MSI)
  tt.oneSide <- wilcox.test(MSI,MSS, alternative = "greater")
  if (as.numeric(tt$p.value) < 0.05){
    boxplot(log2(x$exp) ~ x$status, ylab="log2 normalized TE expression [log2 CPM]", main=paste0(te,"\np-val=",signif(as.numeric(tt$p.value)),3))
  }
    dd <- rbind(dd, data.frame("id"=te,"mss_median_exp"=median(MSS),"msi_median_exp"=median(MSI), "pvalue" = as.numeric(tt$p.value),"pvalue_one_side"=as.numeric(tt.oneSide$p.value)))
}
dev.off()




dd$qvalue <- qvalue(dd$pvalue)$qvalues
dd$qvalue_one_side <- qvalue(dd$pvalue_one_side)$qvalues 
dd[which(dd$qvalue <= 0.05),]
write.table(dd, file="mss_msi_wilcox_analysis.txt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

#dd <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/revisions/analysis_MSS_MSI_TE_expression/mss_msi_wilcox_analysis.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))


pdf("pvalueDistribution.pdf",4,4.5)
hist(dd$pvalue, main="MSI vs MSS samples\nTE expression",xlab="p-value distribution")
hist(dd$qvalue, main="MSI vs MSS samples\nTE expression",xlab="q-value distribution")
dev.off()


pdf("mss_msi_tePlots_Log2Exp_SignificantOnesFDR005.pdf",height=4.5,width=13)
par(mfrow=c(1,3))
rownames(dd) <- dd$id 
FDR_IDs <- dd[which(dd$qvalue <= 0.05),]$id
for (rrr in FDR_IDs)
{
  #cat("Processed",rrr,"-",nrow(TMP), "\n")
  #te <- rownames(TMP[rrr,])
  te <- rrr
  x <- as.data.frame(t(TMP[rrr,]),stringsAsFactors = FALSE)
  colnames(x) <- c("exp")
  x$samples <- rownames(x)
  
  x <- merge(x,SS, by="samples")
  
  boxplot(log2(x$exp) ~ x$status, ylab="log2 normalized TE expression [log2 CPM]", main=paste0(te,"\nq-value=",format.pval(as.numeric(dd[rrr,]$qvalue))), xlab="patient status")
  
}
dev.off()
