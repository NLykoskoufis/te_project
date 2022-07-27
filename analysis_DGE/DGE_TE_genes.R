#!/usr/bin/env Rscript 
library(DESeq2)
library(dplyr)
library(BiocParallel)
register(MulticoreParam(6))


COL <- RColorBrewer::brewer.pal(8,"Dark2")
col = list("TE"=COL[3], "GENE"=COL[1])


bed <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/raw/SYSCOL_normal_275_tumor_276_gene_TE.raw.counts.bed.gz",header=TRUE,stringsAsFactors=FALSE,sep="\t")) 
rownames(bed) <- bed$id 



## COVARIATES ## 
QC <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_covariates/master_qc_230814.txt", header=TRUE,stringsAsFactors = FALSE,sep="\t"))

QC <- QC[QC$sample %in% names(bed)[-c(1:6)],]
QC <- QC %>% dplyr::select(c("sample", "date", "GC_mean", "INSERT_SIZE_MODE"))
rownames(QC) <- QC$sample


coldata <- data.frame("sample"=names(bed)[-(1:6)])
coldata$tissue <- "normal"
coldata$tissue[!grepl("N", coldata$sample)] <- "tumor"
coldata <- merge(coldata, QC, by="sample")
rownames(coldata) <- coldata$sample
coldata <- coldata[,-1,drop=FALSE]
coldata$tissue <- as.factor(coldata$tissue)

countdata <- bed[,-c(1:6)]
countdata <- countdata %>% select(rownames(coldata))

dds <- DESeqDataSetFromMatrix(countData = countdata, 
                              colData = coldata,
                              design= ~ date + GC_mean + INSERT_SIZE_MODE + tissue)

dds <- DESeq(dds,parallel=TRUE)

norm.counts <- as.data.frame(counts(dds, normalized = TRUE))
norm.counts$id <- rownames(norm.counts)
norm.counts <- merge(bed[,c(1:6),], norm.counts, by="id")
norm.counts <- dplyr::select(norm.counts, c("#chr","start","end","id","info","strand",names(bed)[-c(1:6)]))
norm.counts <- norm.counts[with(norm.counts, order(`#chr`,start,end)),]
write.table(norm.counts, file="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/analysis_DGE/SYSCOL_normal_275_tumor_276.norm.counts.bed", sep="\t", col.names=TRUE,row.names=FALSE,quote=FALSE)

res <- results(dds, alpha = 0.05, parallel=TRUE)
summary(res)


res <- as.data.frame(res)
res$type <- "TE"
res$type[grepl("ENSG", rownames(res))] <- "Gene"
res$sig <- res$padj <= 0.05
res$id <- rownames(res)
res <- res[,c(ncol(res), c(1:ncol(res)-1))]
res <- merge(bed[,c(1:6),],res, by="id")
#res <- dplyr::select(res, c("#chr","start","end","id","info","strand","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","type","sig"))

### COMPUTE MEDIAN EXP VALUES FOR NORMAL AND TUMOR ###

n.samples <- as.character(read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/stat/syscol_normal_samples.txt",header=FALSE,sep="\t", stringsAsFactors = FALSE))
t.samples <- as.character(read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/stat/syscol_tumor_samples.txt",header=FALSE,sep="\t", stringsAsFactors = FALSE))

med <- data.frame()
for (rrr in 1:nrow(norm.counts))
{
  cat(rrr,"-",nrow(norm.counts),"\r")
  x <- norm.counts[rrr,]
  n.median <- median(as.numeric(x[,n.samples]))
  t.median <- median(as.numeric(x[,t.samples]))
  med <- rbind(med, data.frame("id"=x$id, "n_median"=n.median,"t_median"=t.median))
}

res <- merge(res,med, by="id")
res <- dplyr::select(res, c("#chr","start","end","id","info","strand","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","type","sig","n_median","t_median"))

write.table(res, file="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/analysis_DGE/DGE_results_DESEq2_all_TE_genes.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

gg1 <- ggplot(data=res, aes(x=log2FoldChange, y= -log10(padj), colour=type))+
  geom_point()+
  geom_point(data=res[which(res$padj > 0.05),], aes(x=log2FoldChange, y= -log10(padj)), colour="grey")+
  geom_vline(xintercept = 0, lty=2)+
  scale_colour_manual(values=c("Gene"=col$GENE, "TE"=col$TE))+
  theme_classic()+
  theme(legend.position="top",legend.title=element_blank())
ggsave(gg1, filename="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/analysis_DGE/TE_gene_volcano_plot_deseq2.pdf",height=5,width=6.5)


