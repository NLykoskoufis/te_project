#!/usr/bin/env Rscript 


library(qvalue)
library(ggplot2)
library(dplyr)
library(foreach)
library(RColorBrewer)
COL = brewer.pal(9,"Set1")
library(doMC)
library(ggpubr)

TSPE <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/specific/raw/tissueSpecific_eQTLs.txt",sep=" ",header=T,stringsAsFactors=F))
SHARED <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/specific/raw/shared_eQTLs.txt",sep=" ", header=T,stringsAsFactors=F))

TSPE <- TSPE[which(TSPE$spe_hit == 1),]
TSPE.gene <- TSPE[grepl("ENSG",TSPE$phe_id),]
TSPE.te <- TSPE[!grepl("ENSG",TSPE$phe_id),]

SHARED.gene <- SHARED[grepl("ENSG",SHARED$phe_id),]
SHARED.te <- SHARED[!grepl("ENSG",SHARED$phe_id),]

ME <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_methylation/data_intersection/syscol_variants_methylations_2.5kbwindow.txt.gz",header=FALSE,stringsAsFactors=F,sep="\t"))
colnames(ME) <- c("var_chr","var_start","var_end","var_id","var_strand","meth_chr","meth_start","meth_end","meth_id","dummy","meth_strand")
tmp <- dplyr::select(ME, c("var_id","meth_id"))
tmp <- aggregate(meth_id ~ var_id, tmp, paste, collapse = ",")

bed <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_methylation/data_methylation/syscol_Norm_Betas.3-6-2015.bed.gz",header=T,stringsAsFactors=F,sep="\t"))

#### TESTING FOR ONE CASE ####
t.samples <- as.character(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/stat/syscol_tumor_samples.txt",header=F,stringsAsFactors=F,sep="\t"))[-c(1:6)]


n.samples <- as.character(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/stat/syscol_normal_samples.txt",header=F,stringsAsFactors=F,sep="\t"))[-c(1:6)]


methanalysis <- function(file){
    df <- merge(file,tmp,by="var_id")
    df$key <- paste(df$phe_id, df$var_id, sep=";")

    res <- foreach(rrr=unique(df$key), .combine = rbind, .multicombine = T) %do% {
    cat(rrr,"\n")
    x <- df[which(df$key == rrr),]
    meth <- unlist(strsplit(x$meth_id,"\\,"))

    d <- bed[bed$id %in% meth,]
    info <- d[,c(1:6)]
    n <- names(d)[grepl("N",names(d))]
    t <- names(d)[!grepl("N",names(d))][-c(1:6)]
    n <- intersect(n,n.samples)
    t <- intersect(t,t.samples)

    normal <- d %>% select(n)
    tumor <- d %>% select(t)
    
    n.med <- as.numeric(apply(normal,2,median,na.rm=TRUE)) 
    t.med <- as.numeric(apply(tumor, 2, median,na.rm=TRUE))
    test <- try(wilcox.test(t.med,n.med))

    if(class(test) == "try-error"){
        pval <- NA
    }else{
        pval <- test$p.value
    }
    return(data.frame("var_id"=x$var_id,"phe_id"=x$phe_id,"wilcox_pval"=pval,"abs_med_dif"=abs(median(n.med,na.rm=TRUE)-median(t.med,na.rm=TRUE)),"median_normal"=median(n.med,na.rm=TRUE),"median_tumor"=median(t.med,na.rm=TRUE),"n_med"=paste(n.med,collapse=","),"t_med"=paste(t.med,collapse=",")))
}
    res$qval <- qvalue(res$wilcox_pval)$qvalues
    return(res)
}

tspe.gene.meth <- methanalysis(TSPE.gene)
tspe.te.meth <- methanalysis(TSPE.te)
shared.gene.meth <- methanalysis(SHARED.gene)
shared.te.meth <- methanalysis(SHARED.te)

tspe.gene.meth$type <- "tumor-specific"
tspe.te.meth$type <- "tumor-specific"
shared.gene.meth$type <- "shared"
shared.te.meth$type <- "shared"

data.table::fwrite(tspe.gene.meth, file="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/methylation/tspe_gene_meth.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
data.table::fwrite(tspe.te.meth, file="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/methylation/tspe_te_meth.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
data.table::fwrite(shared.gene.meth, file="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/methylation/shared_gene_meth.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
data.table::fwrite(shared.te.meth, file="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/methylation/shared_te_meth.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

ALL <- rbind(rbind(rbind(tspe.gene.meth, tspe.te.meth),shared.gene.meth),shared.te.meth)
data.table::fwrite(ALL, file="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/methylation/methylation_results_specificityeQTLs.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)


pdf("pvalue_qvalue_distribution_plots.pdf",8,8)
hist(tspe.gene.meth$wilcox_pval,main="tumor-specific gene methylation\npvalue distribution")
plot(qvalue(tspe.gene.meth$wilcox_pval),main="tumor-specific gene methylation\nqvalue")
hist(tspe.te.meth$wilcox_pval,main="tumor-specific TE methylation\npvalue distribution")
plot(qvalue(tspe.te.meth$wilcox_pval),main="tumor-specific TE methylation\nqvalue")
hist(shared.gene.meth$wilcox_pval,main="shared gene methylation\npvalue distribution")
plot(qvalue(shared.gene.meth$wilcox_pval),main="shared gene methylation\nqvalue")
hist(shared.te.meth$wilcox_pval,main="shared TE methylation\npvalue distribution")
plot(qvalue(shared.te.meth$wilcox_pval),main="shared gene methylation\nqvalue")
dev.off()

tspe.gene.meth.fdr <- tspe.gene.meth[which(tspe.gene.meth$qval <= 0.05),]
tspe.te.meth.fdr <- tspe.te.meth[which(tspe.te.meth$qval <= 0.05),]
shared.gene.meth <- shared.gene.meth[which(shared.gene.meth$qval <= 0.05),]
shared.te.meth <- shared.te.meth[which(shared.te.meth$qval <= 0.05),]


gene <- rbind(tspe.gene.meth.fdr,shared.gene.meth)
te <- rbind(tspe.te.meth.fdr, shared.te.meth)

gg_gene <- ggplot(data=gene, aes(x=type, y=abs_med_dif,fill=type))+
    geom_boxplot(color="black",lwd=1.15,fatten=0.75,outlier.color="black")+
    stat_compare_means(method="wilcox.test",size=7)+
    labs(x="",y="Absolute value difference between\nnormal and tumor medians",title="Distribution of absolute value difference\nof median methylation betas in normal vs. tumors\n for gene eQTLs")+
    scale_fill_manual(values=c("tumor-specific"=COL[1],"shared"=COL[2]))+
    theme_classic(base_size=20)+
    theme(plot.title=element_text(hjust=0.5,size=18,face="bold"),
        legend.position=0)
ggsave(gg_gene, filename="absolute_value_difference_methylation_tumorSpecific_gene_eQTLs.pdf",device="pdf",height=8,width=8)

gg_te <- ggplot(data=te, aes(x=type, y=abs_med_dif,fill=type))+
    geom_boxplot(color="black",lwd=1.15,fatten=0.75,outlier.color="black")+
    stat_compare_means(method="wilcox.test",size=7)+
    labs(x="",y="Absolute value difference between\nnormal and tumor medians",title="Distribution of absolute value difference\nof median methylation betas in normal vs. tumors\n for TE eQTLs")+
    scale_fill_manual(values=c("tumor-specific"=COL[1],"shared"=COL[2]))+
    theme_classic(base_size=20)+
    theme(plot.title=element_text(hjust=0.5,size=18,face="bold"),
        legend.position=0)
ggsave(gg_te,filename="absolute_value_difference_methylation_tumorSpecific_TE_eQTLs.pdf",device="pdf",height=8,width=8)



te$difNT <- te$median_normal > te$median_tumor
gene$difNT <- gene$median_normal > gene$median_tumor
