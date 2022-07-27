#!/usr/bin/env Rscript 

setwd("/Users/nikolaoslykoskoufis/Documents/PROJECTS/te_project/paper/rerun/Figures/SuppFigure7")

library(ggplot2)
library(ggpubr)
library(RColorBrewer)
COL = brewer.pal(9,"Set1")

tissueSpe <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/specific/raw/tissueSpecific_eQTLs.txt",header=T,stringsAsFactors=F,sep=" ")
shared <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/specific/raw/shared_eQTLs.txt",header=T,stringsAsFactors=F,sep=" ")

gg_gene <- ggplot()+
    geom_point(data=tissueSpe[grepl("ENSG",tissueSpe$phe_id),], aes(x=dist_phe_var.y,y=-log10(bwd_pval),colour="tissue-specific"))+
    geom_point(data=shared[grepl("ENSG",shared$phe_id),], aes(x=dist_phe_var,y=-log10(t_eqtl_bwd_pval),colour="shared"))+
    labs(x="distance (bp)",y="-log10(P-value)")+
    theme_classic(base_size=20)+
    scale_colour_manual(values=c("tissue-specific"="red","shared"="black"))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.position=c(0.8,0.9),
          legend.title=element_blank(),
          legend.box.background = element_rect(color="black", size=0.5),
          legend.box.margin = margin(2, 4, 4, 4),
          plot.title=element_text(hjust=0.5,face="bold"))
ggsave(gg_gene, filename="gene_tissue_specific_shared_eQTLs_dist.pdf",device="pdf",height=8,width=8)




### Diff in distance between shared and tissueSPe for TE-eQTLs
wilcox.test(abs(shared[grepl("ENSG",shared$phe_id),]$dist_phe_var),abs(tissueSpe[!grepl("ENSG",tissueSpe$phe_id),]$dist_phe_var.y)) 
# Wilcox rank sum test P<2.2e-16
median(abs(tissueSpe[grepl("ENSG",shared$phe_id),]$dist_phe_var))
# median abs shared distance gene: 11,671
median(abs(tissueSpe[grepl("ENSG",tissueSpe$phe_id),]$dist_phe_var.x))
# median abs tissueSPe distance gene: 39,981

## mosaic plot 

t_gene <- tissueSpe[which(tissueSpe$tissue == "tumor"),]
t_gene <- t_gene[grepl("ENSG",t_gene$phe_id),]
n_gene <- tissueSpe[which(tissueSpe$tissue == "normal"),]
n_gene <- n_gene[grepl("ENSG",n_gene$phe_id),]

data2 <- matrix(c(nrow(t_gene),nrow(n_gene),nrow(shared[grepl("ENSG",shared$phe_id),]),nrow(shared[grepl("ENSG",shared$phe_id),])),byrow=T,2,2)
colnames(data2) <- c("Tumor","Normal")
rownames(data2) <- c("Unique to cell type","shared")
pdf("Gene_eqtl_level_sharing_mosaic_plot.pdf",8,8)
mosaicplot(t(data2),cex.axis=1.5,main="")
dev.off()


#### effect size plot #### 


gg_effectSize_gene <- ggplot()+
    geom_point(data=shared[grepl("ENSG",shared$phe_id),], aes(x=n_eqtl_bwd_slope,y=t_eqtl_bwd_slope,colour="Shared"))+
    geom_point(data=t_gene, aes(x=nominal_slope, y=bwd_slope, colour="Tumor specific"))+
    geom_point(data=n_gene, aes(x=bwd_slope, y=nominal_slope, colour="Normal specific"))+
    labs(x="Normal Gene-eQTL slope",y="Tumor Gene-eQTL slope")+
    theme_classic()+
    scale_colour_manual(values=c("Shared"="black","Tumor specific"="#e41a1c","Normal specific"="#377eb8"))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.position="top",
    legend.direction="vertical",
    legend.justification="left",
    legend.title=element_blank())
ggsave(gg_effectSize_gene, filename="Gene_tissue_specific_shared_slopes.pdf",device="pdf",height=5, width=5.5)



##### Methylation levels ##### 

tspe.gene.meth <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/methylation/tspe_gene_meth.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t"))
shared.gene.meth <- as.data.frame(data.table::fread("//Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/methylation/shared_gene_meth.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t"))

gene <- rbind(tspe.gene.meth[which(tspe.gene.meth$qval <= 0.05),], shared.gene.meth[which(shared.gene.meth$qval <= 0.05),])


gg_gene <- ggplot(data=gene, aes(x=type, y=abs_med_dif,fill=type))+
    geom_boxplot(color="black",lwd=1.15,fatten=0.75,outlier.color="black")+
    stat_compare_means(method="wilcox.test",size=7)+
    labs(x="",y="Absolute value difference between\nnormal and tumor medians")+
    scale_fill_manual(values=c("tumor-specific"=COL[1],"shared"=COL[2]))+
    theme_classic(base_size=20)+
    theme(plot.title=element_text(hjust=0.5,size=18,face="bold"),
          legend.position=0)
ggsave(gg_gene, filename="absolute_value_difference_methylation_tumorSpecific_gene_eQTLs.pdf",device="pdf",height=8,width=8)
