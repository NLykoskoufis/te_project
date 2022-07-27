#!/usr/bin/env Rscript 

library(ggplot2)
library(ggpubr)


tissueSpe <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/specificity/tissueSpecific_eQTLs.txt.gz",header=T,stringsAsFactors=F,sep="\t")
shared <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/specificity/shared_eQTLs.txt.gz",header=T,stringsAsFactors=F,sep=" ")

gg_te <- ggplot()+
  geom_point(data=tissueSpe[!grepl("ENSG",tissueSpe$phe_id),], aes(x=dist_phe_var.y,y=-log10(bwd_pval),colour="tissue-specific"))+
  geom_point(data=shared[!grepl("ENSG",shared$phe_id),], aes(x=dist_phe_var,y=-log10(t_eqtl_bwd_pval),colour="shared"))+
  #geom_point(data=tissueSpe[grepl("ENSG",tissueSpe$phe_id),], aes(x=dist_phe_var,y=log10(bwd_pval),colour="tissue-specific",shape="gene"))+
  #geom_point(data=shared[grepl("ENSG",shared$phe_id),], aes(x=dist_phe_var,y=log10(t_eqtl_bwd_pval),colour="shared",shape="gene"))+
  labs(x="distance (bp)",y="-log10(P-value)")+
  theme_classic(base_size=20)+
  scale_colour_manual(values=c("tissue-specific"="red","shared"="black"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.position=c(0.8,0.9),
        legend.title=element_blank(),
        legend.box.background = element_rect(color="black", size=0.5),
        legend.box.margin = margin(2, 4, 4, 4),
        plot.title=element_text(hjust=0.5,face="bold"))
ggsave(gg_te, filename="TE_tissue_specific_shared_eQTLs_dist.pdf",device="pdf",height=8,width=8)


### Diff in distance between shared and tissueSpe?
wilcox.test(abs(shared$dist_phe_var),abs(tissueSpe$dist_phe_var)) 
#Wilcox rank sum test Pvalue < 2.2e-16
# median abs shared distance 6'557 bp
# median abs tissueSpe distance 51'868 bp

### Diff in distance between shared and tissueSPe for TE-eQTLs
wilcox.test(abs(shared[!grepl("ENSG",shared$phe_id),]$dist_phe_var),abs(tissueSpe[!grepl("ENSG",tissueSpe$phe_id),]$dist_phe_var.y)) 
# Wilcox rank sum test P<2.2e-16
# median abs shared distance TE: 5'580.5
# median abs tissueSPe distance TE: 60'870

## mosaic plot 

### mosaic plot for TE-eQTLs 

t_te <- tissueSpe[which(tissueSpe$tissue == "tumor"),]
t_te <- t_te[!grepl("ENSG",t_te$phe_id),]

n_te <- tissueSpe[which(tissueSpe$tissue == "normal"),]
n_te <- n_te[!grepl("ENSG",n_te$phe_id),]


data1 <- matrix(c(nrow(t_te),nrow(n_te),nrow(shared[!grepl("ENSG",shared$phe_id),]),nrow(shared[!grepl("ENSG",shared$phe_id),])),byrow=T,2,2)
colnames(data1) <- c("Tumor","Normal")
rownames(data1) <- c("Unique to cell type","shared")
pdf("TE_eqtl_level_sharing_mosaic_plot.pdf",8,8)
mosaicplot(t(data1),cex.axis=1.5,main="")
dev.off()

#                    Tumor Normal
#Unique to cell type   376   1685
#shared                524    524




#### effect size plot #### 
gg_effectSize_te <- ggplot()+
  geom_point(data=shared[!grepl("ENSG",shared$phe_id),], aes(x=n_eqtl_bwd_slope,y=t_eqtl_bwd_slope,colour="Shared"))+
  geom_point(data=t_te, aes(x=nominal_slope, y=bwd_slope, colour="Tumor specific"))+
  geom_point(data=n_te, aes(x=bwd_slope, y=nominal_slope, colour="Normal specific"))+
  labs(x="Normal TE-eQTL slope",y="Tumor TE-eQTL slope")+
  theme_classic()+
  scale_colour_manual(values=c("Shared"="black","Tumor specific"="#e41a1c","Normal specific"="#377eb8"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.position="top",
        legend.direction="vertical",
        legend.justification="left",
        legend.title=element_blank())
ggsave(gg_effectSize_te, filename="TE_tissu e_specific_shared_slopes.pdf",device="pdf",height=5, width=5.5)


##### Methylation levels ##### 

tspe.te.meth <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_methylation/tspe_te_meth.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t"))
shared.te.meth <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_methylation/shared_te_meth.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t"))

te <- rbind(tspe.te.meth[which(tspe.te.meth$qval <= 0.05),], shared.te.meth[which(shared.te.meth$qval <= 0.05),])

gg_te <- ggplot(data=te, aes(x=type, y=abs_med_dif,fill=type))+
  geom_boxplot(color="black",lwd=1.15,fatten=0.75,outlier.color="black")+
  stat_compare_means(method="wilcox.test",size=7)+
  labs(x="",y="Absolute value difference between\nnormal and tumor medians")+
  scale_fill_manual(values=c("tumor-specific"=COL[1],"shared"=COL[2]))+
  theme_classic(base_size=20)+
  theme(plot.title=element_text(hjust=0.5,size=18,face="bold"),
        legend.position=0)
ggsave(gg_te,filename="absolute_value_difference_methylation_tumorSpecific_TE_eQTLs.pdf",device="pdf",height=8,width=8)

