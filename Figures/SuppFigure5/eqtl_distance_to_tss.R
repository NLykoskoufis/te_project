
library(qvalue)
library(RColorBrewer)
COL <- RColorBrewer::brewer.pal(8, "Dark2")
col = list("TEs"=COL[3], "Genes"=COL[1])
library(ggplot2)
library(ggpubr)
library(ggthemes)
setwd("/Users/nikolaoslykoskoufis/Documents/PROJECTS/te_project/paper/rerun/Figures/SuppFigure5")
# normal 
n <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/permute/normal/NORMAL.PC30.All.txt.gz",header=T,stringsAsFactors = F,sep=" "))
n <- n[!is.na(n$adj_beta_pval),]
QP.n <- qvalue(n$adj_beta_pval)

n.fdr <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/permute/normal/NORMAL.PC30.All.significant.txt", header=T,stringsAsFactors = F,sep=" "))
n.fdr$type <- "TEs"
n.fdr$type[grepl("ENSG",n.fdr$phe_id)] <- "Genes"

# tumor
t <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/permute/tumor/TUMOR.PC30.All.txt.gz",header=T,stringsAsFactors = F,sep=" "))
t <- t[!is.na(t$adj_beta_pval),]
QP.t <- qvalue(t$adj_beta_pval)

t.fdr <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/permute/tumor/TUMOR.PC30.All.significant.txt", header=T,stringsAsFactors = F,sep=" "))
t.fdr$type <- "TEs"
t.fdr$type[grepl("ENSG",t.fdr$phe_id)] <- "Genes"


#### Distance to tss boxplot ## 

dist_boxplot_norm_log10 <- ggplot(data=n.fdr, aes(x=type, y=log10(abs(dist_phe_var)),fill=type))+
  geom_boxplot(color="black",lwd=1.5,fatten=0.75)+
  labs(y="Absolute distance to TSS (log10 scaled)", title="Normal")+
  stat_compare_means(method="wilcox.test",size= 5,label.x=1.2)+
  theme_base(base_size=20)+
  labs(x="")+
  scale_fill_manual(values = col)+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        legend.position=0,
        plot.title=element_text(hjust=0.5,face="bold"),
        plot.background = element_blank()
        )
ggsave(dist_boxplot_norm_log10,filename="dist_to_tss_normal_boxplotplot_nfdr_TE_vs_GENE_log10scaled.pdf",height=7,width=7,device="pdf")


dist_boxplot_tumor_log10 <- ggplot(data=t.fdr, aes(x=type, y=log10(abs(dist_phe_var)),fill=type))+
  geom_boxplot(color="black",lwd=1.5,fatten=0.75)+
  labs(y="Absolute distance to TSS (log10 scaled)",title="Tumor")+
  stat_compare_means(method="wilcox.test",size= 5,label.x=1.2)+
  theme_base(base_size=20)+
  labs(x="")+
  scale_fill_manual(values = col)+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        legend.position=0,
        plot.title=element_text(hjust=0.5,face="bold"),
        plot.background = element_blank()
  )
ggsave(dist_boxplot_tumor_log10,filename="dist_to_tss_tumor_boxplotplot_tfdr_TE_vs_GENE_log10scaled.pdf",height=7,width=7,device="pdf")

