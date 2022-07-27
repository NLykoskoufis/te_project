#/usr/bin/env Rscript 

library(gplots)
library(qvalue)
library(RColorBrewer)
COL <- RColorBrewer::brewer.pal(8, "Dark2")
col = list("TEs"=COL[3], "Genes"=COL[1])
library(ggplot2)
library(ggpubr)

setwd("~/Documents/PROJECTS/te_project/paper/rerun/Figures/Figure2/")

# Distance to TSS plots - Figure 1A-B 

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


## DISTANCE TO TSS ## 

cor.test(-log10(n.fdr[which(n.fdr$type == "TEs"),]$adj_beta_pval),abs(n.fdr[which(n.fdr$type == "TEs"),]$dist_phe_var),method="spearman")
cor.test(-log10(t.fdr[which(t.fdr$type == "TEs"),]$adj_beta_pval),abs(t.fdr[which(t.fdr$type == "TEs"),]$dist_phe_var),method="spearman")



dist_tss_norm <- ggplot(data=n.fdr, aes(x=dist_phe_var, y=-log10(adj_beta_pval),colour=type))  +
  geom_point()+
  labs(x="distance to TSS (bp)", y="-log10(pvalue)",title="eQTLs in normal")+
  theme_bw(base_size=20)+
  scale_colour_manual(values = col)+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        legend.title=element_blank(),
        legend.position=c(0.9,0.9),legend.background = element_rect(colour="black",size=0.5),legend.margin=margin(l=0.1,t=0,r=0.1,b=0.1, unit='cm'),
        plot.title=element_text(hjust=0.5,face="bold"))

ggsave(dist_tss_norm,filename="Figure1A.pdf",height=7,width=7.5,device="pdf")


dist_tss_tumor <- ggplot(data=t.fdr, aes(x=dist_phe_var, y=-log10(adj_beta_pval),colour=type))  +
  geom_point()+
  labs(x="distance to TSS (bp)", y="-log10(pvalue)",title="eQTLs in tumor")+
  theme_bw(base_size=20)+
  scale_colour_manual(values = col)+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        legend.title=element_blank(),
        legend.position=c(0.9,0.9),legend.background = element_rect(colour="black",size=0.5),legend.margin=margin(l=0.1,t=0,r=0.1,b=0.1, unit='cm'),
        plot.title=element_text(hjust=0.5,face="bold")
        )

ggsave(dist_tss_tumor,filename="figure2B.pdf",height=7,width=7.5,device="pdf")

## Independent eQTLs - Figure 2C-D

tumor <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/mapping/tumor/TUMOR.conditional_bestPerRank.txt.gz",header=T,stringsAsFactors=F))

normal <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/mapping/normal/NORMAL.conditional_bestPerRank.txt.gz",header=T,stringsAsFactors=F,sep=" "))

COL <- RColorBrewer::brewer.pal(8, "Dark2")
col = list("TE"=COL[3], "GENE"=COL[1])


N.te <- normal[!grepl("ENSG",normal$phe_id),]
N.gene <- normal[grepl("ENSG",normal$phe_id),]



df<-merge(data.frame(table(N.gene$rank)),data.frame(table(N.te$rank)),by=1,all=T)
rownames(df) <- df$Var1
df<-df[-1]
names(df) <- c("genes","TEs")
df[,1] <- df[,1] / df[1,1]
df[,2] <- df[,2] / df[1,2]
df<-df[-1,]

pdf(file="figure2CD.pdf",bg="white",8,4.5)
par(mfrow=c(1,2))
barplot(as.matrix(t(df)), main = "Independent eQTLs in normal", xlab = "Number of secondary eQTLs", ylab = "Proportion of genes/TEs with secondary eQTLs", beside = T,col=c(col$GENE, col$TE),legend.text =T)

T.te <- tumor[!grepl("ENSG",tumor$phe_id),]
T.gene <- tumor[grepl("ENSG",tumor$phe_id),]


df.t<-merge(data.frame(table(T.gene$rank)),data.frame(table(T.te$rank)),by=1,all=T)
rownames(df.t) <- df.t$Var1
df.t<-df.t[-1]
names(df.t) <- c("genes","TEs")
df.t[,1] <- df.t[,1] / df.t[1,1]
df.t[,2] <- df.t[,2] / df.t[1,2]
df.t<-df.t[-1,]

barplot(as.matrix(t(df.t)), main = "Independent eQTLs in tumor", xlab = "Number of secondary eQTLs", ylab = "Proportion of genes/TEs with secondary eQTLs", beside = T,col=c(col$GENE, col$TE),legend.text =T)
dev.off()



