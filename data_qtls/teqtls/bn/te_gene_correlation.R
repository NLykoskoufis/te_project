#!/usr/bin/env Rscript 


suppressMessages(source("/home/users/l/lykoskou/bin/tools.R"))
suppressMessages(library(ggplot2))
library(gridExtra)
library(qvalue)
library(ggpubr)
library(RColorBrewer)
COL = brewer.pal(9,"Set1")
COL2 = brewer.pal(9,"Set3")


PERM = as.data.frame(data.table::fread("zcat < /srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/emods/nominal/tumor/TUMOR.TE.GENE.nominal_All.txt.gz",header=T,stringsAsFactors = F,sep=" "))
PERM$key <- paste(PERM$phe_id, PERM$var_id, sep="_")

BN_TUMOR <- nfread("/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/teqtls/bn/tumor/TUMOR.BNlearn.ALL.extraINFO.txt.gz")
BN_TUMOR$key <- paste(BN_TUMOR$gene, BN_TUMOR$te,sep="_")

TMP <- merge(PERM,BN_TUMOR, by="key")
TMP$best <- gsub("m1","Causal",TMP$best)
TMP$best <- gsub("m2","Reactive",TMP$best)
TMP$best <- gsub("m3","Independent",TMP$best)

gg_t <- ggplot(data=TMP, aes(x=slope))+
    geom_histogram(breaks=brx(TMP$slope),color="black",fill="white")+
    labs(x="Effect Size (regression slope)",y="Counts",title="TE-Gene association in tumor")+
    theme_bw(base_size=20)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.title=element_text(hjust=0.5,face="bold"),
          legend.position=0)

ggsave(gg_t, filename="TUMOR.histogram_te_gene_association_slopes.pdf",device="pdf",height=5.5,width=6)


PERM.N = nfread("/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/emods/nominal/normal/NORMAL.TE.GENE.nominal_All.txt.gz",separator=" ")
PERM.N$key <- paste(PERM.N$phe_id, PERM.N$var_id, sep="_")

pdf("te_gene_association_slope.pdf",8,5)
par(mfrow=c(1,2))
hist(PERM.N$slope,breaks=100,main="TE-gene association in normal",xlab="Effect size (regression slope)")
hist(PERM$slope,breaks=100,main="TE-gene association in tumor",xlab="Effect size (regression slope)")
dev.off()

BN_NORMAL <- nfread("/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/teqtls/bn/normal/NORMAL.BNlearn.ALL.extraINFO.txt")
BN_NORMAL$key <- paste(BN_NORMAL$gene, BN_NORMAL$te,sep="_")

TMP.n <- merge(PERM.N,BN_NORMAL, by="key")
TMP.n$best <- gsub("m1","Causal",TMP.n$best)
TMP.n$best <- gsub("m2","Reactive",TMP.n$best)
TMP.n$best <- gsub("m3","Independent",TMP.n$best)

gg_n <- ggplot(data=TMP, aes(x=slope))+
    geom_histogram(breaks=brx(TMP$slope),color="black",fill="white")+
    labs(x="Effect Size (regression slope)",y="Counts",title="TE-Gene association in normal")+
    theme_bw(base_size=20)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.title=element_text(hjust=0.5,face="bold"),
          legend.position=0)

ggsave(gg_n, filename="NORMAL.histogram_te_gene_association_slopes.pdf",device="pdf",height=5.5,width=10)


TMP$tissue <- "Tumor"
TMP.n$tissue <- "Normal"

df <- rbind(TMP, TMP.n)

gg<- ggplot(data=df, aes(x=slope, fill=tissue))+
    geom_histogram(breaks=brx(df$slope),alpha=0.7,color="black")+
    facet_wrap(~best)+
    labs(x="",y="Effect size (slope)\nabsolute value ",title="TE-GENE association")+
    scale_fill_manual(values=c("Normal"=COL[2],"Tumor"=COL[1]))+
    theme_bw(base_size=20)+
     theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.title=element_text(hjust=0.5,face="bold"),
          legend.position="top",legend.title=element_blank())

ggsave(gg, filename="te_gene_slopes_per_model.pdf",device="pdf",height=5,width=12)




### SHARED ### 

SHARED <- nfread("/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/teqtls/shared/NORMAL_TUMOR_shared_BNlearn.txt.gz")
SHARED$key <- paste(SHARED$gene, SHARED$te, sep="_")


PERM = nfread("/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/emods/nominal/tumor/TUMOR.TE.GENE.nominal_All.txt.gz",separator=" ")
PERM$key <- paste(PERM$phe_id, PERM$var_id, sep="_")

PERM.N = nfread("/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/emods/nominal/normal/NORMAL.TE.GENE.nominal_All.txt.gz",separator=" ")
PERM.N$key <- paste(PERM.N$phe_id, PERM.N$var_id, sep="_")


PERM.T <- PERM[,c(12,14,16)]
PERM.T$tissue <- "tumor"
PERM.N <- PERM.N[,c(12,14,16)]
PERM.N$tissue <- "normal"


TMP <- merge(merge(SHARED,PERM.T, by="key"),PERM.N, by="key")
TMP$tumor_best <- gsub("m1","Causal",TMP$tumor_best)
TMP$tumor_best <- gsub("m2","Reactive",TMP$tumor_best)
TMP$tumor_best <- gsub("m3","Independent",TMP$tumor_best)
TMP$normal_best <- gsub("m1","Causal",TMP$normal_best)
TMP$normal_best <- gsub("m2","Reactive",TMP$normal_best)
TMP$normal_best <- gsub("m3","Independent",TMP$normal_best)
TMP$switch <- paste(TMP$normal_best, TMP$tumor_best, sep="_")

test <- TMP[,c(1,24,27,29)]
tmp <- melt(test, id.vars=c("key","switch"))
tmp$variable <- gsub("slope.x","tumor",tmp$variable)
tmp$variable <- gsub("slope.y","normal",tmp$variable)

ggt <- ggplot(data=TMP, aes(x=slope.x))+
    geom_histogram(breaks=brx(TMP$slope.x),color="black",fill="white")+
    labs(x="Effect Size (regression slope)",y="Counts",title="TE-Gene association in tumor")+
    theme_classic(base_size=20)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.title=element_text(hjust=0.5,face="bold"),
          legend.position=0)
ggsave(ggt,filename="SHARED.te_gene_effectSize_tumor.pdf",device="pdf",height=5.5,width=6)

ggn <- ggplot(data=TMP, aes(x=slope.y))+
    geom_histogram(breaks=brx(TMP$slope.y),color="black",fill="white")+
    labs(x="Effect Size (regression slope)",y="Counts",title="TE-Gene association in normal")+
    theme_classic(base_size=20)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.title=element_text(hjust=0.5,face="bold"),
          legend.position=0)
ggsave(ggn,filename="SHARED.te_gene_effectSize_normal.pdf",device="pdf",height=5.5,width=6)


gg <- ggplot(data=tmp, aes(x=variable,y=abs(value), fill=variable))+
    geom_boxplot()+
    stat_compare_means(method="wilcox.test")+
    labs(x="",y="Effect size (slope)\nabsolute value ",title="TE-GENE association\nshared triplets (N=1592)")+
    scale_fill_manual(values=c("normal"=COL[2],"tumor"=COL[1]))+
    theme_bw(base_size=20)+
     theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.title=element_text(hjust=0.5,face="bold"),
          legend.position=0)

ggsave(gg, filename="SHARED.boxplot_te_gene_slopes_n_vs_t.pdf",device="pdf",height=5.5,width=6)

gg_dens <- ggplot(data=tmp, aes(x=value, fill=variable))+
    geom_density(alpha=0.6, color="black")+
    facet_wrap(~tumor_best)+
    labs(x="Effect size (slope)\nabsolute value ",title="TE-GENE association\nshared triplets (N=1592)")+
    scale_fill_manual(values=c("normal"=COL[2],"tumor"=COL[1]))+
    theme_bw(base_size=20)+
     theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.title=element_text(hjust=0.5,face="bold"),
          legend.position=0)
ggsave(gg_dens, filename="SHARED.density_plot_effect_size_models.pdf",device="pdf",height=5.5,width=6)



gg_shared_switch <- ggplot(data=tmp, aes(x=variable,y=abs(value), fill=variable))+
    geom_boxplot()+
    facet_wrap(~switch)+
    stat_compare_means(method="wilcox.test",size=5)+
    labs(x="",y="Effect size (slope)\nabsolute value ",title="TE-GENE association\nshared triplets (N=1592)")+
    scale_fill_manual(values=c("normal"=COL[2],"tumor"=COL[1]))+
    theme_bw(base_size=20)+
     theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.title=element_text(hjust=0.5,face="bold"),
          legend.position=0)
ggsave(gg_shared_switch,filename="SHARED.boxplot_te_gene_slope_switches.pdf",device="pdf",height=16.5, width=18)


tmp$tumor_best <- sapply(strsplit(tmp$switch,"\\_"),"[[",2)
tmp$normal_best <- sapply(strsplit(tmp$switch,"\\_"),"[[",1)


gg_tumor <- ggplot(data=TMP, aes(x=slope.x, fill=tumor_best))+
    geom_histogram(breaks=brx(TMP$slope.x),color="black")+ 
    facet_wrap(~tumor_best)+
    labs(x="Effect size (slope)\nabsolute value ",title="TE-GENE association\nshared triplets (N=1592) in tumor")+
    scale_fill_manual(values=c("Causal"=COL2[3],"Reactive"=COL2[4],"Independent"=COL2[6]))+
    theme_bw(base_size=20)+
     theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.title=element_text(hjust=0.5,face="bold"))
ggsave(gg_tumor, filename="SHARED.tumor_tegene_slope_permodel.pdf",device="pdf",height=5.5,width=18)

gg_normal <- ggplot(data=TMP, aes(x=slope.y, fill=normal_best))+
    geom_histogram(breaks=brx(TMP$slope.y),color="black")+ 
    facet_wrap(~normal_best)+
    labs(x="Effect size (slope)\nabsolute value ",title="TE-GENE association\nshared triplets (N=1592) in normal")+
    scale_fill_manual(values=c("Causal"=COL2[3],"Reactive"=COL2[4],"Independent"=COL2[6]))+
    theme_bw(base_size=20)+
     theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.title=element_text(hjust=0.5,face="bold"))
ggsave(gg_normal, filename="SHARED.normal_tegene_slope_permodel.pdf",device="pdf",height=5.5,width=18)


### UNION ### 

UNION <- nfread("/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/teqtls/shared/NORMAL.TUMOR_ALL_union.txt.gz")
UNION$key <- paste(UNION$gene, UNION$te, sep="_")


PERM = nfread("/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/emods/nominal/tumor/TUMOR.TE.GENE.nominal_All.txt.gz",separator=" ")
PERM$key <- paste(PERM$phe_id, PERM$var_id, sep="_")

PERM.N = nfread("/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/emods/nominal/normal/NORMAL.TE.GENE.nominal_All.txt.gz",separator=" ")
PERM.N$key <- paste(PERM.N$phe_id, PERM.N$var_id, sep="_")



TMP <- merge(merge(UNION,PERM.T, by="key"),PERM.N, by="key")
TMP$tumor_best <- gsub("m1","Causal",TMP$tumor_best)
TMP$tumor_best <- gsub("m2","Reactive",TMP$tumor_best)
TMP$tumor_best <- gsub("m3","Independent",TMP$tumor_best)
TMP$normal_best <- gsub("m1","Causal",TMP$normal_best)
TMP$normal_best <- gsub("m2","Reactive",TMP$normal_best)
TMP$normal_best <- gsub("m3","Independent",TMP$normal_best)
TMP$switch <- paste(TMP$normal_best, TMP$tumor_best, sep="_")

test <- TMP[,c(1,9,24,26)]
tmp <- melt(test, id.vars=c("key","switch"))
tmp$variable <- gsub("slope.x","tumor",tmp$variable)
tmp$variable <- gsub("slope.y","normal",tmp$variable)

ggt <- ggplot(data=TMP, aes(x=slope.x))+
    geom_histogram(breaks=brx(TMP$slope.x),color="black",fill="white")+
    labs(x="Effect Size (regression slope)",y="Counts",title="TE-Gene association in tumor")+
    theme_classic(base_size=20)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.title=element_text(hjust=0.5,face="bold"),
          legend.position=0)
ggsave(ggt,filename="UNION.te_gene_effectSize_tumor.pdf",device="pdf",height=5.5,width=6)

ggn <- ggplot(data=TMP, aes(x=slope.y))+
    geom_histogram(breaks=brx(TMP$slope.y),color="black",fill="white")+
    labs(x="Effect Size (regression slope)",y="Counts",title="TE-Gene association in normal")+
    theme_classic(base_size=20)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.title=element_text(hjust=0.5,face="bold"),
          legend.position=0)
ggsave(ggn,filename="UNION.te_gene_effectSize_normal.pdf",device="pdf",height=5.5,width=6)



gg <- ggplot(data=tmp, aes(x=variable,y=abs(value), fill=variable))+
    geom_boxplot()+
    stat_compare_means(method="wilcox.test",size=7)+
    labs(x="",y="Effect size (slope)\nabsolute value ",title=paste("TE-GENE association\nunion triplets (N=",nrow(UNION),")",sep=""))+
    scale_fill_manual(values=c("normal"=COL[2],"tumor"=COL[1]))+
    theme_bw(base_size=20)+
     theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.title=element_text(hjust=0.5,face="bold"),
          legend.position=0)

ggsave(gg, filename="UNION.boxplot_te_gene_slopes_n_vs_t.pdf",device="pdf",height=5.5,width=6)

gg_model_union <- ggplot(data=TMP, aes(x=slope.x,fill=tumor_best))+
    geom_histogram(breaks=brx(TMP$slope.x),color="black")+
    facet_wrap(~tumor_best)+
    labs(x="Effect Size (regression slope)",y="Counts",title="TE-Gene association in tumor (Union triplets)")+
    theme_bw(base_size=20)+
    scale_fill_manual(values=c("Causal"=COL2[3],"Reactive"=COL2[4],"Independent"=COL2[6]))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.title=element_text(hjust=0.5,face="bold"),
          legend.position=0)
ggsave(gg_model_union, filename="UNION.tumor_tegene_slope_permodel.pdf",device="pdf",height=5.5,width=18)
