#!/usr/bin/env Rscript 

library(ggpubr)
library(ggplot2)
library(ggthemes)
setwd("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/revisions/analysis_NeqtlsTEage")

QTL.n <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/mapping/normal/NORMAL.conditional_bestPerRank.TEs.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep=" "))
AGE <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/revisions/analysis_NeqtlsTEage/TE_dfam_diverg.txt",header=TRUE,stringsAsFactors = FALSE))
TMP <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/stat/features_family.bed.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
TMP$family <- sapply(strsplit(TMP$fam, "\\/"),"[[",1)

tab <- as.data.frame(table(QTL.n$phe_id),stringsAsFactors=FALSE)
tab <- tab[order(tab$Freq,decreasing = TRUE),]

tab <- merge(tab,TMP[,c(4:7)],by.x="Var1",by.y="id")

toPlot <- merge(tab, AGE, by.x="name",by.y="Repeat")



tt <- data.frame()
for(n in unique(toPlot$`Divergence_(mya)`)){
  x <- toPlot[which(toPlot$`Divergence_(mya)` == n),]
  mean_eqtls <- sum(x$Freq) / nrow(x)
  xx <- data.frame("Divergence_mya"=n, "N"=nrow(x),"sum"=sum(x$Freq),"mean"=mean_eqtls)
  tt <- rbind(tt,xx)
}
#tt <- merge(tt, unique(toPlot[,c(4,5,6,7,8)]), by="fam")

gg <- ggscatter(data=tt, x="Divergence_mya", y="mean")+
  stat_cor(method="spearman", label.x=200,label.y=1.3,size=5, cor.coef.name = "rho")+
  labs(x="Estimated age of transposable elements [Mya]", y = "Mean number of independent eQTLs")+
  theme_base()+
  theme(axis.title = element_text(face="bold"))
gg
ggsave(gg, filename="Mya_meanEqtls.pdf", height=5.5,width=6,device="pdf")


tt <- data.frame()
for(n in unique(toPlot$fam)){
  x <- toPlot[which(toPlot$fam == n),]
  mean_eqtls <- sum(x$Freq) / nrow(x)
  xx <- data.frame("Divergence_mya"=x$`Divergence_(mya)`, "N"=nrow(x),"sum"=sum(x$Freq),"mean"=mean_eqtls, "fam"=n)
  tt <- rbind(tt,xx)
}
#tt <- merge(tt, unique(toPlot[,c(4,5,6,7,8)]), by="fam")
gg_pearson <- ggscatter(data=tt, x="Divergence_mya", y="mean")+
  stat_cor(method="spearman", label.x=200,label.y=1.5,size=5)+
  labs(x="Estimated age of transposable elements [Mya]", y = "Mean number of independent eQTLs")+
  theme_base()+
  theme(axis.title = element_text(face="bold"))
ggsave(gg_pearson, filename="Mya_meanEqtls_pearson.pdf", height=6,width=6.5,device="pdf")
