#!/usr/bin/env Rscript 

library(RColorBrewer)
library(ggplot2)
library(dplyr)



##========================================================##
##========= CHECKING TUMOR SPECIFIC AND TRIPLETS =========##
##========================================================##


UNION <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/union/UNION_normal_tumor_with_posteriors.txt.gz",header=T,stringsAsFactors=F,sep="\t"))
UNION$switch <- "switch"
UNION[which(UNION$normal_best == UNION$tumor_best),]$switch <- "no_switch"
UNION[which(UNION$normal_best != UNION$tumor_best & UNION$tumor_best == "m1"),]$switch <- "switch_causal"


UNION$sub <- paste(UNION$normal_best, UNION$tumor_best,sep="_")


bn <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/tumor/TUMOR.BNlearn.ALL.extraINFO.txt.gz",header=TRUE,stringsAsFactors=F,sep="\t"))
tspe <- read.table("//Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/specific/tumor_specific/TUMOR.specific_TE_eQTLs.txt.gz",header=T,stringsAsFactors = F,sep=" ")

bn$triplet <- paste(bn$variant, bn$gene, bn$te, sep=";")
bn$key <- paste(bn$variant, bn$te,sep="_")
tspe$key <- paste(tspe$var_id, tspe$phe_id,sep="_")
UNION$triplet <- paste(UNION$var_id, UNION$gene, UNION$te, sep=";")
d.bn <- merge(bn,tspe, by="key")


d.un <- UNION[UNION$triplet %in% d.bn$triplet,]
nrow(d.un)

d.bn$type = "tumor_specific"
d.bn$key2 <- paste(d.bn$var_id, d.bn$gene,sep="_")
d.un$key2 <- paste(d.un$var_id, d.un$gene,sep="_")




#Eat(" * Processing [",file, "]\n")
#x <- as.data.frame(data.table::fread("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/nominal/NORMAL.All.txt.gz",header=FALSE, stringsAsFactors = FALSE, sep=" "))
#x$key <- paste(x$V8, x$V1, sep="_")
#df <- x[x$key %in% d.bn$key2,]

#write.table(df, file="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/union/tspe_te_snp_association_gene.txt",sep="\t", row.names=F,col.names=TRUE,quote=F)




df <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/union/tspe_te_snp_association_gene.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t")
head(df)
df.n <- df[,c(12,14,16)]
head(df.n)
colnames(df.n) <- c("snp_gene_pval","snp_gene_slope","key2")
df.n$snp_gene_sig <- df.n$snp_gene_pval < 0.05
table(df.n$snp_gene_sig)

INDEPN <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/mapping/normal/NORMAL.conditional_All.txt.gz",header=TRUE,stringsAsFactors=FALSE,sep=" "))
INDEPN$key2 <- paste(INDEPN$var_id, INDEPN$phe_id, sep="_")

indep <- merge(d.bn, INDEPN, by="key2")

df$indep_level <- 0 
for(i in indep$key2)
{
  x <- indep[which(indep$key2 == i),]
  print(x)
  df$indep_level[which(df$key2 == i)] <- x$bwd_sig
}

df.un <- merge(d.un, df.n,by="key2")
head(df.un)


df.un$indep_level <- 0 
for(i in indep$key2)
{
  x <- indep[which(indep$key2 == i),]
  print(x)
  df.un$indep_level[which(df.un$key2 == i)] <- x$bwd_sig
}


table(df.un[which(df.un$indep_level == 0),]$normal_best)
table(df.un[which(df.un$indep_level == 1),]$normal_best)

# Figure 6C

pdf("inactive_TE_SNP_tumor_specific_TE_eqtls_INDEP_level.pdf",width=5.5,height=5)
barplot(table(df.un$indep_level),names.arg = c("Inactive eQTL\nin normal","Active eQTL\nin normal (5% FDR)"),ylab="Triplets with tumor-specific TE-eQTLs",main="Tumor-specific TE eQTL activity\n in normal for other TEs/genes")
bp = table(df.un$indep_level)
bp[1] <- bp[1]+ 2 # two eQTLs were not tested because exceeding the 1Mb window. I just add them as inactive. 
text(c(0.7,1.9),as.numeric(bp)/2, labels=as.numeric(bp))  
dev.off()



#### SHARED #####
shared <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/specific/raw/shared_eQTLs.txt",header=T,stringsAsFactors = F,sep=" ")
shared$key <- paste(shared$var_id, shared$phe_id, sep="_")

sh.bn <- merge(bn, shared, by="key")
sh.un <- UNION[UNION$triplet %in% sh.bn$triplet,]

sh.bn$type = "shared"

emod.t <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/emods/nominal/tumor/TUMOR.TE.GENE.nominal.chrAll.txt.gz",header=FALSE,stringsAsFactors=F,sep=" "))
emod.n <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/emods/nominal/normal/NORMAL.TE.GENE.nominal1.chrALL.txt.gz",header=FALSE,stringsAsFactors=F,sep=" "))

te_gene <- paste(d.bn$te,d.bn$gene, sep="_")
emod.t$key <- paste(emod.t$V8, emod.t$V1,sep="_")
emod.n$key <- paste(emod.n$V8, emod.n$V1, sep="_")

dn <- emod.n[emod.n$key %in% te_gene,]
dt <- emod.t[emod.t$key %in% te_gene,]

d.bn.table <- as.data.frame(table(d.bn$best),stringsAsFactors = F)
sh.bn.table <- as.data.frame(table(sh.bn$best),stringsAsFactors = F)
d.bn.table$type <- "Tumor-specific"
sh.bn.table$type <- "Shared"

col=brewer.pal(9,"Set1")
df <- rbind(d.bn.table,sh.bn.table)

#Figure 6A 

gg1<-ggplot(data=df, aes(x=Var1,y=Freq,fill=type))+
  geom_bar(stat="identity",color="black",size=1.15,position="dodge")+
  labs(y="Frequency",x="")+
  scale_x_discrete(labels=c("m1"="Causal","m2"="Reactive","m3"="Independent"))+
  scale_fill_manual(values=c("Shared"=col[2],"Tumor-specific"=col[1]),labels=c("Shared"="Shared TE eQTLs","Tumor-specific"="Tumor-specific TE eQTLs"))+
  theme_classic(base_size=20)+
  theme(legend.position=c(0.8,0.9),legend.title=element_blank())
ggsave(gg1,filename="tissue_specificity_TE_eQTLs_bestmodel_tumor.pdf",device="pdf",height=6.5,width=7)

d.un.table <- as.data.frame(table(d.un$switch),stringsAsFactors = F)
sh.un.table <- as.data.frame(table(sh.un$switch),stringsAsFactors = F)

d.un.table$type <- "Tumor-specific"
d.un.table$prop <- d.un.table$Freq / sum(d.un.table$Freq)
sh.un.table$type <- "Shared"
sh.un.table$prop <- sh.un.table$Freq / sum(sh.un.table$Freq)
df <- rbind(d.un.table,sh.un.table)

# Figure 6 B

gg2 <- ggplot(data=df, aes(x=Var1,y=prop,fill=type))+
  geom_bar(stat="identity",color="black",size=1.15,position="dodge")+
  labs(y="Proportion",x="")+
  scale_x_discrete(labels=c("no_switch"="no model change","switch"="model change","switch_causal"="model change\nto causal"))+
  scale_fill_manual(values=c("Shared"=col[2],"Tumor-specific"=col[1]),labels=c("Shared"=paste("Shared TE eQTLs\n(N=",sum(sh.bn.table$Freq),")",sep=""),"Tumor-specific"=paste("Tumor-specific TE eQTLs\n(N=",sum(d.bn.table$Freq),")",sep="")))+
  theme_classic(base_size=20)+
  theme(legend.position="top",legend.title=element_blank())
ggsave(gg2,filename="tissue_specificity_TE_eQTLs_model_switch.pdf",device="pdf",height=6.5,width=7)

#### Enrichment analysis ####

mat <- matrix(c(77,140,51,223),byrow=TRUE,ncol=2)
fisher.test(mat)
#data:  mat
#p-value = 3.147e-05
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  1.560965 3.715096
#sample estimates:
#  odds ratio 
#2.400463



#Fisher's Exact Test for Count Data
mat <- matrix(c(80,136,53,184),byrow=TRUE,ncol=2)
fisher.test(mat)
mat
# New results --> Tumor specific switching to causal vs shared switch to causal. 
#data:  mat
#p-value = 0.0006588
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 1.326012 3.153268
#sample estimates:
#odds ratio 
#  2.038923