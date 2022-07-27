#!/usr/bin/env Rscript 

library(RColorBrewer)
library(ggplot2)
library(dplyr)



##========================================================##
##========= CHECKING TUMOR SPECIFIC AND TRIPLETS =========##
##========================================================##


UNION <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/teqtls/shared/UNION_normal_tumor_with_posteriors.txt",header=T,stringsAsFactors=F,sep="\t"))
UNION$switch <- "switch"
UNION[which(UNION$normal_best == UNION$tumor_best),]$switch <- "no_switch"
UNION[which(UNION$normal_best != UNION$tumor_best & UNION$tumor_best == "m1"),]$switch <- "switch_causal"


#UNION$sub <- paste(UNION$normal_best, UNION$tumor_best,sep="_")


bn <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/teqtls/bn/tumor/TUMOR.BNlearn.ALL.extraINFO.txt.gz",header=TRUE,stringsAsFactors=F,sep="\t"))
tspe <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/specificity/TUMOR.specific_TE_eQTLs.txt.gz",header=T,stringsAsFactors = F,sep=" ")
head(tspe)
bn$triplet <- paste(bn$var_id, bn$gene, bn$te, sep="_")
bn$key <- paste(bn$var_id, bn$te,sep="_")
tspe$key <- paste(tspe$var_id, tspe$phe_id,sep="_")
#UNION$key <- paste(UNION$var_id, SHARED$te,sep="_")
UNION$triplet <- paste(UNION$var_id, UNION$gene, UNION$te, sep="_")
d.bn <- merge(bn,tspe, by="key")
head(d.bn)

d.un <- UNION[UNION$triplet %in% d.bn$triplet,]
nrow(d.un)


d.bn$type = "tumor_specific"
d.bn$key2 <- paste(d.bn$var_id.x, d.bn$gene,sep="_")
d.un$key2 <- paste(d.un$var_id, d.un$gene,sep="_")
##### CHECKING WHETHER SNP IS ALSO EQTL FOR GENE IN NORMAL OR NOT ##### 
# THE ABOVE WILL LET US KNOW WHETHER ITS THE SNP THAT GETS ACTIVATED DURING TUMORIGENESIS OR THE TE.


nom_files <- Sys.glob("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/nominal/normal/tmp/*.txt.gz")

df <- data.frame()
for(file in nom_files)
{
  cat(" * Processing [",file, "]\n")
  x <- as.data.frame(data.table::fread(paste0(file),header=FALSE, stringsAsFactors = FALSE, sep=" "))
  x$key <- paste(x$V8, x$V1, sep="_")
  
  y <- x[x$key %in% d.bn$key2,]
  if(nrow(y) == 0){next}
  df <- rbind(df, y)
  print(df)
}
write.table(df, file="tspe_te_snp_association_gene.txt",sep="\t", row.names=F,col.names=TRUE,quote=F)

df <- read.table("",header=TRUE,stringsAsFactors = FALSE,sep="\t")
head(df)
df.n <- df[,c(12,14,16)]
head(df.n)
colnames(df.n) <- c("snp_gene_pval","snp_gene_slope","key2")
df.n$snp_gene_sig <- df.n$snp_gene_pval < 0.05
table(df.n$snp_gene_sig)

INDEPN <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/mapping/normal/NORMAL.conditional_All.txt.gz",header=TRUE,stringsAsFactors=FALSE,sep=" "))
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



pdf("inactive_TE_SNP_tumor_specific_TE_eqtls_INDEP_level.pdf",width=5.5,height=5)
barplot(table(df.un$indep_level),names.arg = c("Inactive eQTL\nin normal","Active eQTL\nin normal (5% FDR)"),ylab="Triplets with tumor-specific TE-eQTLs",main="Tumor-specific TE eQTL activity\n in normal for other TEs/genes")
bp = table(df.un$indep_level)
bp[1] <- bp[1]+ 2 # two eQTLs were not tested because exceeding the 1Mb window. I just add them as inactive. 
text(c(0.7,1.9),as.numeric(bp)/2, labels=as.numeric(bp))  
dev.off()



a <- df.un[which(df.un$indep_level == 0),]
a.tot <- nrow(a)
a.causal <- sum(a$M1_n) / a.tot
a.reactive <- sum(a$M2_n) / a.tot 
a.indep <- sum(a$M3_n) / a.tot

tot <- data.frame("model"=c("m1","m2","m3"),"mean_probability"=c(a.causal,a.reactive,a.indep))


b <- df.un[which(df.un$indep_level == 1),]
b.tot <- nrow(b)
b.causal <- sum(b$M1_n) / b.tot
b.reactive <- sum(b$M2_n) / b.tot 
b.indep <- sum(b$M3_n) / b.tot

totb <- data.frame("model"=c("m1","m2","m3"),"mean_probability"=c(b.causal,b.reactive,b.indep))


tot$type <- "Inactive SNP"
totb$type <- "Inactive TE"

all <- rbind(tot,totb)
head(all)

COL <- brewer.pal(9,"Set3")
gg <- ggplot(data=all, aes(x=type,y=mean_probability,fill=model,label=round(mean_probability*100,0)))+
  geom_bar(stat="identity",size=1.25,color="black")+
  geom_text(position=position_stack(vjust=0.5),size=7)+
  scale_fill_manual(values=c(COL[3],COL[4],COL[6]),labels=c("m1" = "Causal (V->T->G)","m2"="Reactive (V->G->T)","m3"="Independent (V->T & V->G)"))+
  labs(y="Mean probability",x="")+
  theme_classic(base_size=20)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="top",
        legend.direction="vertical",
        legend.justification="left",
        legend.title=element_blank())



ggsave(gg, filename="model_probabilities_te_snp_activity_tumor_specific_INDEP_level.pdf",device="pdf",height=7,width=7.5)


BED <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_expression/normal/raw/CPM_trimmed_Genes_proteinCoding_lincRNA_TEs_filtered_50_percent_and_more_275_normal.raw.bed.gz",header=TRUE, stringsAsFactors = FALSE,sep="\t"))

a_te <- unique(a$te)

aa <- BED[BED$gene %in% a_te,]
aa <- aa[,-c(1:6)]
aaa <- as.numeric(apply(aa, 1, median))

b_te <- unique(b$te)
bb <- BED[BED$gene %in% b_te,]
bb <- bb[,-c(1:6)]
bbb <- as.numeric(apply(bb, 1, median))

d <- data.frame("case"=c(replicate(length(aaa),"inactive_snp"),replicate(length(bbb),"active_snp")),"median_exp"=c(aaa,bbb))
head(d)
ggplot(data=d, aes(x=case, y=median_exp))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_classic(base_size=20)



#### SHARED #####


pdf("tumor_specific_TE_eqtls_model.pdf",8,8)
par(mfrow=c(2,2))
barplot(table(d.bn$best),main="Tumor-specific TE-eQTLs model",names.arg = c("Causal","Reactive","Independent"))
barplot(table(d.un$sub),main="Tumor-specific TE-eQTLs substitutions",las=2)
barplot(table(d.un$switch),main="Tumor-specific TE-eQTLs model switch")
plot(tspe[tspe$key %in% d.bn$key,]$nominal_slope, tspe[tspe$key %in% d.bn$key,]$bwd_slope,xlab="Effect size in normal",ylab="effect size in tumor")
dev.off()

shared <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/specificity/SHARED.TE_eQTLs.txt.gz",header=T,stringsAsFactors = F,sep=" ")
shared$key <- paste(shared$var_id, shared$phe_id, sep="_")

sh.bn <- merge(bn, shared, by="key")
sh.un <- UNION[UNION$triplet %in% sh.bn$triplet,]

sh.bn$type = "shared"

pdf("Shared_specific_TE_eQTLs_models.pdf")
par(mfrow=c(2,2))
barplot(table(sh.bn$best),main="Shared TE-eQTLs model",names.arg=c("Causal","Reactive","Independent"))
barplot(table(sh.un$sub),main="Shared TE-eQTLs model substitution",las=2)
barplot(table(sh.un$switch),main="Tumor-specific TE-eQTLs model switch")
plot(shared[shared$key %in% sh.bn$key,]$n_eqtl_bwd_slope, shared[shared$key %in% sh.bn$key,]$t_eqtl_bwd_slope,xlab="Effect size in normal",ylab="Effect size in tumor")
dev.off()

emod.t <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/emods/nominal/THRESHOLD/tumor/TUMOR.TE.GENE.nominal_All.txt.gz",header=TRUE,stringsAsFactors=F,sep=" "))
emod.n <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/emods/nominal/FULL/normal/NORMAL.TE.GENE.nominal1_All.txt.gz",header=TRUE,stringsAsFactors=F,sep=" "))

te_gene <- paste(d.bn$te,d.bn$gene, sep="_")
emod.t$key <- paste(emod.t$var_id, emod.t$phe_id,sep="_")
emod.n$key <- paste(emod.n$var_id, emod.n$phe_id, sep="_")

dn <- emod.n[emod.n$key %in% te_gene,]
dt <- emod.t[emod.t$key %in% te_gene,]

pdf("TE_gene_effect_size-normal_tumor_for_tumor_specific_eQTLs.pdf",8,4.5)
par(mfrow=c(1,2))
plot(dn$slope, dt$slope,main="TE gene association effect size\ncomparison in normal and tumor\nfor the tumor-specific TE-eQTLs",xlab="Effect size (regression slope) in normal",ylab="Effect size (regression slope) in tumor",ylim=c(0,1),xlim=c(min(dn$slope),1))
abline(v = 0,lty=2)
abline(h=0,lty=2)


d.bn.table <- as.data.frame(table(d.bn$best),stringsAsFactors = F)
sh.bn.table <- as.data.frame(table(sh.bn$best),stringsAsFactors = F)
d.bn.table$type <- "Tumor-specific"
sh.bn.table$type <- "Shared"

col=brewer.pal(9,"Set1")
df <- rbind(d.bn.table,sh.bn.table)
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

gg2 <- ggplot(data=df, aes(x=Var1,y=prop,fill=type))+
  geom_bar(stat="identity",color="black",size=1.15,position="dodge")+
  labs(y="Proportion",x="")+
  scale_x_discrete(labels=c("no_switch"="no model change","switch"="model change","switch_causal"="model change\nto causal"))+
  scale_fill_manual(values=c("Shared"=col[2],"Tumor-specific"=col[1]),labels=c("Shared"=paste("Shared TE eQTLs\n(N=",sum(sh.bn.table$Freq),")",sep=""),"Tumor-specific"=paste("Tumor-specific TE eQTLs\n(N=",sum(d.bn.table$Freq),")",sep="")))+
  theme_classic(base_size=20)+
  theme(legend.position="top",legend.title=element_blank())
ggsave(gg2,filename="tissue_specificity_TE_eQTLs_model_switch.pdf",device="pdf",height=6.5,width=7)

d.bn.slope <- d.bn[,c(41,42)]
sh.bn.slope <- sh.bn[,c(41,43)]

colnames(d.bn.slope) <- c("normal_slope","tumor_slope")
colnames(sh.bn.slope) <- c("normal_slope","tumor_slope")
d.bn.slope$type <- "tumor_specific"
sh.bn.slope$type <- "shared"

d <- rbind(d.bn.slope,sh.bn.slope)
gg3 <- ggplot(data=d,aes(x=normal_slope,y=tumor_slope,colour=type))+
  geom_point()+
  labs(x="Effect size (regression slope)\nin normal",y="Effect size (regression slope)\nin tumor")+
  scale_colour_manual(values=c("shared"=col[2],"tumor_specific"=col[1]),labels=c("shared"=paste("Shared TE eQTLs\n(N=",sum(sh.bn.table$Freq),")",sep=""),"tumor_specific"=paste("Tumor-specific TE eQTLs\n(N=",sum(d.bn.table$Freq),")",sep="")))+
  theme_classic(base_size=20)+
  theme(legend.position="top",legend.title=element_blank())
ggsave(gg3,filename="tissue_specificity_TE_eQTLs_model_effectSize.pdf",device="pdf",height=6,width=6.5)
gg2













tes <- unique(d.bn$triplet)
pos <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_expression/stat/features_genomic_positions.bed.gz",header=TRUE,stringsAsFactors = F))

TEs <- union(tes,te_cgc)
head(pos)
pos <- pos[pos$id %in% TEs,]
head(pos)
nrow(pos)

write.table(pos, file="check_TcTs.bed",col.names=TRUE,row.names=F,quote=F)






