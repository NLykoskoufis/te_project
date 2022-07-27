#!/bin/usr/env Rscript 

library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggthemes)
library(ggsci)

setwd("~/Documents/PROJECTS/te_project/paper/rerun/Figures/Figure1")

#Figure 1A 
ere <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/info_files/ereinfoTable.txt.gz", header=TRUE,stringsAsFactors = FALSE,sep="\t"))
ere$teID <- paste(ere$unMerged, ere$start, ere$end, ere$strand, sep="_")

BED <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/normal/SYSCOL_normal_275_gene_TE.FM05.resid.cpm.bed.gz", header=TRUE,stringsAsFactors = FALSE,sep="\t"))
BED <- BED[,c(1:6)]

tmp <- merge(ere[,c(2,3,4,7,8,9)],BED,by.x="teID",by.y="id")
tmp <- tmp[,c(2,3,4,1,5,6)]
colnames(tmp) <- c("#chr","start","end","id","fam","name")
write.table(tmp, file="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/stat/features_family.bed.gz",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

tab <- as.data.frame(prop.table(table(tmp$fam)))
tab$Var1 <- as.character(tab$Var1)
tab$family <- sapply(strsplit(tab$Var1, "\\/"),"[[",1)
tab <- tab[!grepl("\\?",tab$family),]

gg1 <- ggplot(data=tab, aes(x=Var1, y=Freq,fill=family))+
  geom_bar(stat="identity",color="black")+
  labs(x="TE subfamilies",y="Proportion",title="Expressed TE subfamilies", fill="TE family")+
  scale_fill_npg()+
  theme_base(base_size=12)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.title=element_text(face="bold"),
        plot.title=element_text(hjust=0.5), axis.title = element_text(face="bold"),
        plot.background=element_blank())
ggsave(gg1, filename="figure1A.pdf",height=6,width=8,device="pdf")

#Figure 1A.2 ALTERNATIVE 
## Instead of taking the proportion of TEs subfamilies expressed in the whole do it per subfamily. 

ALL <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_TEannotation/hg19_TE_repmask_LTRm_s_20140131.maptable.bed.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
ALL$id <- paste(ALL$id, ALL$start+1, ALL$end, ALL$strand, sep="_")

tmp <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/stat/features_family.bed.gz", header=TRUE,stringsAsFactors = FALSE,sep="\t"))

allFAM <- as.data.frame(table(ALL$gid), stringsAsFactors=FALSE)
expFAM <- as.data.frame(table(tmp$fam), stringsAsFactors = FALSE)

colnames(allFAM) <- c("Var1","allFreq")
colnames(expFAM) <- c("Var1","expFreq")

df <- merge(allFAM, expFAM, by="Var1")
df$prop <- df$expFreq / df$allFreq
df$Var1 <- as.character(df$Var1)
df$family <- sapply(strsplit(df$Var1, "\\/"),"[[",1)
df <- df[!grepl("\\?", df$family),]


gg2 <- ggplot(data=df, aes(x=Var1, y=prop,fill=family))+
  geom_bar(stat="identity",color="black")+
  labs(x="TE subfamilies",y="Proportion",title="Proportion of expressed TEs per subfamily", fill="TE family")+
  scale_fill_npg()+
  theme_base(base_size=12)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.title=element_text(face="bold"),
        plot.title=element_text(hjust=0.5), axis.title = element_text(face="bold"),
        plot.background=element_blank())
ggsave(gg2, file="figure1A2.pdf",height=6,width=8,device="pdf")

colnames(tab) <- c("Var1","propvsall","family")

df <- merge(df, tab, by=c("Var1","family"))

df.melt <- reshape2::melt(df[,c(1,5,6)], id.vars="Var1")
df.melt$family <- sapply(strsplit(df.melt$Var1, "\\/"),"[[",1)
df.melt$variable <- as.character(df.melt$variable)
df.melt$variable[which(df.melt$variable == "prop")] <- "Proportion of expressed\nTEs per subfamily"
df.melt$variable[which(df.melt$variable == "propvsall")] <- "Proportion of expressed TE\nsubfamilies compared to all \nexpressed TEs"

ggALT <- ggplot(data=df.melt, aes(x=Var1, y=value, fill=family))+
  geom_bar(stat="identity",color="black")+
  labs(x="TE subfamilies", y="Proportion")+
  facet_grid(~ variable, scales="free")+
  coord_flip()+
  scale_fill_npg()+
  theme_base()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.title=element_text(face="bold"),
        plot.title=element_text(hjust=0.5), axis.title = element_text(face="bold"),
        plot.background=element_blank(),
        strip.text.x = element_text(face="bold"))

ggsave(ggALT, file="figure1ALT.pdf",height=9,width=9, device="pdf")



#Figure 1B
tmp <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/stat/TE_ensemblRegulatoryBuild_miRBase_overlap.bed.gz",header = FALSE,sep="\t",stringsAsFactors = FALSE))
length(unique((tmp[which(tmp$V10 != "."),]$V4)))

p <- table(tmp$V10)
p <- p[-1]
p <- p/ sum(p)


lbls <- paste(names(p), "\n" , round(as.numeric(p)*100,2),"%",sep="")
pdf("regulatory_element_overlap_TE.pdf",7,6.5)
pie(p,col=RColorBrewer::brewer.pal(8, "Dark2"),labels=lbls, main="TEs with overlaping regulatory elements\n(N=13,590)")
dev.off()




### CHECKING ENRICHMENT OF EXPRESSED VS NON EXPRESSED TE FOR REG ELEMENTS ### 

ALL <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_TEannotation/TE_ensemblRegBuild_miRBase_overlap_ALL_TEs.bed.gz",header=FALSE,stringsAsFactors = FALSE,sep="\t"))
colnames(ALL) <- c("phe_chr","phe_from","phe_to", "phe_id","phe_gid","phe_strd","reg_chr","reg_from","reg_to","reg_id","reg_gid","reg_strd")
ALL$phe_id <- paste(ALL$phe_id,ALL$phe_from+1, ALL$phe_to, ALL$phe_strd, sep="_")

exp <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/stat/TE_ensemblRegulatoryBuild_miRBase_overlap.bed.gz",header = FALSE,sep="\t",stringsAsFactors = FALSE))
colnames(exp) <- c("phe_chr","phe_from","phe_to", "phe_id","phe_gid","phe_strd","reg_chr","reg_from","reg_to","reg_id","reg_gid","reg_strd")
exp <- exp[,-ncol(exp)]

REG_EL <- unique(ALL$reg_id)

dataFrame <- data.frame()


for (rrr in REG_EL)
{
  NON_EXP <- ALL[! ALL$phe_id %in% exp$phe_id,]
  
  ###       | EXP  | NO EXP | 
  ### OV    | A    |   B    | 
  ### NO OV | C    |   D    | 
  
  A <- length(unique(exp[which(exp$reg_id == rrr),]$phe_id))
  B <- length(unique(NON_EXP[which(NON_EXP$reg_id == rrr),]$phe_id))
  C <- length(unique(exp[which(exp$reg_id != rrr),]$phe_id))
  D <- length(unique(NON_EXP[which(NON_EXP$reg_id != rrr),]$phe_id))
  mat <- matrix(c(A,B,C,D), ncol=2, byrow=TRUE)
  ft <- fisher.test(mat)
  print(ft)
  
  dataFrame <- rbind(dataFrame, data.frame("reg_id"=rrr,"A"=A, "B"=B, "C"=C, "D"=D, "pvalue"=as.numeric(ft$p.value), "estimate"=as.numeric(ft$estimate), "confInt1"=as.numeric(ft$conf.int[1]), "confInt2"=as.numeric(ft$conf.int[2])))
}
dataFrame <- dataFrame[-2,]
dataFrame$padj <- p.adjust(dataFrame$pvalue, method = "fdr")

write.table(dataFrame, file="ExpressedTEsVSnonExpressed_enrichment_RegEl.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

ggEnrichment <- ggplot(data=dataFrame, aes(x=reg_id, y=estimate,fill=padj))+
  geom_point(size=3,shape=21, colour="black")+
  labs(x="Regulatory elements", y = "Odds-ratio",fill="adjusted\np-value")+
  coord_flip()+
  theme_bw(base_size=12)
ggsave(ggEnrichment, filename="reg_el_TE_enrichment.pdf",height=4,width=5)



