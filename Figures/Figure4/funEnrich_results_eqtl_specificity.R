#!/usr/bin/env Rscript 

setwd("/Users/nikolaoslykoskoufis/Documents/PROJECTS/te_project/paper/rerun/Figures/Figure4")
##LIBRARIES##
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)


### DGE FILES ### 
tf <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/analysis_DGE/DGE_results_DESEq2_all_TE_genes.txt",header=T,stringsAsFactors = F,sep="\t"))
tf$expRatio <- tf$t_median / tf$n_median

convert <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_expression/stat/ensemblToGeneName.txt",header=F,sep=" ",stringsAsFactors = F)
colnames(convert) <- c("id","geneName")
tf <- merge(convert,tf,by="id")

#### changing names of tf genes to match ddte. 
tf[which(tf$id == "ENSG00000085276"),]$geneName <- "EVI1"
tf[which(tf$id == "ENSG00000125820"),]$geneName <- "NKX2.2"
tf[which(tf$id == "ENSG00000157557"),]$geneName <- "Ets-2"
tf <- tf[which(tf$padj <= 0.05),]

### USEFUL FUNCTIONS #### 
combine_results <- function(x)
{
  #### x -> list of files obtained with fenrich to combine together
  
  df <- data.frame()
  for(file in x)
  {
    r <- as.data.frame(data.table::fread(paste0(file), header=TRUE,stringsAsFactors = FALSE, sep="\t"))
    tmp <- gsub("_fenrich_.*","",basename(file))
    if(length(strsplit(tmp,"\\_")[[1]]) == 2)
    {
      mark <- sapply(strsplit(basename(tmp),"\\_"), "[[",2)
      cellType <- sapply(strsplit(basename(tmp),"\\_"), "[[",1)
    }else{
      mark <- sapply(strsplit(basename(tmp),"\\_"), "[[",3)
      cellType <- paste(sapply(strsplit(basename(tmp),"\\_"), "[[",1),sapply(strsplit(basename(tmp),"\\_"), "[[",2),sep="_")
    }
    
    df <- rbind(df, cbind(data.frame("mark"=mark,"celltype"=cellType), r))
  }
  return(df)
}

filterEnrichment <- function(t)
{
  
  t.sig.both <- t[which(t$sig.x == TRUE & t$sig.y == TRUE),]
  t.sig.x.not.y <- t[which(t$sig.x == TRUE & t$sig.y == FALSE),]
  t.sig.y.not.x <- t[which(t$sig.y == TRUE & t$sig.x == FALSE),]
  
  pb1 <- t.sig.y.not.x[t.sig.y.not.x$OddsRatio.x > t.sig.y.not.x$OddsRatio.y,]
  pb2 <- t.sig.x.not.y[t.sig.x.not.y$OddsRatio.y > t.sig.x.not.y$OddsRatio.x,]
  print(c(pb1$mark,pb2$mark))
  t.sig.y.not.x <- t.sig.y.not.x[!t.sig.y.not.x$OddsRatio.x > t.sig.y.not.x$OddsRatio.y,]
  t.sig.x.not.y <- t.sig.x.not.y[!t.sig.x.not.y$OddsRatio.y > t.sig.x.not.y$OddsRatio.x,]
  toKeep <- unique(c(t.sig.both$mark, t.sig.x.not.y$mark, t.sig.y.not.x$mark))
  return(toKeep)
}

tspe_files_te <- Sys.glob("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/specific/tumor_specific/*TE.results_k100_w2500.txt")
tspe_files_gene <- Sys.glob("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/specific/tumor_specific/*GENE.results_k100_w2500.txt")

shared_files_te <- Sys.glob("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/specific/shared/*TE.results_k100_w2500.txt")
shared_files_gene <- Sys.glob("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/specific/shared/*GENE.results_k100_w2500.txt")



TSPE_TE <- combine_results(tspe_files_te)
SHARED_TE <- combine_results(shared_files_te)
TSPE_GENE <- combine_results(tspe_files_gene)
SHARED_GENE <- combine_results(shared_files_gene)



performEnrichment <- function(df){
  data <- data.frame()
  for(rrr in c(1:nrow(df)))
  {
    
    x <- df[rrr,]
    if (x$A == 0 & x$B == 0){
      next;
    }else{
      
      if((x$A != 0 & x$B == 0) | (x$A == 0 & x$B != 0))
      {
        x$A = x$A + 1
        x$B = x$B + 1
        mat <- matrix(c(x$A,x$B,x$C,x$D),byrow=TRUE,ncol=2)
        f = fisher.test(mat)
        or <- as.numeric(f$estimate)
        pvalue <- as.numeric(f$p.value)
        data <- rbind(data, cbind(x,"OddsRatio"=or,"pvalue"=pvalue))
      }else{
        mat <- matrix(c(x$A,x$B,x$C,x$D),byrow=TRUE,ncol=2)
        f = fisher.test(mat)
        or <- as.numeric(f$estimate)
        pvalue <- as.numeric(f$p.value)
        data <- rbind(data, cbind(x,"OddsRatio"=or,"pvalue"=pvalue))
      }
    }
  }
  return(data)
  
}
tspe_te <- performEnrichment(TSPE_TE)
shared_te <- performEnrichment(SHARED_TE)
tspe_gene <- performEnrichment(TSPE_GENE)
shared_gene <- performEnrichment(SHARED_GENE)




# FDR correction
tspe_te$fdr <- p.adjust(tspe_te$pvalue, method="fdr")
shared_te$fdr <- p.adjust(shared_te$pvalue, method="fdr")
tspe_gene$fdr <- p.adjust(tspe_gene$pvalue, method="fdr")
shared_gene$fdr <- p.adjust(shared_gene$pvalue, method="fdr")


## Histogram of pvalue distribution ## 
pdf("pvalue_distribution_Fenrich_eqtl_specificity.pdf",height=4, width=4.5)
hist(tspe_te$pvalue, main="Functional enrichment tumor-specific\nTE-eQTLspvalue distribution", xlab="P-value")
hist(tspe_gene$pvalue, main="Functional enrichment tumor-specific\ngene-eQTLs pvalue distribution", xlab="P-value")
hist(shared_te$pvalue, main="Functional enrichment shared\nTE-eQTLs pvalue distribution", xlab="P-value")
hist(shared_gene$pvalue, main="Functional enrichment shared\ngene-eQTLs pvalue distribution", xlab="P-value")
dev.off()

tspe_te$mark <- paste(tspe_te$celltype, tspe_te$mark, sep=" ") 
shared_te$mark <- paste(shared_te$celltype, shared_te$mark, sep=" ")
tspe_gene$mark <- paste(tspe_gene$celltype, tspe_gene$mark, sep=" ") 
shared_gene$mark <- paste(shared_gene$celltype, shared_gene$mark, sep=" ")


write.table(tspe_te, file="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/specific/tumor_specific/tumor_specific_TEeQTL_enrichment_results_maf002_w2500_k100.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(shared_te, file="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/specific/shared/shared_TEeQTL_enrichment_results_maf002_w2500_k100.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(tspe_gene, file="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/specific/tumor_specific/tumor_specific_GENEeQTL_enrichment_results_maf002_w2500_k100.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(shared_gene, file="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/specific/shared/shared_GENEeQTL_enrichment_results_maf002_w2500_k100.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


tspe_te <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/specific/tumor_specific/tumor_specific_TEeQTL_enrichment_results_maf002_w2500_k100.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
shared_te <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/specific/shared/shared_TEeQTL_enrichment_results_maf002_w2500_k100.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
tspe_gene <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/specific/tumor_specific/tumor_specific_GENEeQTL_enrichment_results_maf002_w2500_k100.txt",sep="\t",header=TRUE, stringsAsFactors = FALSE)
shared_gene <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/specific/shared/shared_GENEeQTL_enrichment_results_maf002_w2500_k100.txt",sep="\t",header=TRUE, stringsAsFactors = FALSE)

#************************************#
#    Functional enrichment for TEs   # 
#************************************#

# Individual plots
tspe_te$sig <- tspe_te$fdr <= 0.05
tspe_te$sig10 <- tspe_te$fdr <= 0.1

tspe_te_gg <- ggplot(data=tspe_te[tspe_te$sig,], aes(x=reorder(mark,-OddsRatio),y=OddsRatio,label=round(-log10(fdr),1)))+
  geom_bar(stat="identity")+
  geom_hline(yintercept=1, lty=2)+
  theme_classic(base_size=20)+
  scale_fill_brewer(palette="Set2")+
  labs(x="",y="Enrichment over the null",title="Functional enrichment tumor-specific TE-eQTLs")+
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1,vjust=0.5),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.85,0.9),legend.title=element_blank(),legend.background = element_rect(color="black",size=0.5),legend.margin=margin(t = 0,b=0.1,r=0.1,l=0.1, unit='cm'))
ggsave(tspe_te_gg, filename="functional_enrichment_tspe_te_LoVo_fdr5_k100_w2500.pdf",device="pdf",height=7,width=15)


shared_te$sig <- shared_te$fdr <= 0.05
shared_te$sig10 <- shared_te$fdr <= 0.1
shared_te_gg <- ggplot(data=shared_te[shared_te$sig,], aes(x=reorder(mark,-OddsRatio),y=OddsRatio,label=round(-log10(fdr),1)))+
  geom_bar(stat="identity")+
  geom_hline(yintercept=1, lty=2)+
  theme_classic(base_size=20)+
  scale_fill_brewer(palette="Set2")+
  labs(x="",y="Enrichment over the null",title="Functional enrichment shared TE-eQTLs")+
  theme(axis.text.x = element_text(size=9,angle=90,hjust=1,vjust=0.5),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.85,0.9),legend.title=element_blank(),legend.background = element_rect(color="black",size=0.5),legend.margin=margin(t = 0,b=0.1,r=0.1,l=0.1, unit='cm'))
ggsave(shared_te_gg, filename="functional_enrichment_shared_te_LoVo_fdr5_k100.pdf",device="pdf",height=7,width=15)


# combined tumor-spe and shared plot
comb_te <- merge(tspe_te, shared_te, by="mark")
t <- comb_te 

toKeep_te <- filterEnrichment(t)
toKeep <- t[which(t$sig.x | t$sig.y),]$mark




# if tumor specific is not significant but has higher OR than normal and exclude. If shared is not significant and has higher OR than tumor-specific, exclude. 

TMP_tspe_te <- tspe_te[tspe_te$mark %in% toKeep_te,c(1,7,8)]
TMP_tspe_te$type <- "tumor-specific"

TMP_shared_te <- shared_te[shared_te$mark %in% toKeep_te, c(1,7,8)]
TMP_shared_te$type <- "shared"

TMP_te <- rbind(TMP_shared_te,TMP_tspe_te)
TMP_te$OddsRatio[which(TMP_te$OddsRatio > 30)] <- 30
TMP_te$ypos <- 0
for(i in unique(TMP_te$mark)){
  x <- TMP_te[which(TMP_te$mark == i),]
  mm <- max(x$OddsRatio)
  TMP_te[which(TMP_te$mark == i & TMP_te$type == "tumor-specific"),]$ypos <- mm + 1.5
  TMP_te[which(TMP_te$mark == i & TMP_te$type == "shared"),]$ypos <- mm + 4
}


subset_to_order <- subset(TMP_te,type == "tumor-specific")
subset_to_order$mark <- with(subset_to_order, reorder(mark,-OddsRatio))
TMP_te$mark = factor(TMP_te$mark, levels = levels(subset_to_order$mark))


gg_combined_te <- ggplot(data=TMP_te, aes(x=mark, y=OddsRatio, fill=type,label=round(-log10(pvalue),1)))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_hline(yintercept=1, lty=2, colour="grey")+
  geom_text(aes(x=mark, y=ypos),angle=90,size=3)+
  theme_classic(base_size=20)+
  labs(x="",y="Enrichment over the null",title="Functional enrichment for tumor-specific and shared TE-eQTLs")+
  scale_fill_manual(values=c("tumor-specific"="#e41a1c", "shared"="black"))+
  theme(axis.text.x = element_text(size=9,angle=90,hjust=1,vjust=0.5),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.9,0.9), legend.title=element_blank())
ggsave(gg_combined_te, filename="functional_enrichment_combined_tspe_shared_TE_eQTLs_k100_maf002_w2500.pdf",device="pdf",height=7,width=20)

### EXCLUDED ###

toExclude <- comb_te[!comb_te$mark %in% toKeep,]$mark

exclude_tspe_te <- tspe_te[tspe_te$mark %in% toExclude,c(1,7,8)]
exclude_tspe_te$type <- "tumor-specific"

exclude_shared_te <- shared_te[shared_te$mark %in% toExclude, c(1,7,8)]
exclude_shared_te$type <- "shared"

exclude_te <- rbind(exclude_shared_te,exclude_tspe_te)

exclude_te$ypos <- 0
for(i in unique(exclude_te$mark)){
  x <- exclude_te[which(exclude_te$mark == i),]
  mm <- max(x$OddsRatio)
  exclude_te[which(exclude_te$mark == i & exclude_te$type == "tumor-specific"),]$ypos <- mm + 1.5
  exclude_te[which(exclude_te$mark == i & exclude_te$type == "shared"),]$ypos <- mm + 4
}


subset_to_order <- subset(exclude_te,type == "tumor-specific")
subset_to_order$mark <- with(subset_to_order, reorder(mark,-OddsRatio))
exclude_te$mark = factor(exclude_te$mark, levels = levels(subset_to_order$mark))


gg_excluded_te <- ggplot(data=exclude_te, aes(x=mark, y=OddsRatio, fill=type,label=round(-log10(pvalue),1)))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_hline(yintercept=1, lty=2, colour="grey")+
  geom_text(aes(x=mark, y=ypos),angle=90,size=3)+
  theme_classic(base_size=20)+
  labs(x="",y="Enrichment over the null",title="Functional enrichment for tumor-specific and shared TE-eQTLs\nEXCLUDED CASES")+
  scale_fill_manual(values=c("tumor-specific"="#e41a1c", "shared"="black"))+
  theme(axis.text.x = element_text(size=9,angle=90,hjust=1,vjust=0.5),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.9,0.9), legend.title=element_blank())
ggsave(gg_excluded_te, filename="functional_enrichment_excluded_cases_combined_tspe_shared_TE_eQTLs.pdf",device="pdf",height=7,width=20)

gg_pb_te <- ggplot(data=exclude_te[exclude_te$mark %in% pb$mark, ], aes(x=mark, y=OddsRatio, fill=type,label=round(-log10(pvalue),1)))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_hline(yintercept=1, lty=2, colour="grey")+
  geom_text(aes(x=mark, y=ypos),angle=90,size=3)+
  theme_classic(base_size=20)+
  labs(x="",y="Enrichment over the null",title="Functional enrichment for tumor-specific and shared TE-eQTLs\nPROBLEMATIC CASES")+
  scale_fill_manual(values=c("tumor-specific"="#e41a1c", "shared"="black"))+
  theme(axis.text.x = element_text(size=9,angle=90,hjust=1,vjust=0.5),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.9,0.9), legend.title=element_blank())
ggsave(gg_pb_te, filename="functional_enrichment_problematic_cases_combined_tspe_shared_TE_eQTLs.pdf",device="pdf",height=7,width=20)

### END OF EXCLUDED ###


# comparison tumor-specific and shared enrichment


comb_te$ratio <- log2(comb_te$OddsRatio.x / comb_te$OddsRatio.y)
gg_comb_te1 <- ggplot(data=comb_te[comb_te$mark %in% toKeep_te,], aes(x=reorder(mark,-ratio), y=ratio))+
  geom_bar(stat="identity",color="black",fill="grey")+
  geom_hline(yintercept=0,lty=2,size=0.2)+
  theme_classic(base_size=20)+
  labs(x="",y="Log2(tumor-specific enrichment/shared enrichment)")+
  scale_fill_brewer(palette="Set2",direction = -1)+
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1,vjust=0.5),
        axis.title.y=element_text(size=15),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.85,0.9),legend.title=element_blank(),legend.background = element_rect(color="black",size=0.5),legend.margin=margin(t = 0,b=0.1,r=0.1,l=0.1, unit='cm'))
gg_comb_te1


ggsave(gg_comb_te1, filename="functional_enrichment_ratio_tspe_shared_TE_eQTLs_tumor_spe_FDR5_only_k100_maf002_w2500.pdf",device="pdf",height=7.5,width=20)  


### ADD GENE FOLD CHANGE ###

comb_te$geneName <- sapply(strsplit(comb_te$mark, "\\ "),"[[",2)
comb_te <- comb_te[comb_te$mark %in% toKeep_te,]
df <- merge(comb_te, tf, by="geneName")

sig <- df[which(df$ratio >0),]
nrow(sig[which(log2(sig$expRatio) > 0),])

COL=RColorBrewer::brewer.pal(9,"Set1")
gg_combined_tf_te <- ggplot(data=df[df$mark %in% toKeep_te,], aes(x=reorder(mark, -ratio), y=ratio, fill="Tumor-specific enrichment/shared enrichment"))+
  geom_bar(stat="identity",color="grey",width=0.8)+
  geom_line(aes(y=log2(expRatio), group="mark",colour="TF tumor expression/TF normal expression"))+
  geom_point(aes(y=log2(expRatio),group="mark",colour="TF tumor expression/TF normal expression"),fill="transparent",shape=21,size=2)+
  labs(y="log2(tumor/normal)",x="Transcription factors (TFs)")+
  scale_fill_manual(name="",values=c("Tumor-specific enrichment/shared enrichment"="grey"))+
  scale_colour_manual(name="",values=c("TF tumor expression/TF normal expression"=COL[5]))+
  theme_classic(base_size=20)+
  theme(axis.text.x=element_text(angle=90, vjust=0.5,hjust=1,size = 10),
        legend.position=c(0.2,0.95), legend.title=element_blank(),legend.spacing = unit(-0.5,"cm"),legend.background = element_blank(),
        plot.title = element_text(hjust=0.5,face="bold"))
ggsave(gg_combined_tf_te, filename="functional_enrichment_combined_TE_eQTLs_TFexpression_k100_maf002_w2500.pdf",device="pdf",height=7,width=20)

df$log2expRatio <- log2(df$expRatio)




ggscat_te <- ggscatter(df[df$mark %in% toKeep_te,], x = "ratio", y = "log2expRatio", 
          add = "reg.line", conf.int = FALSE, 
          cor.coef = TRUE, cor.method = "pearson",add.params = list(color="red"),
          xlab = "log2(tumor-specific / shared TE-eQTL enrichment)", ylab = "log2(tumor / normal TE expression)")

ggscat_te <- ggplot(data=df[df$mark %in% toKeep_te,], aes(x= ratio, y = log2expRatio))+
  geom_point()+
  geom_smooth(se=FALSE,method="lm",color="red")+
  stat_cor(method="spearman",cor.coef.name = "rho",size=3)+
  labs(x="log2(tumor-specific / shared TE-eQTL enrichment)", y="log2(tumor / normal TE expression)")+
  ggthemes::theme_base(base_size=12)

ggsave(ggscat_te, filename="tf_exp_vs_enrichment_tspe_shared_TE_eQTL.pdf",height=5,width=5.5,device="pdf")





# CHECKING WHETHER STRONGER ENRICHMENT TUMOR SPECIFIC QTLs ARE INACTIVE / ACTIVE FOR OTHER FEATURES ----------------------------
vtf <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_LoVo/VARIANT_LoVo_ChIPseq_overlap.bed.gz",stringsAsFactors=FALSE,sep="\t",header=FALSE))
vtf$V9 <- gsub("LoVo_|_peak.*","", vtf$V9)
QTL <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/specific/tumor_specific/TUMOR.specific_TE_eQTLs.txt.gz",header=TRUE,stringsAsFactors=FALSE,sep=" "))
NQTL <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/mapping/normal/NORMAL.conditional_All.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep= " "))

comb_te$ratio <- log2(comb_te$OddsRatio.x / comb_te$OddsRatio.y)
tspe_strong_enrich <- gsub("LoVo ","",comb_te[which(comb_te$ratio >0 & comb_te$mark %in% toKeep_te),]$mark)

df <- vtf[vtf$V9 %in% tspe_strong_enrich,]
df <- merge(QTL,df,by.x="var_id",by.y="V4")
df$qtl <- paste(df$var_id, df$phe_id, sep=";")

nn <- data.frame()
for(i in unique(df$var_id))
{
  rrr <- NQTL[which(NQTL$var_id == i),]
  if(nrow(rrr) != 0)
  {
    if(length(unique(rrr$bwd_sig)) > 1)
    {
      nn <- rbind(nn, data.frame("var_id"=i, "bwd_sig"=1))
    }else
      {
        nn <- rbind(nn, data.frame("var_id"=i, "bwd_sig"=unique(rrr$bwd_sig)))
      }
  }else{
    nn <- rbind(nn, data.frame("var_id"=i, "bwd_sig"=0))
  }
}




pdf("functional_enrichment_tumor_spe_eQTL_stronger_eqtl_activity_in_normal.pdf",6,5.5)
bp<-barplot(table(nn$bwd_sig),main="tumor-specific TE-eQTLs with stronger enrichment\nfor TFs than shared TE eQTLs\nactivity in normal for other features",names.arg = c("Inactive\n (genome wide 5% FDR)","Active\n(genome wide 5% FDR)"))
balance <- table(nn$bwd_sig)
text(bp, balance/2, labels = round(balance, digits = 2))
dev.off()

nntf <- merge(nn, vtf[vtf$V9 %in% tspe_strong_enrich ,], by.x="var_id", by.y= "V4")
setdiff(unique(nntf[which(nntf$bwd_sig == 0),]$V9),unique(nntf[which(nntf$bwd_sig == 1),]$V9))
#[1] "FOXG1" "KLF3"  "IRX3"  "SMAD3" "SRF"  ==> PREVIOUS
#[1] "SOX4"  "TCF12" "HES1"

active_tspe_tes <- QTL[QTL$var_id %in% nn[which(nn$bwd_sig == 1),]$var_id ,]$phe_id
inactive_tspe_tes <- QTL[QTL$var_id %in% nn[which(nn$bwd_sig == 0),]$var_id ,]$phe_id











# CHECKING DIFFERENTIAL EXPRESSION FOR THESE TES ---------------------------------------------------------------------------------
dete <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/analysis_DGE/DGE_results_DESEq2_all_TE_genes.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
dete <- dete[which(dete$padj <= 0.05),]
dete$FC <- dete$t_median/ dete$n_median

detetf_act <- data.frame()
detetf_inact <- data.frame()
x <- dete[dete$id %in% active_tspe_tes ,]
y <- dete[dete$id %in% inactive_tspe_tes ,]

bed <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/analysis_DGE/SYSCOL_normal_275_tumor_276.norm.counts.bed",header=TRUE,stringsAsFactors = FALSE,sep="\t"))


for(i in 1:nrow(x))
{
  r <- x[i,]
  n <- unlist(strsplit(r$normal_exp, "\\,"))
  t <- unlist(strsplit(r$tumor_exp, "\\,"))
  ne <- data.frame("id"=replicate(length(n),r$id),"exp"=n,"tissue"=replicate(length(n),"normal"))
  te <- data.frame("id"=replicate(length(t),r$id), "exp"=t, "tissue"=replicate(length(t), "tumor"))
  detetf_act <- rbind(detetf_act, rbind(ne,te))
}
detetf_act$exp <- as.numeric(as.character(detetf_act$exp))

detetf_act <- merge(detetf_act, QTL[,c(1,8)],by.x="id",by.y="phe_id")
detetf_act <- merge(detetf_act, nn, by="var_id")

for(i in 1:nrow(y))
{
  r <-y[i,]
  n <- unlist(strsplit(r$normal_exp, "\\,"))
  t <- unlist(strsplit(r$tumor_exp, "\\,"))
  ne <- data.frame("id"=replicate(length(n),r$id),"exp"=n,"tissue"=replicate(length(n),"normal"))
  te <- data.frame("id"=replicate(length(t),r$id), "exp"=t, "tissue"=replicate(length(t), "tumor"))
  detetf_inact <- rbind(detetf_inact, rbind(ne,te))
}
detetf_inact$exp <- as.numeric(as.character(detetf_inact$exp))

detetf_inact <- merge(detetf_inact, QTL[,c(1,8)],by.x="id",by.y="phe_id")
detetf_inact <- merge(detetf_inact, nn, by="var_id")



gg_act <- ggplot(data=detetf_act[which(detetf_act$bwd_sig == 1),], aes(x=id, y=log2(exp), fill=tissue))+
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means(method="wilcox.test",label = "p.signif")+
  labs(x="TEs with tumor-specific eQTLs\noverlapping stronger enriched TFs compared to shared", y = "log2 normalized expression")+
  facet_wrap(~bwd_sig)+
  scale_fill_manual(values=c("normal"=COL[2], "tumor"=COL[1]))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        legend.position="top",panel.grid=element_line(linetype = "dashed"))
ggsave(gg_act, filename="tes_with_tspe_eQTLs_stronger_enrichment_exp_diff_active_eQTL_normal.pdf",height=8,width=15)


gg_inact <- ggplot(data=detetf_inact[which(detetf_inact$bwd_sig == 0),], aes(x=id, y=log2(exp), fill=tissue))+
  geom_boxplot(outlier.shape = NA)+
  stat_compare_means(method="wilcox.test",label = "p.signif")+
  labs(x="TEs with tumor-specific eQTLs\noverlapping stronger enriched TFs compared to shared", y = "log2 normalized expression")+
  facet_wrap(~bwd_sig)+
  scale_fill_manual(values=c("normal"=COL[2], "tumor"=COL[1]))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        legend.position="top",panel.grid=element_line(linetype = "dashed"))
ggsave(gg_inact, filename="tes_with_tspe_eQTLs_stronger_enrichment_exp_diff_inactive_eQTL_normal.pdf",height=8,width=20)










# CHECKING METHYLATION LEVELS FOR STRONGER ENRICHMENT FOR TUMOR-SPECIFIC EQTLs --------------------------------------------------

vmeth <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_methylation/data_intersection/syscol_variants_methylations_2.5kbwindow.txt.gz",header=FALSE,sep="\t",stringsAsFactors = FALSE,fill=TRUE))
meth <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_methylation/data_methylation/syscol_Norm_Betas.3-6-2015.bed.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
vmeth[210274,]$V5 <- 1
vmeth[210274,]$V6 <- 62117263
vmeth[210274,]$V7 <- 62122264
vmeth[210274,]$V8 <- "cg00786618"
vmeth[210274,]$V9 <- "." 
vmeth[210274,]$V10 <- "+" 

inact_tspe_var <- nn[which(nn$bwd_sig == 0),]$var_id
act_tspe_var <- nn[which(nn$bwd_sig == 1),]$var_id

inact <- vmeth[vmeth$V4 %in% inact_tspe_var ,]
act <- vmeth[vmeth$V4 %in% act_tspe_var ,]

n_samples <- as.character(read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_expression/stat/syscol_normal_samples.txt",header=F,stringsAsFactors = FALSE))
t_samples <- as.character(read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_expression/stat/syscol_tumor_samples.txt",header=F,stringsAsFactors = FALSE))

res_act <- foreach(s=unique(act_tspe_var), .combine = rbind, .multicombine = T) %do% 
{
  r <- act[which(act$V4 == s),]
  x <- meth[meth$id %in% r$V8,]
  n <- names(x)[grepl("N",names(x))]
  t <- names(x)[!grepl("N",names(x))]
  n <- intersect(n, n_samples)
  t <- intersect(t, t_samples)
  
  normal <- x %>% select(n)
  tumor <- x %>% select(t)
  
  n.med <- as.numeric(apply(normal, 2, median, na.rm=TRUE))
  t.med <- as.numeric(apply(tumor, 2, median, na.rm=TRUE))
  test <- try(wilcox.test(t.med, n.med))
  if(class(test) == "try-error"){
    pval <- NA
  }else{
    pval <- test$p.value
  }
  return(data.frame("snp"=s,
                    "wilcox_pval"=pval,
                    "abs_med_dif"=abs(median(n.med,na.rm=TRUE)-median(t.med,na.rm=TRUE)),
                    "median_normal"=median(n.med,na.rm=TRUE),
                    "median_tumor"=median(t.med,na.rm=TRUE),
                    "n_med"=paste(n.med,collapse=","),
                    "t_med"=paste(t.med,collapse=",")))
}

res_inact <- foreach(s=unique(inact_tspe_var), .combine = rbind, .multicombine = T) %do% 
  {
    r <- inact[which(inact$V4 == s),]
    x <- meth[meth$id %in% r$V8,]
    n <- names(x)[grepl("N",names(x))]
    t <- names(x)[!grepl("N",names(x))]
    n <- intersect(n, n_samples)
    t <- intersect(t, t_samples)
    
    normal <- x %>% select(n)
    tumor <- x %>% select(t)
    
    n.med <- as.numeric(apply(normal, 2, median, na.rm=TRUE))
    t.med <- as.numeric(apply(tumor, 2, median, na.rm=TRUE))
    test <- try(wilcox.test(t.med, n.med))
    if(class(test) == "try-error"){
      pval <- NA
    }else{
      pval <- test$p.value
    }
    return(data.frame("snp"=s,
                      "wilcox_pval"=pval,
                      "abs_med_dif"=abs(median(n.med,na.rm=TRUE)-median(t.med,na.rm=TRUE)),
                      "median_normal"=median(n.med,na.rm=TRUE),
                      "median_tumor"=median(t.med,na.rm=TRUE),
                      "n_med"=paste(n.med,collapse=","),
                      "t_med"=paste(t.med,collapse=",")))
  }

par(mfrow=c(1,2))
hist(res_act$wilcox_pval,xlim=c(0,1))
hist(res_inact$wilcox_pval,xlim=c(0,1))

res_act$activity <- "active"
res_inact$activity <- "inactive"
res_both <- rbind(res_act, res_inact)

gg <- ggplot(data=res_both, aes(x=activity, y=abs_med_dif))+
  geom_boxplot()+
  geom_jitter(width=0.3,shape=21)+
  stat_compare_means(method="wilcox.test")+
  labs(x="tumor specific TE-eQTL activity in normal for other features", y="Absolute median methylation difference\nbetween normal and tumor")+
  theme_bw()+
  theme(panel.grid=element_blank())
ggsave(gg, filename="tumor_specific_TE-eQTL_activity_abs_methylation_dif.pdf",height=6,width=6.5,device="pdf")


#************************************#
#   Functional enrichment for GENES  # 
#************************************#

# Individual plots
tspe_gene$sig <- tspe_gene$fdr <= 0.05
tspe_gene$sig10 <- tspe_gene$fdr <= 0.1
tspe_gene$nom_sig <- tspe_gene$pvalue < 0.05
tspe_gene_gg <- ggplot(data=tspe_gene[tspe_gene$nom_sig,], aes(x=reorder(mark,-OddsRatio),y=OddsRatio,label=round(-log10(pvalue),1)))+
  geom_bar(stat="identity")+
  geom_hline(yintercept=1, lty=2)+
  scale_fill_brewer(palette="Set2")+
  theme_classic(base_size=20)+
  labs(x="",y="Enrichment over the null",title="Functional enrichment tumor-specific gene-eQTLs")+
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1,vjust=0.5),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.85,0.9),legend.title=element_blank(),legend.background = element_rect(color="black",size=0.5),legend.margin=margin(t = 0,b=0.1,r=0.1,l=0.1, unit='cm'))
ggsave(tspe_gene_gg, filename="functional_enrichment_tspe_gene_LoVo_CS_CM_nominal_k100_maf002_w2500.pdf",device="pdf",height=7,width=15)

shared_gene$sig <- shared_gene$fdr <= 0.05
shared_gene$sig10 <- shared_gene$fdr <= 0.1
shared_gene$nom_sig <- shared_gene$pvalue < 0.05
shared_gene_gg <- ggplot(data=shared_gene[shared_gene$nom_sig,], aes(x=reorder(mark,-OddsRatio),y=OddsRatio,label=round(-log10(fdr),1)))+
  geom_bar(stat="identity")+
  geom_hline(yintercept=1, lty=2)+
  theme_classic(base_size=20)+
  scale_fill_brewer(palette="Set2")+
  labs(x="",y="Enrichment over the null",title="Functional enrichment shared gene-eQTLs")+
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1,vjust=0.5),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.85,0.9),legend.title=element_blank(),legend.background = element_rect(color="black",size=0.5),legend.margin=margin(t = 0,b=0.1,r=0.1,l=0.1, unit='cm'))
ggsave(shared_gene_gg, filename="functional_enrichment_shared_gene_LoVo_CS_CM_nominal_k100_maf002_w2500.pdf",device="pdf",height=7,width=15)



# combined tumor-spe and shared plot
comb_gene <- merge(tspe_gene, shared_gene, by="mark")

t <- comb_gene

toKeep2 <- filterEnrichment(t)



TMP_tspe_gene <- tspe_gene[tspe_gene$mark %in% toKeep2,c(1,7,8)]
TMP_tspe_gene$type <- "tumor-specific"

TMP_shared_gene <- shared_gene[shared_gene$mark %in% toKeep2, c(1,7,8)]
TMP_shared_gene$type <- "shared"

TMP_gene <- rbind(TMP_shared_gene,TMP_tspe_gene)
head(TMP_gene)
TMP_gene$OddsRatio[which(TMP_gene$OddsRatio >30)] <- 30
TMP_gene$ypos <- 0
for(i in unique(TMP_gene$mark)){
  x <- TMP_gene[which(TMP_gene$mark == i),]
  mm <- max(x$OddsRatio)
  TMP_gene[which(TMP_gene$mark == i & TMP_gene$type == "tumor-specific"),]$ypos <- mm + 1.5
  TMP_gene[which(TMP_gene$mark == i & TMP_gene$type == "shared"),]$ypos <- mm + 4
}


subset_to_order <- subset(TMP_gene,type == "tumor-specific")
subset_to_order$mark <- with(subset_to_order, reorder(mark,-OddsRatio))
TMP_gene$mark = factor(TMP_gene$mark, levels = levels(subset_to_order$mark))

head(TMP_gene)
gg_combined_gene <- ggplot(data=TMP_gene, aes(x=mark, y=OddsRatio, fill=type,label=round(-log10(pvalue),1)))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_hline(yintercept=1, lty=2, colour="grey")+
  geom_text(aes(x=mark, y=ypos),angle=90,size=3)+
  theme_classic(base_size=20)+
  labs(x="",y="Enrichment over the null",title="Functional enrichment for tumor-specific and shared gene-eQTLs")+
  scale_fill_manual(values=c("tumor-specific"="#e41a1c", "shared"="black"))+
  theme(axis.text.x = element_text(size=9,angle=90,hjust=1,vjust=0.5),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.9,1), legend.title=element_blank())
ggsave(gg_combined_gene, filename="functional_enrichment_combined_tspe_shared_gene_eQTLs_k100_maf002_w2500.pdf",device="pdf",height=7,width=20)


# comparison tumor-specific and shared enrichment

comb_gene$ratio <- log2(comb_gene$OddsRatio.x / comb_gene$OddsRatio.y)

gg_comb_gene1 <- ggplot(data=comb_gene[comb_gene$mark %in% toKeep2,], aes(x=reorder(mark,-ratio), y=ratio))+
  geom_bar(stat="identity",color="black",fill="grey")+
  theme_classic(base_size=20)+
  labs(x="",y="Log2(Tumor-specific enrichment/shared enrichment)")+
  theme(axis.text.x = element_text(size=9,angle=90,hjust=1,vjust=0.5),
        axis.title.y = element_text(size=15),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.9,0.9), legend.title=element_blank())
ggsave(gg_comb_gene1, filename="functional_enrichment_ratio_tspe_shared_gene_eQTLs_k100_maf002_w2500.pdf",device="pdf",height=7,width=20)  



comb_gene$geneName <- sapply(strsplit(comb_gene$mark, "\\ "),"[[",2)

df <- merge(comb_gene, tf, by="geneName")
COL=RColorBrewer::brewer.pal(9,"Set1")
gg_combined_tf_gene <- ggplot(data=df[df$mark %in% toKeep2,], aes(x=reorder(mark, -ratio), y=ratio, fill="Tumor-specific enrichment/shared enrichment"))+
  geom_bar(stat="identity",color="grey",width=0.8)+
  geom_line(aes(y=log2(expRatio), group="mark",colour="TF tumor expression/TF normal expression"))+
  geom_point(aes(y=log2(expRatio),group="mark",colour="TF tumor expression/TF normal expression"),fill="transparent",shape=21,size=2)+
  labs(y="log2(tumor/normal)",x="Transcription factors (TFs)")+
  scale_fill_manual(name="",values=c("Tumor-specific enrichment/shared enrichment"="grey"))+
  scale_colour_manual(name="",values=c("TF tumor expression/TF normal expression"=COL[5]))+
  theme_classic(base_size=20)+
  theme(axis.text.x=element_text(angle=90, vjust=0.5,hjust=1,size = 10),
        legend.position=c(0.2,0.95), legend.title=element_blank(),legend.spacing = unit(-0.5,"cm"),legend.background = element_blank(),
        plot.title = element_text(hjust=0.5,face="bold"))
ggsave(gg_combined_tf_gene, filename="functional_enrichment_combined_gene_eQTLs_TFexpression_k100_maf002_w2500.pdf",device="pdf",height=7,width=20)




#### TESTING null distribution --------------------------------------------------------
pdf("null_distribution_random_seeds_1000permutation.pdf",height=4.5, width=12)
par(mfrow=c(1,3))
all <- data.frame()
for(i in c(10,20,30,40,50,60,70,80,90,100,150,200,300,400,500)){
ran <- Sys.glob(paste0("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/enrichment/analysis_functionalEnrichment/V1/specific/tumor_specific/dif_seed/rand",i,"/*.txt"))

res <- data.frame()
for(rrr in 1:length(ran))
{
  cat(rrr,"\t")
  x <- ran[rrr]
  seed <- gsub(".txt","",gsub(".*_","",basename(x)))
  df <- read.table(x,header=TRUE,stringsAsFactors = FALSE,sep="\t")
  df <- performEnrichment(df)
  df$seed <- seed
  res <- rbind(res, df)
}
res$k <- i

hist(res$B,main=paste("Random hits",i,sep="\n"), xlab="Random hits")
hist(res$OddsRatio, main="OddsRatio", xlab="OddsRatio")
hist(res$pvalue, main="Pvalue", xlab="Pvalue")
all <- rbind(all, res)
}
dev.off()

df <- data.frame()
for(i in c(10,20,30,40,50,60,70,80,90,100,150,200,300,400,500))
{
  x <- all[which(all$k == i),]
  mean_rhits <- mean(x$B)
  mean_or <- mean(x$OddsRatio)
  mean_pv <- mean(x$pvalue)
  
  ymin_rhits <- min(x$B)
  ymax_rhits <- max(x$B)
  
  ymin_or <- min(x$OddsRatio)
  ymax_or <- max(x$OddsRatio)
  
  ymin_pv <- min(x$pvalue)
  ymax_pv <- max(x$pvalue)
  
  df <- rbind(df, data.frame("k"=i,"mean_rhits"=mean_rhits, "mean_or"=mean_or, "mean_pv"=mean_pv, "ymin_rhits"=ymin_rhits,"ymax_rhits"=ymax_rhits, "ymin_or"=ymin_or, "ymax_or"=ymax_or, "ymin_pv"=ymin_pv, "ymax_pv"=ymax_pv))  
}

gg1<-ggplot()+
  geom_pointrange(data=df, aes(x=k, y=mean_rhits,ymin=ymin_rhits, ymax=ymax_rhits), width=0.2, fill="white",shape=21)+
  labs(x="Random SNPS picked per eQTL", y="Random Hits", title="Random hits")+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5, face="bold"))

gg2<- ggplot()+
  geom_pointrange(data=df, aes(x=k, y=mean_or,ymin=ymin_or, ymax=ymax_or), width=0.2, fill="white",shape=21)+
  labs(x="Random SNPS picked per eQTL", y="Odds ratio", title="Odds ratio")+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5, face="bold"))

gg3<- ggplot()+
  geom_pointrange(data=df, aes(x=k, y=mean_pv,ymin=ymin_pv, ymax=ymax_pv), width=0.2, fill="white",shape=21)+
  labs(x="Random SNPS picked per eQTL", y="Pvalue", title="Pvalue")+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5, face="bold"))

gg<-ggarrange(gg1,gg2,gg3, ncol=3,nrow=1)
gg <- annotate_figure(gg, top=text_grob("Functional enrichment with different N SNPs picked (1000 times with different seed)",face="bold"))
gg
ggsave(gg, filename="null_distribution_accuracy_fenrich.pdf",height=5,width=15,device="pdf")




#### N random used 
all_n <- data.frame()
for(i in c(10,20,30,40,50,60,70,80,90,100,150,200,300,400,500))
{
  xxx <- Sys.glob(paste0("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/enrichment/analysis_functionalEnrichment/V1/specific/tumor_specific/dif_seed/rand",i,"/*used"))

  
  rrr <- 1
  x <- xxx[rrr]
  df <- read.table(x,header=TRUE,stringsAsFactors = FALSE,sep="\t")
    
  mean_n <- mean(df$n_random)
  min_n <- min(df$n_random)
  max_n <- max(df$n_random)
  res <- data.frame("k"=i, "mean_n"=mean_n, "min_n"=min_n, "max_n"=max_n)
  all_n <- rbind(all_n, res)
}  



k300 <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/enrichment/analysis_functionalEnrichment/V1/specific/tumor_specific/dif_seed/LoVo_ADNP_fenrich_tumor_specific_TE.results_seed_13976_k100_w5000.txt_null_distribution_used",header=TRUE,stringsAsFactors = FALSE,sep="\t"))

par(mfrow=c(2,2))
hist(k300$n_random,main="number of random SNPs per eQTL")
hist(k300$all_random, main="number of all matching SNPs per eQTL")
nrow(k300[which(k300$n_random < 100),])

k300 <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/enrichment/analysis_functionalEnrichment/V1/specific/tumor_specific/dif_seed/rand100/LoVo_ADNP_fenrich_tumor_specific_TE.results_seed_13976.txt_null_distribution_used",header=TRUE,stringsAsFactors = FALSE,sep="\t"))

hist(k300$n_random,main="number of random SNPs per eQTL")
hist(k300$all_random, main="number of all matching SNPs per eQTL")
nrow(k300[which(k300$n_random<100),])



dev.off()

nrow(k300)
81/344
