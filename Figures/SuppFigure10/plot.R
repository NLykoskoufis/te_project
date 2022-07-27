#!/usr/bin/env Rscript 

library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
setwd("/Users/nikolaoslykoskoufis/Documents/PROJECTS/te_project/paper/rerun/Figures/SuppFigure10")
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

shared_gene <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/specific/shared/shared_GENEeQTL_enrichment_results_maf002_w2500_k100.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
tspe_gene <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/specific/tumor_specific/tumor_specific_GENEeQTL_enrichment_results_maf002_w2500_k100.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))

# Individual plots
tspe_gene$sig <- tspe_gene$fdr <= 0.05
tspe_gene$sig10 <- tspe_gene$fdr <= 0.1
tspe_gene$nom_sig <- tspe_gene$pvalue < 0.05
shared_gene$sig <- shared_gene$fdr <= 0.05
shared_gene$sig10 <- shared_gene$fdr <= 0.1
shared_gene$nom_sig <- shared_gene$pvalue < 0.05


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

