#!/usr/bin/env Rscript 

setwd("/Users/nikolaoslykoskoufis/Documents/PROJECTS/te_project/paper/rerun/Figures/SuppFigure9")
### USEFUL FUNCTIONS #### 
COL <- RColorBrewer::brewer.pal(8,"Dark2")
col = list("TE"=COL[3], "GENE"=COL[1])

filterEnrichment <- function(t)
{
  
  t.sig.both <- t[which(t$sig.x == TRUE & t$sig.y == TRUE),]
  t.sig.x.not.y <- t[which(t$sig.x == TRUE & t$sig.y == FALSE),]
  t.sig.y.not.x <- t[which(t$sig.y == TRUE & t$sig.x == FALSE),]
  
  pb1 <- t.sig.y.not.x[t.sig.y.not.x$OddsRatio.x > t.sig.y.not.x$OddsRatio.y,]
  pb2 <- t.sig.x.not.y[t.sig.x.not.y$OddsRatio.y > t.sig.x.not.y$OddsRatio.x,]
  
  t.sig.y.not.x <- t.sig.y.not.x[!t.sig.y.not.x$OddsRatio.x > t.sig.y.not.x$OddsRatio.y,]
  t.sig.x.not.y <- t.sig.x.not.y[!t.sig.x.not.y$OddsRatio.y > t.sig.x.not.y$OddsRatio.x,]
  toKeep <- unique(c(t.sig.both$mark, t.sig.x.not.y$mark, t.sig.y.not.x$mark))
  return(toKeep)
}

n_gene <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/genes/normal/NORMAL.All.gene.functionalEnrichment.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
t_gene <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/genes/tumor/TUMOR.All.gene.functionalEnrichment.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
n_te <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/tes/normal/NORMAL.All.TE.functionalEnrichment.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
t_te <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/tes/tumor/TUMOR.All.TE.functionalEnrichment.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))

n_gene$sig <- n_gene$fdr <= 0.05
t_gene$sig <- t_gene$fdr <= 0.05
n_te$sig <- n_te$fdr <= 0.05
t_te$sig <- t_te$fdr <= 0.05

comb_t <- merge(t_gene, t_te,by="mark")
t <- comb_t

toKeep_t <- filterEnrichment(t)
#toKeep <- comb_t[which(comb_t$sig.y | comb_t$sig.x),]$mark

TMP_t_gene <- t_gene[t_gene$mark %in% toKeep_t,c(1,7,8)]
TMP_t_gene$type <- "gene"

TMP_t_te <- t_te[t_te$mark %in% toKeep_t, c(1,7,8)]
TMP_t_te$type <- "te"

TMP_t <- rbind(TMP_t_gene,TMP_t_te)

TMP_t$ypos <- 0
for(i in unique(TMP_t$mark)){
  x <- TMP_t[which(TMP_t$mark == i),]
  mm <- max(x$OddsRatio)
  TMP_t[which(TMP_t$mark == i & TMP_t$type == "gene"),]$ypos <- mm + 1.5
  TMP_t[which(TMP_t$mark == i & TMP_t$type == "te"),]$ypos <- mm + 6
}

subset_to_order <- subset(TMP_t,type == "te")
subset_to_order$mark <- with(subset_to_order, reorder(mark,-OddsRatio))
TMP_t$mark = factor(TMP_t$mark, levels = levels(subset_to_order$mark))


gg_combined_t <- ggplot(data=TMP_t, aes(x=mark, y=OddsRatio, fill=type,label=round(-log10(pvalue),1)))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_hline(yintercept=1, lty=2, colour="grey")+
  geom_text(aes(x=mark, y=ypos),angle=90,size=3)+
  theme_classic(base_size=20)+
  labs(x="",y="Enrichment over the null",title="Functional enrichment for tumor gene and TE eQTLs")+
  scale_fill_manual(values=c("gene"=col$GENE, "te"=col$TE))+
  theme(axis.text.x = element_text(size=9,angle=90,hjust=1,vjust=0.5),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.5,0.95), legend.title=element_blank())
ggsave(gg_combined_t, filename="functional_enrichment_combined_tumor_gene_TE.pdf",device="pdf",height=7,width=24)
