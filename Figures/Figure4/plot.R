#!/usr/bin/env Rscript 


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

n_gene <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/te_project/Figures/OUTLINE/B/n_gene_enrichments.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
t_gene <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/te_project/Figures/OUTLINE/B/t_gene_enrichments.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
n_te <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/te_project/Figures/OUTLINE/B/n_te_enrichments.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
t_te <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/te_project/Figures/OUTLINE/B/t_te_enrichments.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))



comb_n <- merge(n_gene, n_te,by="mark")
t <- comb_n
toKeep_n <- filterEnrichment(t)


comb_n$ratio <- log2(comb_n$OddsRatio.y/comb_n$OddsRatio.x)
toKeep_n2 <- toKeep_n[!toKeep_n %in% comb_n[which(comb_n$OddsRatio.x < 1 & comb_n$OddsRatio.y < 1),]$mark]

ggcomb_n <- ggplot(data=comb_n[comb_n$mark %in% toKeep_n2,], aes(x=reorder(mark, -ratio), y=ratio))+
  geom_bar(stat="identity", fill="grey", color="black")+
  labs(x="Transcription factors (TFs)", y = "log2(TE enrichment / gene enrichment)", title="Functional enrichment of TE and gene eQTLs in normal")+
  theme_classic(base_size=20)+
  theme(axis.text.x = element_text(size=9,angle=90,hjust=1,vjust=0.5),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.5,0.95), legend.title=element_blank(),
        panel.grid=element_blank())
ggsave(ggcomb_n, filename="te_gene_eqtl_functional_enrichment_in_normal.pdf",device="pdf", height=7, width=24)

## TUMOR ## 

comb_t <- merge(t_gene, t_te,by="mark")
t <- comb_t

toKeep_t <- filterEnrichment(t)

comb_t$ratio <- log2(comb_t$OddsRatio.y / comb_t$OddsRatio.x)
toKeep_t2 <- toKeep_t[!toKeep_t %in% comb_t[which(comb_t$OddsRatio.x < 1 & comb_t$OddsRatio.y < 1),]$mark]

ggcomb_t <- ggplot(data=comb_t[comb_t$mark %in% toKeep_t2,], aes(x=reorder(mark, -ratio), y=ratio))+
  geom_bar(stat="identity", fill="grey", color="black")+
  labs(x="Transcription factors (TFs)", y = "log2(TE enrichment / gene enrichment)", title="Functional enrichment of TE and gene eQTLs in tumor")+
  theme_classic(base_size=20)+
  theme(axis.text.x = element_text(size=9,angle=90,hjust=1,vjust=0.5),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.5,0.95), legend.title=element_blank(),
        panel.grid=element_blank())
ggsave(ggcomb_t, filename="te_gene_eqtl_functional_enrichment_in_tumor.pdf",device="pdf",height=7,width=24)

##### TISSUE SPECIFIC #####


tspe_te <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/te_project/Figures/OUTLINE/B/tumor_specific_TEeQTL_enrichment_results_maf002_w2500_k100.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
shared_te <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/te_project/Figures/OUTLINE/B/shared_TEeQTL_enrichment_results_maf002_w2500_k100.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
tspe_gene <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/te_project/Figures/OUTLINE/B/tumor_specific_GENEeQTL_enrichment_results_maf002_w2500_k100.txt",sep="\t",header=TRUE, stringsAsFactors = FALSE)
shared_gene <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/te_project/Figures/OUTLINE/B/shared_GENEeQTL_enrichment_results_maf002_w2500_k100.txt",sep="\t",header=TRUE, stringsAsFactors = FALSE)

#************************************#
#    Functional enrichment for TEs   # 
#************************************#

# Individual plots
tspe_te$sig <- tspe_te$fdr <= 0.05
tspe_te$sig10 <- tspe_te$fdr <= 0.1
shared_te$sig <- shared_te$fdr <= 0.05
shared_te$sig10 <- shared_te$fdr <= 0.1

# combined tumor-spe and shared plot
comb_te <- merge(tspe_te, shared_te, by="mark")
t <- comb_te 

toKeep_te <- filterEnrichment(t)
toKeep <- t[which(t$sig.x | t$sig.y),]$mark

comb_te$ratio <- log2(comb_te$OddsRatio.x / comb_te$OddsRatio.y)
gg_comb_te1 <- ggplot(data=comb_te[comb_te$mark %in% toKeep_te,], aes(x=reorder(mark,-ratio), y=ratio))+
  geom_bar(stat="identity",color="black",fill="grey")+
  geom_hline(yintercept=0,lty=2,size=0.2)+
  theme_classic(base_size=20)+
  labs(x="",y="Log2(tumor-specific enrichment/shared enrichment)", title="functional enrichment of tumor-specific and shared TE-eQTLs")+
  scale_fill_brewer(palette="Set2",direction = -1)+
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1,vjust=0.5),
        axis.title.y=element_text(size=15),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.85,0.9),legend.title=element_blank(),legend.background = element_rect(color="black",size=0.5),legend.margin=margin(t = 0,b=0.1,r=0.1,l=0.1, unit='cm'))

ggsave(gg_comb_te1, filename="functional_enrichment_ratio_tspe_shared_TE_eQTLs_tumor_spe_FDR5_only_k100_maf002_w2500.pdf",device="pdf",height=7.5,width=20)  
