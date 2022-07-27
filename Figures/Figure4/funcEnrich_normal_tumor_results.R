#!/usr/bin/env Rscript 

setwd("~/Documents/PROJECTS/te_project/paper/rerun/Figures/Figure4/")
### USEFUL FUNCTIONS #### 
COL <- RColorBrewer::brewer.pal(8,"Dark2")
col = list("TE"=COL[3], "GENE"=COL[1])
combine_results <- function(x)
{
  #### x -> list of files obtained with fenrich to combine together
  
  df <- data.frame()
  for(file in x)
  {
    r <- as.data.frame(data.table::fread(paste0(file), header=TRUE,stringsAsFactors = FALSE, sep="\t"))
    mark <- sapply(strsplit(basename(file),"\\_"), "[[",1)
    df <- rbind(df, cbind("mark"=mark, r))
  }
  return(df)
}

### ENSEMBL REGULATORY BUILD ### 

combine_results <- function(x)
{
  #### x -> list of files obtained with fenrich to combine together
  fileCount <- 0
  df <- data.frame()
  for(file in x)
  {
    fileCount <- fileCount + 1
    cat(fileCount, "-",length(x),"\n")
    
    r <- as.data.frame(data.table::fread(paste0(file), header=TRUE,stringsAsFactors = FALSE, sep="\t"))
    mark <- sapply(strsplit(basename(file),"\\_"), "[[",1)
    tmp <- as.data.frame(data.table::fread(Sys.glob(paste0("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_Ensembl_Regulatory_Build/combined_peaks/",mark,".*.filtered.bed.gz")),header=FALSE,sep="\t",stringsAsFactors = FALSE))
    celltypes <- paste(unique(unlist(strsplit(unique(tmp$V5),"\\,"))),collapse=",")
    df <- rbind(df, cbind("mark"=mark,"celltype"=celltypes, r))
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
  
  t.sig.y.not.x <- t.sig.y.not.x[!t.sig.y.not.x$OddsRatio.x > t.sig.y.not.x$OddsRatio.y,]
  t.sig.x.not.y <- t.sig.x.not.y[!t.sig.x.not.y$OddsRatio.y > t.sig.x.not.y$OddsRatio.x,]
  toKeep <- unique(c(t.sig.both$mark, t.sig.x.not.y$mark, t.sig.y.not.x$mark))
  return(toKeep)
}


#### LIST OF FILES #### 
normal_gene <- Sys.glob("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/genes/normal/*enrichment.txt")
normal_te <- Sys.glob("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/tes/normal/*enrichment.txt")

tumor_gene <- Sys.glob("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/genes/tumor/*enrichment.txt")
tumor_te <- Sys.glob("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/tes/tumor/*enrichment.txt")

# COMBINE RESULTS #
n_gene <- combine_results(normal_gene)
n_te <- combine_results(normal_te)
t_gene <- combine_results(tumor_gene)
t_te <- combine_results(tumor_te)



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

n_gene <- performEnrichment(n_gene)
n_te <- performEnrichment(n_te)
t_gene <- performEnrichment(t_gene)
t_te <- performEnrichment(t_te)



 # MULTIPLE TEST CORRECTION # 
n_gene$fdr<- p.adjust(n_gene$pvalue, method="fdr")
n_te$fdr <- p.adjust(n_te$pvalue, method="fdr")
t_gene$fdr <- p.adjust(t_gene$pvalue, method="fdr")
t_te$fdr <- p.adjust(t_te$pvalue, method="fdr")

### WRITE RESULTS TO FILES ### 

write.table(n_gene,"/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/genes/normal/NORMAL.All.gene.functionalEnrichment.txt", sep="\t", col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(n_te,"/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/genes/normal/NORMAL.All.TE.functionalEnrichment.txt", sep="\t", col.names=TRUE,row.names=FALSE,quote=FALSE)

write.table(t_gene,"/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/genes/tumor/TUMOR.All.gene.functionalEnrichment.txt", sep="\t", col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(t_te,"/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/enrichment/genes/tumor/TUMOR.All.TE.functionalEnrichment.txt", sep="\t", col.names=TRUE,row.names=FALSE,quote=FALSE)

save(n_gene, n_te, t_gene, t_te, file="combined_results_normal_tumor_fenrich.RData")
load("combined_results_normal_tumor_fenrich.RData")



pdf("funcEnrich_normal_tumor_te_gene.pdf",height = 8.5,width=8)
par(mfrow=c(2,2))
hist(n_gene$pvalue)
hist(n_te$pvalue)
hist(t_gene$pvalue)
hist(t_te$pvalue)
dev.off()

n_gene$sig <- n_gene$fdr <= 0.05
n_te$sig <- n_te$fdr <= 0.05
t_gene$sig <- t_gene$fdr <= 0.05
t_te$sig <- t_te$fdr <= 0.05


comb_n <- merge(n_gene, n_te,by="mark")
t <- comb_n
toKeep_n <- filterEnrichment(t)


TMP_n_gene <- n_gene[n_gene$mark %in% toKeep_n,c(1,7,8)]
TMP_n_gene$type <- "gene"

TMP_n_te <- n_te[n_te$mark %in% toKeep_n, c(1,7,8)]
TMP_n_te$type <- "te"

TMP_n <- rbind(TMP_n_gene,TMP_n_te)

TMP_n$ypos <- 0
for(i in unique(TMP_n$mark)){
  x <- TMP_n[which(TMP_n$mark == i),]
  mm <- max(x$OddsRatio)
  TMP_n[which(TMP_n$mark == i & TMP_n$type == "gene"),]$ypos <- mm + 1.5
  TMP_n[which(TMP_n$mark == i & TMP_n$type == "te"),]$ypos <- mm + 4
}

subset_to_order <- subset(TMP_n,type == "te")
subset_to_order$mark <- with(subset_to_order, reorder(mark,-OddsRatio))
TMP_n$mark = factor(TMP_n$mark, levels = levels(subset_to_order$mark))


gg_combined_n <- ggplot(data=TMP_n, aes(x=mark, y=OddsRatio, fill=type,label=round(-log10(pvalue),1)))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_hline(yintercept=1, lty=2, colour="grey")+
  geom_text(aes(x=mark, y=ypos),angle=90,size=3)+
  theme_classic(base_size=20)+
  labs(x="",y="Enrichment over the null",title="Functional enrichment for normal gene and TE eQTLs")+
  scale_fill_manual(values=c("gene"=col$GENE, "te"=col$TE))+
  theme(axis.text.x = element_text(size=9,angle=90,hjust=1,vjust=0.5),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.5,0.95), legend.title=element_blank())
ggsave(gg_combined_n, filename="functional_enrichment_combined_normal_gene_TE_k100_maf002_w2500.pdf",device="pdf",height=7,width=24)

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

comb_nf <- comb_n[comb_n$mark %in% toKeep_n,]
comb_nf[which(comb_nf$ratio > 0),]$mark


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

ggarrange(gg_combined_t, ggcomb_t, ncol=1, nrow=2)

##### strongher enrichment TFs for TEs.
union(comb_n[which(comb_n$ratio > 0 & comb_n$mark %in% toKeep_n),]$mark,comb_t[which(comb_t$ratio > 0 & comb_t$mark %in% toKeep_t),]$mark)
intersect(comb_n[which(comb_n$ratio > 0 & comb_n$mark %in% toKeep_n),]$mark,comb_t[which(comb_t$ratio > 0 & comb_t$mark %in% toKeep_t),]$mark)

#18 TFs in total 
# 5 in normal #Gata2 NFATC1 PHB2 RXRA ZNF274
# 15 in tumor #Bdp1 Brf1 FOSL2 Gata2 HES1 NFE2L2 PAX8 RUNX1 SOX13 STAT5A TBX3 TCF7 TRIM28 XRCC4 ZBTB33
# 1 in common Gata2

########## COMBINED TUMOR AND NORMAL PLOTS TOGETHER ###########
# Allows to better compare results # 


comb_n <- merge(n_gene, n_te,by="mark")
comb_n$ratio <- log2(comb_n$OddsRatio.y/comb_n$OddsRatio.x)
comb_n$tissue <- "normal"
comb_t <- merge(t_gene, t_te, by="mark")
comb_t$ratio <- log2(comb_t$OddsRatio.y / comb_t$OddsRatio.x)
comb_t$tissue <- "tumor"

combToKeep <- union(toKeep_n2,toKeep_t2)

combinedNT <- rbind(comb_n, comb_t)

ggcombinedNT <- ggplot(data=combinedNT[combinedNT$mark %in% combToKeep,], aes(x=reorder(mark, -ratio), y = ratio, fill=tissue))+
  geom_bar(stat="identity", color="white",position="dodge")+
  labs(x="Transcription factors (TFs)", y = "log2(TE enrichment / gene enrichment)", title="Functional enrichment of TE and gene eQTLs in normal and tumor")+
  scale_fill_manual(values=c("normal"="#377eb8", "tumor"="#e41a1c"))+
  theme_base(base_size=14)+
  theme(axis.text.x = element_text(size=9,angle=90,hjust=1,vjust=0.5),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.5,0.95), legend.title=element_blank(),legend.background = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor= element_blank(), panel.grid.major.x= element_line(linetype="dotted",color="lightgrey",size=0.3))

ggsave(ggcombinedNT, filename="NORMAL.TUMOR_TE_genes_log2enrichments.pdf",device="pdf", height=7, width=24)




### NORMAL VS TUMOR ### 

comb_nt_gene <- merge(n_gene, t_gene,by="mark")
t <- comb_nt_gene
toKeep <- filterEnrichment(t)
#toKeep <- comb_nt_gene[which(comb_nt_gene$sig.y | comb_nt_gene$sig.x),]$mark


TMP_n_gene <- n_gene[n_gene$mark %in% toKeep,c(1,7,8)]
TMP_n_gene$type <- "normal"

TMP_t_gene <- t_gene[t_gene$mark %in% toKeep, c(1,7,8)]
TMP_t_gene$type <- "tumor"

TMP_nt_gene <- rbind(TMP_n_gene,TMP_t_gene)

TMP_nt_gene$ypos <- 0
for(i in unique(TMP_nt_gene$mark)){
  x <- TMP_nt_gene[which(TMP_nt_gene$mark == i),]
  mm <- max(x$OddsRatio)
  TMP_nt_gene[which(TMP_nt_gene$mark == i & TMP_nt_gene$type == "normal"),]$ypos <- mm + 1.5
  TMP_nt_gene[which(TMP_nt_gene$mark == i & TMP_nt_gene$type == "tumor"),]$ypos <- mm + 4
}

subset_to_order <- subset(TMP_nt_gene,type == "tumor")
subset_to_order$mark <- with(subset_to_order, reorder(mark,-OddsRatio))
TMP_nt_gene$mark = factor(TMP_nt_gene$mark, levels = levels(subset_to_order$mark))

COL2 <- RColorBrewer::brewer.pal(9, "Set1")
gg_combined_nt_gene <- ggplot(data=TMP_nt_gene, aes(x=mark, y=OddsRatio, fill=type,label=round(-log10(pvalue),1)))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_hline(yintercept=1, lty=2, colour="grey")+
  geom_text(aes(x=mark, y=ypos),angle=90,size=3)+
  theme_classic(base_size=20)+
  labs(x="",y="Enrichment over the null",title="Functional enrichment for normal/tumor gene eQTLs")+
  scale_fill_manual(values=c("tumor"=COL2[1], "normal"=COL2[2]))+
  theme(axis.text.x = element_text(size=9,angle=90,hjust=1,vjust=0.5),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.5,0.95), legend.title=element_blank(),legend.background = element_rect(color="black"),legend.margin = margin(t=-0.2,b=0.1,r=0.1,l=0.1, "cm"))
ggsave(gg_combined_nt_gene, filename="functional_enrichment_combined_normal_tumor_gene.pdf",device="pdf",height=7,width=20)


comb_nt_te <- merge(n_te, t_te,by="mark")
t <- comb_nt_te
toKeep <- filterEnrichment(t)
#toKeep <- comb_nt_te[which(comb_nt_te$sig.y | comb_nt_te$sig.x),]$mark



TMP_n_te <- n_te[n_te$mark %in% toKeep,c(1,7,8)]
TMP_n_te$type <- "normal"

TMP_t_te <- t_te[t_te$mark %in% toKeep, c(1,7,8)]
TMP_t_te$type <- "tumor"

TMP_nt_te <- rbind(TMP_n_te,TMP_t_te)

TMP_nt_te$ypos <- 0
for(i in unique(TMP_nt_te$mark)){
  x <- TMP_nt_te[which(TMP_nt_te$mark == i),]
  mm <- max(x$OddsRatio)
  TMP_nt_te[which(TMP_nt_te$mark == i & TMP_nt_te$type == "normal"),]$ypos <- mm + 1.5
  TMP_nt_te[which(TMP_nt_te$mark == i & TMP_nt_te$type == "tumor"),]$ypos <- mm + 7
}

subset_to_order <- subset(TMP_nt_te,type == "tumor")
subset_to_order$mark <- with(subset_to_order, reorder(mark,-OddsRatio))
TMP_nt_te$mark = factor(TMP_nt_te$mark, levels = levels(subset_to_order$mark))


gg_combined_nt_te <- ggplot(data=TMP_nt_te, aes(x=mark, y=OddsRatio, fill=type,label=round(-log10(pvalue),1)))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_hline(yintercept=1, lty=2, colour="grey")+
  geom_text(aes(x=mark, y=ypos),angle=90,size=3)+
  theme_classic(base_size=20)+
  labs(x="",y="Enrichment over the null",title="Functional enrichment for normal/tumor TE eQTLs")+
  scale_fill_manual(values=c("tumor"=COL2[1], "normal"=COL2[2]))+
  theme(axis.text.x = element_text(size=9,angle=90,hjust=1,vjust=0.5),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.5,0.95), legend.title=element_blank(),legend.background = element_rect(color="black"),legend.margin = margin(t=-0.2,b=0.1,r=0.1,l=0.1, "cm"))
ggsave(gg_combined_nt_te, filename="functional_enrichment_combined_normal_tumor_TE.pdf",device="pdf",height=7,width=20)


TMP_nt_gene$feature <- "gene"
TMP_nt_te$feature <- "te"


TMP_nt_all <- rbind(TMP_nt_te,TMP_nt_gene)
TMP_nt_all$feature <- factor(TMP_nt_all$feature,levels=c("te","gene"))
TMP_nt_all$OddsRatio[which(TMP_nt_all$OddsRatio >30)] <- 30 
head(TMP_nt_all)
gg_nt_all <- ggplot(data=TMP_nt_all,aes(x=mark, y=OddsRatio, fill=type))+
  geom_bar(stat="identity",position=position_dodge(),width=0.7)+
  geom_hline(yintercept=1, lty=2, colour="grey")+
  facet_grid(feature ~ .)+
  theme_bw(base_size=20)+
  labs(x="",y="Enrichment over the null",title="Functional enrichment for normal/tumor eQTLs")+
  scale_fill_manual(values=c("tumor"=COL2[1], "normal"=COL2[2]))+
  theme(axis.text.x = element_text(size=9,angle=90,hjust=1,vjust=0.5),
        panel.grid.minor=element_blank(),panel.grid.major.y =element_blank(),panel.grid.major.x= element_line(linetype = "dashed",size = 0.5),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.5,0.95), legend.title=element_blank(),legend.background = element_rect(color="black"),legend.margin = margin(t=-0.2,b=0.1,r=0.1,l=0.1, "cm"))
ggsave(gg_nt_all, filename="functional_enrichment_combined_normal_tumor_TE_GENE.pdf",device="pdf",height=10,width=24)
   

#### TF and TF expression --------------------------------------------------------------------------

tf <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/analysis_DEG/GENE.DEG.txt.gz",header=T,stringsAsFactors = F,sep="\t"))
tf$expRatio <- tf$tumor_median / tf$normal_median
convert <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_expression/stat/ensemblToGeneName.txt",header=F,sep=" ",stringsAsFactors = F)
colnames(convert) <- c("ensemblID","geneName")
tf <- merge(convert,tf,by.x="ensemblID",by.y="id")

#### changing names of tf genes to match ddte. 
tf[which(tf$ensemblID == "ENSG00000085276"),]$geneName <- "EVI1"
tf[which(tf$ensemblID == "ENSG00000125820"),]$geneName <- "NKX2.2"
tf <- tf[which(tf$qvalue <= 0.05),]
tf$geneName <- toupper(tf$geneName)


comb_nt_te$mark <- toupper(comb_nt_te$mark)
comb_nt_te <- merge(comb_nt_te, tf, by.x="mark", by.y="geneName")

comb_nt_te$orRatio <- comb_nt_te$OddsRatio.y / comb_nt_te$OddsRatio.x

gg_comb_nt_te <- ggplot(data=comb_nt_te, aes(x=reorder(mark, -log2(orRatio)), y=log2(orRatio),fill="TF tumor enrichment/TF normal enrichment"))+
  geom_bar(stat="identity",colour="black")+
  geom_line(aes(x=reorder(mark, -log2(orRatio)), y=log2(expRatio),colour="TF tumor expression/TF normal expression"),group="mark")+
  geom_point(aes(x=reorder(mark, -log2(orRatio)), y=log2(expRatio),colour="TF tumor expression/TF normal expression"),group="mark",shape=21,fill="white",size=1.5)+
  labs(y="log2(tumor/normal)",x="Transcription factors (TFs)", title="TF enrichment for TE-eQTLs in normal and tumor")+
  scale_fill_manual(name="",values=c("TF tumor enrichment/TF normal enrichment"="grey"))+
  scale_colour_manual(name="",values=c("TF tumor expression/TF normal expression"=COL[5]))+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90, vjust=0.5,hjust=1),
        legend.position=c(0.2,0.95), legend.title=element_blank(),legend.spacing = unit(-0.5,"cm"),legend.background = element_blank(),
        plot.title = element_text(hjust=0.5,face="bold"))
ggsave(gg_comb_nt_te, filename="functional_enrichment_tf_normal_tumor_te_tf_expression.pdf",device="pdf", height=8,width=20)


cor.test(log2(comb_nt_te$orRatio), log2(comb_nt_te$expRatio), method="pearson")
plot(log2(comb_nt_te$orRatio), log2(comb_nt_te$expRatio))



comb_nt_gene$mark <- toupper(comb_nt_gene$mark)
comb_nt_gene <- merge(comb_nt_gene, tf, by.x="mark", by.y="geneName")

comb_nt_gene$orRatio <- comb_nt_gene$OddsRatio.y / comb_nt_gene$OddsRatio.x

gg_comb_nt_gene <- ggplot(data=comb_nt_gene, aes(x=reorder(mark, -log2(orRatio)), y=log2(orRatio),fill="TF tumor enrichment/TF normal enrichment"))+
  geom_bar(stat="identity",colour="black")+
  geom_line(aes(x=reorder(mark, -log2(orRatio)), y=log2(expRatio),colour="TF tumor expression/TF normal expression"),group="mark")+
  geom_point(aes(x=reorder(mark, -log2(orRatio)), y=log2(expRatio),colour="TF tumor expression/TF normal expression"),group="mark",shape=21,fill="white",size=1.5)+
  labs(y="log2(tumor/normal)",x="Transcription factors (TFs)", title="TF enrichment for gene-eQTLs in normal and tumor")+
  scale_fill_manual(name="",values=c("TF tumor enrichment/TF normal enrichment"="grey"))+
  scale_colour_manual(name="",values=c("TF tumor expression/TF normal expression"=COL[5]))+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90, vjust=0.5,hjust=1),
        legend.position=c(0.2,0.95), legend.title=element_blank(),legend.spacing = unit(-0.5,"cm"),legend.background = element_blank(),
        plot.title = element_text(hjust=0.5,face="bold"))
ggsave(gg_comb_nt_gene, filename="functional_enrichment_tf_normal_tumor_gene_tf_expression.pdf",device="pdf", height=8,width=20)








#### Common peaks between LoVo and Ensembl Reg ####

normal_gene <- Sys.glob("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/enrichment/analysis_functionalEnrichment/fenrich_cpp/gene_fenrich/normal/ensembl_regulatory/*enrichment.txt")
normal_te <- Sys.glob("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/enrichment/analysis_functionalEnrichment/fenrich_cpp/te_fenrich/normal/ensembl_regulatory/*enrichment.txt")

tumor_gene <- Sys.glob("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/enrichment/analysis_functionalEnrichment/fenrich_cpp/gene_fenrich/tumor/LoVo/*_fenrich.results")
tumor_te <- Sys.glob("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/enrichment/analysis_functionalEnrichment/fenrich_cpp/te_fenrich/tumor/LoVo/*_fenrich.results")

combine_LoVo <- function(x)
{
  #### x -> list of files obtained with fenrich to combine together
  
  df <- data.frame()
  for(file in x)
  {
    r <- as.data.frame(data.table::fread(paste0(file), header=TRUE,stringsAsFactors = FALSE, sep="\t"))
    ff <- gsub("_fenrich.results","",basename(file))
    mark <- gsub("LoVo_","",ff)
    df <- rbind(df, cbind("mark"=mark, r))
  }
  return(df)
}


t_gene <- combine_LoVo(tumor_gene)
t_te <- combine_LoVo(tumor_te)
t_gene$fdr <- p.adjust(t_gene$pvalue, method="fdr")
t_te$fdr <- p.adjust(t_te$pvalue, method="fdr")

n_gene$mark <- toupper(n_gene$mark)
t_gene$mark <- toupper(t_gene$mark)
n_te$mark <- toupper(n_te$mark)
t_te$mark <- toupper(t_te$mark)

n_gene$sig <- n_gene$fdr <= 0.05
t_gene$sig <- t_gene$fdr <= 0.05
n_te$sig <- n_te$fdr <= 0.05
t_te$sig <- t_te$fdr <= 0.05

comb_nt_gene <- merge(n_gene, t_gene,by="mark")
t <- comb_nt_gene

toKeep <- filterEnrichment(t)


common_peaks <- intersect(intersect(n_gene$mark, t_gene$mark),toKeep)

ngene <- n_gene[n_gene$mark %in% common_peaks,]
tgene <- t_gene[t_gene$mark %in% common_peaks,]
ngene <- ngene[,-2]
ngene$tissue <- "normal"
tgene$tissue <- "tumor"
nt_gene <- rbind(ngene,tgene)

nt_gene$ypos <- 0
for(i in unique(nt_gene$mark)){
  x <- nt_gene[which(nt_gene$mark == i),]
  mm <- max(x$OddsRatio)
  nt_gene[which(nt_gene$mark == i & nt_gene$tissue == "normal"),]$ypos <- mm + 1
  nt_gene[which(nt_gene$mark == i & nt_gene$tissue == "tumor"),]$ypos <- mm + 3
}

subset_to_order <- subset(nt_gene,tissue == "tumor")
subset_to_order$mark <- with(subset_to_order, reorder(mark,-OddsRatio))
nt_gene$mark = factor(nt_gene$mark, levels = levels(subset_to_order$mark))

nt_gene$mark


gg_nt_gene<-ggplot(data=nt_gene, aes(x=mark,y=OddsRatio,fill=tissue,label=round(-log10(fdr),1)))+
  geom_bar(stat="identity",position=position_dodge(),width=0.7)+
  geom_hline(yintercept = 1,lty=2)+
  geom_text(aes(x=mark, y=ypos),angle=90,size=3)+
  scale_fill_manual(values=c("normal"="blue2","tumor"="red2"))+
  theme_bw(base_size=20)+
  labs(x="",y="Enrichment over the null",title="Functional enrichment for normal/tumor eQTLs")+
  theme(axis.text.x = element_text(size=9,angle=90,hjust=1,vjust=0.5),
        panel.grid.minor=element_blank(),panel.grid.major =element_blank(),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.9,0.9), legend.title=element_blank(),legend.background = element_rect(color="black",size=0.5),legend.margin = margin(t=-0.2,b=0.1,r=0.1,l=0.1, "cm"))
ggsave(gg_nt_gene, filename="functional_enrichment_LoVo_common_peaks_gene_eQTLs.pdf",device="pdf",height=6,width=15)



comb_nt_te <- merge(n_te, t_te,by="mark")
t <- comb_nt_te

toKeep <- filterEnrichment(t)

common_peaks <- intersect(intersect(n_te$mark, t_te$mark),toKeep)

nte <- n_te[n_te$mark %in% common_peaks,]
tte <- t_te[t_te$mark %in% common_peaks,]
nte <- nte[,-2]
nte$tissue <- "normal"
tte$tissue <- "tumor"
nt_te <- rbind(nte,tte)

nt_te$ypos <- 0
for(i in unique(nt_te$mark)){
  x <- nt_te[which(nt_te$mark == i),]
  mm <- max(x$OddsRatio)
  nt_te[which(nt_te$mark == i & nt_te$tissue == "normal"),]$ypos <- mm + 1
  nt_te[which(nt_te$mark == i & nt_te$tissue == "tumor"),]$ypos <- mm + 3
}

subset_to_order <- subset(nt_te,tissue == "tumor")
subset_to_order$mark <- with(subset_to_order, reorder(mark,-OddsRatio))
nt_te$mark = factor(nt_te$mark, levels = levels(subset_to_order$mark))




gg_nt_te<-ggplot(data=nt_te, aes(x=mark,y=OddsRatio,fill=tissue,label=round(-log10(fdr),1)))+
  geom_bar(stat="identity",position=position_dodge(),width=0.7)+
  geom_hline(yintercept = 1,lty=2)+
  geom_text(aes(x=mark, y=ypos),angle=90,size=3)+
  scale_fill_manual(values=c("normal"="blue2","tumor"="red2"))+
  theme_bw(base_size=20)+
  labs(x="",y="Enrichment over the null",title="Functional enrichment for normal/tumor TE-eQTLs")+
  theme(axis.text.x = element_text(size=9,angle=90,hjust=1,vjust=0.5),
        panel.grid.minor=element_blank(),panel.grid.major =element_blank(),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.9,0.9), legend.title=element_blank(),legend.background = element_rect(color="black",size=0.5),legend.margin = margin(t=-0.2,b=0.1,r=0.1,l=0.1, "cm"))
gg_nt_te

ggsave(gg_nt_gene, filename="functional_enrichment_LoVo_common_peaks_TE_eQTLs.pdf",device="pdf",height=6,width=15)



nt_gene$feature <- "gene"
nt_te$feature <- "te"
nt_all <- rbind(nt_te, nt_gene)

nt_all$feature <- factor(nt_all$feature,levels = c("te","gene"))
str(nt_all)
gg_nt_all <- ggplot(data=nt_all,aes(x=mark, y=OddsRatio, fill=tissue,label=round(-log10(fdr),1)))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_hline(yintercept=1, lty=2, colour="grey")+
  geom_text(aes(x=mark, y=ypos),angle=90,size=3)+
  facet_grid(feature ~ .)+
  theme_bw(base_size=20)+
  labs(x="",y="Enrichment over the null",title="Functional enrichment for normal/tumor eQTLs")+
  scale_fill_manual(values=c("tumor"="red2", "normal"="blue2"))+
  theme(axis.text.x = element_text(size=9,angle=90,hjust=1,vjust=0.5),
        panel.grid.minor=element_blank(),panel.grid.major.y =element_blank(),
        plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.9,0.95), legend.title=element_blank(),legend.background = element_rect(color="black"),legend.margin = margin(t=-0.2,b=0.1,r=0.1,l=0.1, "cm"))
ggsave(gg_nt_all, filename="functional_enrichment_combined_normal_tumor_TE_GENE_common_LoVo_Ensembl_peaks.pdf",device="pdf",height=10,width=15)
