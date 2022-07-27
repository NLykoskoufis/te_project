#!/usr/bin/env Rscript 


##########################################################
ANALYSIS_ = "Get all correlations between eQTL-TE-Gene"  #
##########################################################
library(foreach)
library(doMC)
cpus <- 4
registerDoMC(cpus)


if (! file.exists("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/union/UNION_BN_ALL_correlations.txt"))
{
  BED.n <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/normal/SYSCOL_normal_275_gene_TE.FM05.resid.cpm.bed.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
  BED.t <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/tumor/SYSCOL_tumor_276_gene_TE.FM05.resid.cpm.bed.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
  
  UNION <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/union/UNION_normal_tumor_with_posteriors.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
  UNION$triplet <- paste(UNION$var_id, UNION$gene, UNION$te, sep=";")
  
  
  EMOD.T <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/emods/nominal/tumor/TUMOR.TE.GENE.nominal1.chrALL.txt.gz",header=FALSE,stringsAsFactors = FALSE,sep=" "))
  EMOD.N <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/emods/nominal/normal/NORMAL.TE.GENE.nominal1.chrALL.txt.gz",header=FALSE,stringsAsFactors = FALSE,sep=" "))
  emodNames <- c("phe_id","phe_chr","phe_from","phe_to","phe_strd","n_var_in_cis","dist_phe_var","var_id","var_chr","var_from","var_to","nom_pval","r_squared","slope","best_hit")
  colnames(EMOD.N) <- emodNames
  colnames(EMOD.T) <- emodNames
  
  
  QTL.T <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/permute/tumor/TUMOR_teqtls.chrALL.significant.txt",header=TRUE,stringsAsFactors = FALSE, sep=" "))
  QTL.N <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/permute/normal/NORMAL_teqtls.chrALL.significant.txt", header=TRUE,stringsAsFactors = FALSE,sep=" "))
  
  QTL <- rbind(QTL.T, QTL.N)
  
  
  VCF.N <- "/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_genotypes/syscol_275_normalNames_INFOfiltered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
  VCF.T <- "/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_genotypes/syscol_276_cancerNames_INFO_filtered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
  
  
  
  getExp <- function(feature, fbed){
    exp <- as.data.frame(t(fbed[which(fbed$id == feature),-c(1:6)]),stringsAsFactors = FALSE)
    exp$samples <- rownames(exp)
    colnames(exp) <- c("exp","samples")
    return(exp)
  }
  
  getDS <- function(rsid, pos, fvcf)
  {
    df <- as.data.frame(data.table::fread(cmd=paste("/Users/nikolaoslykoskoufis/Documents/Programming/Tools/bcftools-1.12/bin/bcftools view ", fvcf ," -r ", pos, " | grep -v '##'", sep=""),header=T,stringsAsFactors=F,sep="\t",colClasses=c("character")))
    df <- df[which(df$ID == rsid) ,]
    df <- df[,-c(1:9)]
    df <- as.data.frame(t(df),stringsAsFactors=F)
    df$samples <- rownames(df)
    colnames(df) <- c("genotype","samples")
    #df$GT <- sapply(strsplit(df$genotype,"\\:"),"[[",1)
    df$DS <- sapply(strsplit(df$genotype,"\\:"),"[[",2)
    #df$GT <- sub("1\\|0","0\\|1",df$GT)
    return(df[,c(2,3)])
  }
  
  getVarPos <- function(rsid, fqtl)
  {
    tmp <- fqtl[which(fqtl$var_id == rsid),]
    tmp <- tmp[1,]
    return (paste(tmp$var_chr, ":",tmp$var_from, "-",tmp$var_to,sep=""))
  }
  
  
  
  results <- foreach(rrr=1:nrow(UNION), .combine = rbind, .multicombine = T) %dopar% {
    cat(rrr,"-",nrow(UNION),"\n")
    x <- UNION[rrr,]
    
    # TE-gene correlations
    emod.n.tmp <- EMOD.N[which(EMOD.N$phe_id == x$gene & EMOD.N$var_id == x$te),]
    emod.t.tmp <- EMOD.T[which(EMOD.T$phe_id == x$gene & EMOD.T$var_id == x$te),]
    
    #gene 
    n.ds <- getDS(rsid = x$var_id, pos = getVarPos(x$var_id, QTL),fvcf = VCF.N)
    n.exp <- getExp(feature = x$gene, fbed = BED.n)
    t.ds <-  getDS(rsid = x$var_id, pos = getVarPos(x$var_id, QTL),fvcf = VCF.T)
    t.exp <- getExp(feature = x$gene, fbed = BED.t)
    
    df.n <- merge(n.ds, n.exp, by = "samples")
    df.t <- merge(t.ds, t.exp, by = "samples")
    
    lm.n.res <- summary(lm(df.n$exp~ as.numeric(as.character(df.n$DS))))
    lm.t.res <- summary(lm(df.t$exp~ as.numeric(as.character(df.t$DS))))
    
    # te
    n.exp <- getExp(feature = x$te, fbed = BED.n)
    t.exp <- getExp(feature = x$te, fbed = BED.t)
    
    df.n <- merge(n.ds, n.exp, by = "samples")
    df.t <- merge(t.ds, t.exp, by = "samples")
    
    lm.n.res.te <- summary(lm(df.n$exp~ as.numeric(as.character(df.n$DS))))
    lm.t.res.te <- summary(lm(df.t$exp~ as.numeric(as.character(df.t$DS))))
    
    d <- data.frame("teGeneSlope_n"=emod.n.tmp$slope, "teGenePval_n"=emod.n.tmp$nom_pval,
                    "teGeneSlope_t"=emod.t.tmp$slope, "teGenePval_t"=emod.t.tmp$nom_pval,
                    "geneQTL_slope_n"=lm.n.res$coefficients[2,1],"geneQTL_Pval_n"=lm.n.res$coefficients[2,4],
                    "geneQTL_slope_t"=lm.t.res$coefficients[2,1],"geneQTL_Pval_t"=lm.t.res$coefficients[2,4],
                    "teQTL_slope_n"=lm.n.res.te$coefficients[2,1],"teQTL_Pval_n"=lm.n.res.te$coefficients[2,4],
                    "teQTL_slope_t"=lm.t.res.te$coefficients[2,1],"teQTL_Pval_t"=lm.t.res.te$coefficients[2,4])
    return(cbind(x,d))
  }
  write.table(results, file="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/union/UNION_BN_ALL_correlations.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
}else
{
  results <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/union/UNION_BN_ALL_correlations.txt", header=TRUE,stringsAsFactors = FALSE,sep="\t"))
}

nrow(results[which(results$teQTL_Pval_n >= 0.05),]) # 1076
nrow(results[which(results$geneQTL_Pval_n >= 0.05),]) # 420

UNION <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/union/UNION_normal_tumor_with_posteriors.txt.gz",header=T,stringsAsFactors=F,sep="\t"))
UNION$triplet <- paste(UNION$var_id, UNION$gene, UNION$te, sep=";")
UNION$switch <- "switch"
UNION[which(UNION$normal_best == UNION$tumor_best),]$switch <- "no_switch"
UNION[which(UNION$normal_best != UNION$tumor_best & UNION$tumor_best == "m1"),]$switch <- "switch_causal"

ensembl <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_CGC/CancerDriverGenes_list.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE))

cgc <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_CGC/cancer_gene_census.csv",header=TRUE,stringsAsFactors = FALSE,sep=","))

ensembl2GeneName <- read.table("/Users/srv/beegfs/scratch/groups/data/nikos/annotations/bed/ensemblID_to_genename.txt.gz",header=F,stringsAsFactors = FALSE,sep="\t")
colnames(ensembl2GeneName) <- c("EnsemblID", "geneName")
ensembl2GeneName$EnsemblID <- sapply(strsplit(ensembl2GeneName$EnsemblID, "\\."),"[[",1)
cgc <- merge(ensembl2GeneName, cgc, by.y="Gene Symbol", by.x="geneName")


TUMOR <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/tumor/TUMOR.bnlearn_results.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
NORMAL <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/normal/NORMAL.bnlearn_results.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))

TMP <- results[results$triplet %in% TUMOR$triplet ,]
TMP <- TMP[TMP$triplet %in% UNION[which(UNION$switch == "switch_causal"),]$triplet ,]

lab.n <- c(-0.5,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1)
lab.t <- c(-0.4,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1)

TMP.cgc <- TMP[TMP$gene %in% cgc$EnsemblID,]

toRemove <- intersect(TMP.cgc[which(TMP.cgc$teGenePval_n >= 0.05),]$gene,TMP.cgc[which(TMP.cgc$teGenePval_n < 0.05),]$gene)

gg <- ggplot()+
  geom_point(data=TMP[which(TMP$teGenePval_n >= 0.05),], aes(x=teGeneSlope_n,y=teGeneSlope_t,colour="p-value>=0.05 in Normal",shape=normal_best))+
  geom_point(data=TMP[which(TMP$teGenePval_n < 0.05),], aes(x=teGeneSlope_n,y=teGeneSlope_t,colour="p-value<0.05 in Normal",shape=normal_best))+
  geom_point(data=TMP[TMP$gene %in% cgc$EnsemblID,], aes(x=teGeneSlope_n,y=teGeneSlope_t,colour="Cancer driver genes",shape=normal_best))+
  geom_point(data=TMP[TMP$gene %in% cgc[grepl("colorectal",cgc$`Tumour Types(Somatic)`),]$EnsemblID,],aes(x=teGeneSlope_n,y=teGeneSlope_t,colour="CRC genes",shape=normal_best))+
  geom_hline(yintercept=0,lty=2)+
  geom_vline(xintercept=0,lty=2)+
  labs(x="Normal TE-gene effect size (regression slope)",y="Tumor TE-gene effect size (regression slope)",title="Tumor triplets switching to causal in tumor")+
  scale_colour_manual(values=c("p-value>=0.05 in Normal"="grey","Cancer driver genes"="red4","CRC genes"="#fdbe85","p-value<0.05 in Normal"="black"))+
  scale_shape(label=c("m2"="Reactive in normal","m3"="Independent in normal"))+
  scale_x_continuous(breaks=lab.n)+ 
  scale_y_continuous(breaks=lab.t)+
  guides(shape = guide_legend(override.aes = list(size = 1)))+
  guides(color = guide_legend(override.aes = list(size = 1)))+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position=c(0.81,0.13),legend.title=element_blank(),legend.text=element_text(size=9),legend.key.size = unit(0.3,"lines"),legend.background = element_blank(), legend.box.background =  element_rect(colour="black",fill="white"),legend.spacing.y = unit(-0.1, "cm"))
ggsave(gg, filename="CAUSAL_slopes_te_gene_withCDGs.pdf",device="pdf",height=4.5,width=6)
