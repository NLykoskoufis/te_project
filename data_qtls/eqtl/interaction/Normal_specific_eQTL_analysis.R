#!/usr/bin/env Rscript


suppressMessages(library(data.table))
suppressMessages(library(foreach))
suppressMessages(library(doMC))
suppressMessages(library(lmerTest))
suppressMessages(library(GenABEL))
suppressMessages(source("/home/users/l/lykoskou/bin/tools.R"))
suppressMessages(library(ggplot2))
suppressMessages(library(qvalue))
cpus <- 10
registerDoMC(cpus)

### DATA TO IMPORT 

exp.t <- as.data.frame(data.table::fread("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/tumor/SYSCOL_tumor_276_gene_TE.FM05.raw.cpm.bed.gz",header=T,stringsAsFactors=F,sep="\t"))
exp.n <- as.data.frame(data.table::fread("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/normal/SYSCOL_normal_275_gene_TE.FM05.raw.cpm.bed.gz",header=T,stringsAsFactors=F,sep="\t"))


cov.t <- as.data.frame(data.table::fread("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_covariates/tumor/TUMOR.PC30_PC3_geno.pca",header=T,stringsAsFactors=F, sep="\t"))
rownames(cov.t) <- cov.t$SampleID
cov.t <- cov.t[,-c(1)] 
cov.n <- as.data.frame(data.table::fread("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_covariates/normal/NORMAL.PC30_PC3_geno.pca",header=T,stringsAsFactors=F,sep="\t"))
cov.n <- cov.n[,-c(1)] 
cov.nt <- cbind(cov.t,cov.n)
cov.nt <- as.data.frame(t(cov.nt))
cov.nt <- scale(cov.nt,center=T,scale=F)


eqtl <- as.data.frame(data.table::fread("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/mapping/normal/NORMAL.conditional_bestPerRank.txt.gz", header=T,stringsAsFactors=F))


vcf.t <- "/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_genotypes/syscol_276_cancerNames_INFO_filtered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
vcf.n <- "/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_genotypes/syscol_275_normalNames_INFOfiltered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"


results <- foreach(rrr=1:nrow(eqtl), .combine = rbind, .multicombine = T) %do% {

    # Reading eQTL

    x <- eqtl[rrr,]
    gene <- x$phe_id
    phe_chr <- x$phe_chr
    phe_end <- x$phe_to

    var_id <- x$var_id
    region <- paste(x$var_chr,x$var_to,sep=":")

    # Reading gene expression

    cpm.t <- exp.t[which(exp.t$id == gene & exp.t$"#chr" == phe_chr & exp.t$end == phe_end) ,]
    cpm.n <- exp.n[which(exp.n$id == gene & exp.n$"#chr" == phe_chr & exp.n$end == phe_end) ,]
    
    cpm.t <- cpm.t[,-c(1:6)]
    cpm.n <- cpm.n[,-c(1:6)]

    cpm.nt <- as.data.frame(t(cbind(cpm.t,cpm.n)))
    cpm.nt$samples <- rownames(cpm.nt)

    # Reading genotypes

    gen.t <- as.data.frame(data.table::fread(cmd=paste("bcftools view ", vcf.t ," -r ", region, " | grep -v '##'", sep=""),header=T,stringsAsFactors=F,sep="\t",colClasses=c("character")))
    gen.t <- gen.t[which(gen.t$ID == var_id) ,]
    gen.t <- gen.t[,-c(1:9)]
    gen.n <- as.data.frame(data.table::fread(cmd=paste("bcftools view ", vcf.n ," -r ", region, " | grep -v '##'", sep=""),header=T,stringsAsFactors=F,sep="\t",colClasses=c("character")))
    gen.n <- gen.n[which(gen.n$ID == var_id) ,]
    gen.n <- gen.n[,-c(1:9)]
    
    gen.nt <- as.data.frame(t(cbind(gen.t,gen.n)))
    gen.nt$samples <- rownames(gen.nt)
    
    # Process Covariates
    c.nt <- cov.nt

    # Preparing data frame 

    df <- merge(cpm.nt,gen.nt, by="samples")
    colnames(df) <- c("samples","cpm","geno")
    df$GT <- sapply(strsplit(as.character(df$geno), "\\:"),"[[",1)
    df$DS <- sapply(strsplit(as.character(df$geno), "\\:"),"[[",2)
    df$cpm_nt <- rntransform(df$cpm)
    df$ids <- gsub("[TPN].*","",df$samples)
    df$tissue <- "tumor"
    df$tissue[grepl("*N",df$samples)] <- "normal"

    c.nt <- as.matrix(c.nt[df$samples,])

    # Preparing variable for linear mix model 

    Y <- as.numeric(df$cpm_nt)
    X <- as.numeric(df$DS)
    TISSUE <- as.factor(df$tissue)
    ids <- as.factor(df$ids)

    # Linear mix model (lmm)

    lmm <- try(lmer(formula = Y ~ X + TISSUE + (1|ids) + X * TISSUE + c.nt[,c(1:33)]), silent=T)
    
    if(class(lmm)=='try-error'){
        pval_eqtl_int <- NA
        beta_eqtl_int <- NA
        beta_eqtl_normal <- NA
        beta_eqtl_tumor <- NA
    }else{
        pval_eqtl_int <- try(summary(lmm)$coefficients[37,5],silent = T)
        beta_eqtl_int <- try(summary(lmm)$coefficients[37,1],silent = T)
        beta_eqtl_normal <- try(summary(lmm)$coefficients[2,1],silent=T)
        beta_eqtl_tumor <- try(beta_eqtl_normal + beta_eqtl_int,silent=T)
    }

    return(data.frame("gene" = gene, 
                      "phe_chr" = x$phe_chr,
                      "phe_start" = x$phe_from,
                      "phe_end" = x$phe_to,
                      "var_id" = x$var_id,
                      "var_chr" = x$var_chr,
                      "var_start" = x$var_from,
                      "var_end" = x$var_to,
                      "cpm" = paste(df$cpm,collapse=","),
                      "cpm_nt" = paste(df$cpm_nt, collapse=","),
                      "geno" = paste(df$GT, collapse=","),
                      "dosage" = paste(df$DS, collapse=","),
                      "tissue" = paste(df$tissue, collapse=","),
                      "ids" = paste(df$ids, collapse =","),
                      "samples" = paste(df$samples, collapse=","),
                      "eQTL_pval" = x$bwd_pval,
                      "eqtl_int_beta" = beta_eqtl_int,
                      "eqtl_normal_beta" = beta_eqtl_normal,
                      "eqtl_tumor_beta" = beta_eqtl_tumor,
                      "eqtl_int_pval" = pval_eqtl_int
                       ))
}


results$eqtl_int_qval <- qvalue(results$eqtl_int_pval)$qvalue

pdf("diagnostic_plots_ALL_conditional_centeringCov.pdf")
hist(results$eqtl_int_pval,main="pvalue distribution")
p.qqplot(results$eqtl_int_pval)
dev.off()

write.table(results, file="normal_specific_eQTL_analysis_results_ALL_conditional_centeringCOv.txt",col.names=T,row.names=F,sep="\t",quote=F)

