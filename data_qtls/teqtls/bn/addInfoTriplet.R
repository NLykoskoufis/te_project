#!/usr/bin/env Rscript 

suppressMessages(source("/home/users/l/lykoskou/bin/tools.R"))
suppressMessages(library(ggplot2))
suppressMessages(library(qvalue))
library(dplyr)
library(RColorBrewer)
COL = brewer.pal(9,"Set3")
library(reshape2)


# USEFUL FUNCTIONS # 

getBest <- function(x){
    posteriors <- as.numeric(x[4:6])
    max_model <- max(posteriors)
    if(posteriors[1] == max_model){return("m1")}
    if(posteriors[2] == max_model){return("m2")}
    if(posteriors[3] == max_model){return("m3")}
}

nfread <- function(file,Head=TRUE,separator="\t"){
    if(sapply(strsplit(basename(file),"\\."),tail,1) == "gz"){
        return(as.data.frame(data.table::fread(file,header=Head,stringsAsFactors=F,sep=separator)))
    }else{
        return(as.data.frame(data.table::fread(file,header=Head,stringsAsFactors=F,sep=separator)))
    }
}


addInfo <- function(qtl,te,gene,bn){
    ###################### POSITION OF TE COMPARED TO GENE ########################
    cat("Reading ",qtl,"\n")
    QTL = nfread(qtl,separator=" ")

    # PLOT2: EXTRACT TRIPLETS TO BE TESTED
    POS = unique(QTL[, c(8,10)])
    colnames(POS) = c("vid", "vpos")
    TRIPLET = paste(QTL$var_id, QTL$phe_id, sep="_")

    cat("Reading ",bn,"\n")
    #LOAD BN RESULTS 
    BN = nfread(bn)
    #BN$best <- apply(BN,1,function(x) getBest(x))

    cat("Reading ",te,"\n")
    
    #LOAD TEs 
    TES <- nfread(te)
    TES = TES[TES$id %in% BN$te,]
    TES = TES[, c(2,3,4,5,6)]
    TES$name <- sapply(strsplit(TES$gid, "\\;"),"[[",4)
    TES$tstart <- as.numeric(gsub(".*:","",sapply(strsplit(TES$gid, "\\-"),"[[",1)))
    TES$tend <- as.numeric(gsub(";.*","",sapply(strsplit(TES$gid, "\\-"),"[[",2)))
    TES <- TES[,-c(1,4,6)]
    colnames(TES) <- c("ttss","tid","tstrand","tstart","tend")
    
    cat("Reading ",gene,"\n")
    #LOAD GENES 
    GENE = nfread(gene)
    GENE = GENE[GENE$id %in% BN$gene, ]
    GENE = GENE[, c(4,5,6)]
    GENE$gstart <- gsub("R.*:","",sapply(strsplit(sapply(strsplit(GENE$gid,"\\;"),"[[",3),"-"),"[[",1))
    GENE$gend <- sapply(strsplit(sapply(strsplit(GENE$gid,"\\;"),"[[",3),"-"),"[[",2)
    GENE <- GENE[,c(1,3,4,5)]
    colnames(GENE) = c("gid","gstrand","gstart","gend")

    cat("Merging all data together\n")
    #MERGE ALL TOGETHER
    TPM = merge(BN,TES, by.x="te", by.y="tid")
    TPM = merge(TPM, GENE, by.x="gene", by.y="gid")
    TPM = merge(TPM, POS, by.x="variant",by.y="vid")
    
    cat("Computing distances between molecular phenotypes\n")
    TPM$ting = (TPM$ttss > TPM$gstart & TPM$ttss < TPM$gend)
    TPM$std = (TPM$tstrand == TPM$gstrand)
    TPM$vint <- (TPM$tstart <= TPM$vpos & TPM$tend >= TPM$vpos)

    updown <- vector()
    for(i in 1:nrow(TPM)){
        x <- TPM[i,]
        if(!x$ting){
            if(x$gstrand == "+"){
                if(x$gstart >= x$ttss){
                    updown <- c(updown, "up")
                }else{
                    updown <- c(updown,"down")
                }
            } else{
                if(x$gend >= x$ttss){
                        updown <- c(updown, "down")
                }else{
                        updown <- c(updown, "up")}
                }
                
        }else{updown <- c(updown, "in")}
    }

    TPM$updown <- updown 

    dtg <- vector()
    dgv <- vector()
    dtv <- vector()
    for(i in 1:nrow(TPM)){
        x <- TPM[i,]
        if(x$gstrand == "+"){
            gtss <- as.numeric(as.character(x$gstart))+1
            dtg <- c(dtg,abs(gtss-x$ttss))
            dgv <- c(dgv, abs(gtss-x$vpos))
            dtv <- c(dtv, abs(x$ttss-x$vpos))
        }else{
            gtss <- as.numeric(as.character(x$gend))
            dtg <- c(dtg,abs(gtss-x$ttss))
            dgv <- c(dgv, abs(gtss-x$vpos))
            dtv <- c(dtv, abs(x$ttss-x$vpos))
        }

    }
    TPM$dtg <- dtg
    TPM$dgv <- dgv
    TPM$dtv <- dtv

    #TPM$triplet <- paste(TPM$var_id,TPM$gene,TPM$te, sep="_")

    return(TPM)
}

# NORMAL 
BN.NORMAL <- "/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/normal/NORMAL.bnlearn_results.txt"
TE.NORMAL <- "/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/normal/SYSCOL_normal_275_TEs.FM05.resid.cpm.bed.gz"
GENE.NORMAL <- "/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/normal/SYSCOL_normal_275_genes.FM05.resid.cpm.bed.gz"
QTL.NORMAL <- "/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/permute/normal/NORMAL_teqtls.chrALL.significant.txt"

res <- addInfo(QTL.NORMAL,TE.NORMAL,GENE.NORMAL,BN.NORMAL)
write.table(res, file="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/normal/NORMAL.BNlearn.ALL.extraINFO.txt",sep="\t",row.names=F,col.names=T,quote=F)

# TUMOR 
BN.TUMOR <- "/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/tumor/TUMOR.bnlearn_results.txt"
TE.TUMOR <- "/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/tumor/SYSCOL_tumor_276_TEs.FM05.resid.cpm.bed.gz"
GENE.TUMOR <- "/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/tumor/SYSCOL_tumor_276_genes.FM05.resid.cpm.bed.gz"
QTL.TUMOR <- "/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/permute/tumor/TUMOR_teqtls.chrALL.significant.txt"

res.tumor <- addInfo(QTL.TUMOR,TE.TUMOR,GENE.TUMOR,BN.TUMOR)
write.table(res.tumor, file="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/tumor/TUMOR.BNlearn.ALL.extraINFO.txt",sep="\t",row.names=F,col.names=T,quote=F)
