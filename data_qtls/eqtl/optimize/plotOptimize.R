#!/usr/bin/env Rscript 

setwd("~/Documents/PROJECTS/te_project/paper/rerun/data_qtl/")

NORMAL <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/optimize/normal/NORMAL.Neqtls.PC.txt",header=FALSE,stringsAsFactors = FALSE,sep="\t")
TUMOR <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/optimize/tumor/TUMOR.Neqtls.PC.txt",header=FALSE,stringsAsFactors = FALSE,sep="\t")

NORMAL <- NORMAL[-1,]
NORMAL$V1 <- basename(NORMAL$V1)
NORMAL$V1 <- as.numeric(as.character(NORMAL$V1))

TUMOR <- TUMOR[-1,]
TUMOR$V1 <- basename(TUMOR$V1)
TUMOR$V1 <- as.numeric(as.character(TUMOR$V1))

pdf("analysisA_n_eqtls_per_pcs.pdf",8,4)
LPC = c(0,1,2,5,10,20,30,40,50,60,70,80,90,100)
nLPC = length(LPC)
par(mfrow=c(1,2))
plot(LPC,NORMAL$V2, type="b", pch=20, xlab="Number of PCs used as covariates", ylab="Number of eQTLs at 5% FDR", main="A.    eQTL discovery in normal", ylim=c(0,17000))
abline(v=30, col="grey")
points(30, NORMAL[7,2], pch=19, col="red", cex=1.5)
text(30,16500, labels=paste("n=",NORMAL[7,2],sep=""), col="red")


plot(LPC,TUMOR$V2, type="b", pch=20, xlab="Number of PCs used as covariates", ylab="Number of eQTLs at 5% FDR", main="B.    eQTL discovery in tumor", ylim=c(0,17000))
abline(v=30, col="grey")
points(30, TUMOR[7,2], pch=19, col="red", cex=1.5)
text(30,7500, labels=paste("n=",TUMOR[7,2],sep=""), col="red")
dev.off()

