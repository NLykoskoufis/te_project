#!/usr/bin/env Rscript 

library(qvalue)

setwd("/Users/nikolaoslykoskoufis/Documents/PROJECTS/te_project/paper/rerun/Figures/SuppFigure4")

### LOAD DATA
colNames <- c("phe_id","phe_chr","phe_start","phe_end","phe_strand", "var_id","var_chr","var_from","var_to","nom_pval","slope")

NT <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/replicate/tumor/NORMAL_eqtls_replicated_in_TUMOR.txt.gz",header=F,stringsAsFactors=F)
colnames(NT) <- colNames
QP.NT <- qvalue(NT$nom_pval)

TN <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/replicate/normal/TUMOR_eqtls_replicated_in_NORMAL.txt.gz",header=T,stringsAsFactors=F)
colnames(TN) <- colNames
QP.TN <- qvalue(TN$nom_pval)


pdf("eQTL_pval_distribution_pi1s.pdf",8,4.5)
par(mfrow=c(1,2))
hist(NT$nom_pval, breaks=20, xlab="Adjusted P-values", ylab="Frequency", main="A. P-value distribution for eQTLs\ndiscovered in normal\ntested in tumor",ylim=c(0,11000))
legend("topright", legend=c(paste("pi1=", signif(100-QP.NT$pi0*100, 3), "%", sep=""), paste("N=",nrow(NT), sep="")), bty="n")


hist(TN$nom_pval, breaks=20, xlab="Adjusted P-values", ylab="Frequency", main="B. P-value distribution for eQTLs\ndiscovered in tumor\ntested in normal",ylim=c(0,11000))
legend("topright", legend=c(paste("pi1=", signif(100-QP.TN$pi0*100, 3), "%", sep=""), paste("N=",nrow(TN), sep="")), bty="n")
dev.off()


#### TE eqtl pval distribution 

NT.te <- NT[!grepl("ENSG",NT$phe_id),]
QP.NT.te <- qvalue(NT.te$nom_pval)

TN.te <- TN[!grepl("ENSG",TN$phe_id),]
QP.TN.te <- qvalue(TN.te$nom_pval)

pdf("TE_eqtl_pval_distribution_pi1s.pdf",8,4.5)
par(mfrow=c(1,2))
hist(NT.te$nom_pval, breaks=20, xlab="Adjusted P-values", ylab="Frequency", main="E. P-value distribution for TE-eQTLs\ndiscovered in normal\ntested in tumor",ylim=c(0,7000))
legend("topright", legend=c(paste("pi1=", signif(100-QP.NT.te$pi0*100, 3), "%", sep=""), paste("N=",nrow(NT.te), sep="")), bty="n")


hist(TN.te$nom_pval, breaks=20, xlab="Adjusted P-values", ylab="Frequency", main="F. P-value distribution for TE-eQTLs\ndiscovered in tumor\ntested in normal",ylim=c(0,7000))
legend("topright", legend=c(paste("pi1=", signif(100-QP.TN.te$pi0*100, 3), "%", sep=""), paste("N=",nrow(TN.te), sep="")), bty="n")
dev.off()

#### Gene eqtl pval distribution 

NT.gene <- NT[grepl("ENSG",NT$phe_id),]
QP.NT.gene <- qvalue(NT.gene$nom_pval)

TN.gene <- TN[grepl("ENSG",TN$phe_id),]
QP.TN.gene <- qvalue(TN.gene$nom_pval)

pdf("Gene_eqtl_pval_distribution_pi1s.pdf",8,4.5)
par(mfrow=c(1,2))
hist(NT.gene$nom_pval, breaks=20, xlab="Adjusted P-values", ylab="Frequency", main="C. P-value distribution for gene-eQTLs\ndiscovered in normal\ntested in tumor",ylim=c(0,5000))
legend("topright", legend=c(paste("pi1=", signif(100-QP.NT.gene$pi0*100, 3), "%", sep=""), paste("N=",nrow(NT.gene), sep="")), bty="n")


hist(TN.te$nom_pval, breaks=20, xlab="Adjusted P-values", ylab="Frequency", main="D. P-value distribution for gene-eQTLs\ndiscovered in tumor\ntested in normal",ylim=c(0,5000))
legend("topright", legend=c(paste("pi1=", signif(100-QP.TN.gene$pi0*100, 3), "%", sep=""), paste("N=",nrow(TN.gene), sep="")), bty="n")
dev.off()

