
library(qvalue)
library(RColorBrewer)
COL = brewer.pal(9,"Set1")
setwd("~/Documents/PROJECTS/te_project/paper/rerun/Figures/SuppFigure3/")
# normal 
n <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/permute/normal/NORMAL.PC30.All.txt.gz",header=T,stringsAsFactors = F,sep=" "))
n <- n[!is.na(n$adj_beta_pval),]
QP.n <- qvalue(n$adj_beta_pval)

n.fdr <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/permute/normal/NORMAL.PC30.All.significant.txt", header=T,stringsAsFactors = F,sep=" "))
n.fdr$type <- "TE"
n.fdr$type[grepl("ENSG",n.fdr$phe_id)] <- "GENE"

# tumor
t <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/permute/tumor/TUMOR.PC30.All.txt.gz",header=T,stringsAsFactors = F,sep=" "))
t <- t[!is.na(t$adj_beta_pval),]
QP.t <- qvalue(t$adj_beta_pval)

t.fdr <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/permute/tumor/TUMOR.PC30.All.significant.txt", header=T,stringsAsFactors = F,sep=" "))
t.fdr$type <- "TE"
t.fdr$type[grepl("ENSG",t.fdr$phe_id)] <- "GENE"



####################################################################
# PVALUE DISTRIBUTION 

pdf("eQTL_pvalue_distributions.pdf",8,4.5)
par(mfrow=c(1,2))
hist(n$adj_beta_pval, breaks=20, xlab="Adjusted P-values", ylab="Frequency", main="A. normal eQTL association signal",ylim=c(0,20000))
legend("topright", legend=c(paste("pi1=", signif(100-QP.n$pi0*100, 3), "%", sep=""), paste("#hits=", sum(QP.n$qvalue <= 0.05, na.rm=TRUE), " (5% FDR)", sep="")), bty="n")


hist(t$adj_beta_pval, breaks=20, xlab="Adjusted P-values", ylab="Frequency", main="B. tumor eQTL association signal",ylim=c(0,20000))
legend("topright", legend=c(paste("pi1=", signif(100-QP.t$pi0*100, 3), "%", sep=""), paste("#hits=", sum(QP.t$qvalue <= 0.05, na.rm=TRUE), " (5% FDR)", sep="")), bty="n")
dev.off()

