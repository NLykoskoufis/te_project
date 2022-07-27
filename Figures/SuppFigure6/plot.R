#!/usr/bin/env Rscript 
setwd("/Users/nikolaoslykoskoufis/Documents/PROJECTS/te_project/paper/rerun/Figures/SuppFigure6")
##### MAF of eQTL variants ----------------------------------

MAF <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_genotypes/syscol_variants_maf.txt.gz",header=FALSE,sep="\t", stringsAsFactors = FALSE))
colnames(MAF) <- c("var_chr","var_pos", "var_id","maf")

# normal 
n <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/permute/normal/NORMAL.PC30.All.txt.gz",header=T,stringsAsFactors = F,sep=" "))
n <- n[!is.na(n$adj_beta_pval),]
QP.n <- qvalue(n$adj_beta_pval)

n.fdr <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/permute/normal/NORMAL.PC30.All.significant.txt", header=T,stringsAsFactors = F,sep=" "))
n.fdr$type <- "TEs"
n.fdr$type[grepl("ENSG",n.fdr$phe_id)] <- "Genes"

# tumor
t <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/permute/tumor/TUMOR.PC30.All.txt.gz",header=T,stringsAsFactors = F,sep=" "))
t <- t[!is.na(t$adj_beta_pval),]
QP.t <- qvalue(t$adj_beta_pval)

t.fdr <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/permute/tumor/TUMOR.PC30.All.significant.txt", header=T,stringsAsFactors = F,sep=" "))
t.fdr$type <- "TEs"
t.fdr$type[grepl("ENSG",t.fdr$phe_id)] <- "Genes"



n <- merge(n.fdr, MAF, by="var_id")
t <- merge(t.fdr, MAF, by="var_id")

n.gene <- n[grepl("ENSG",n$phe_id),]
n.te <- n[!grepl("ENSG",n$phe_id),]
t.gene <- t[grepl("ENSG",t$phe_id),]
t.te <- t[!grepl("ENSG", t$phe_id),]

COL = brewer.pal(9,"Set1")
pdf("eQTL_frequency_spectrum.pdf",height=5,width=5)
plot(density(MAF$maf), col="grey", xlab="Minor Allele Frequency (MAF)", ylab="Density", lwd=2, main="eQTL allele frequencies")
lines(density(n.gene$maf), col=COL[1], lwd=2)
lines(density(n.te$maf), col=COL[2], lwd=2)
lines(density(t.gene$maf), col=COL[3], lwd=2)
lines(density(t.te$maf), col=COL[4], lwd=2)
legend("topright", legend=c("Tested variants", "Normal gene eQTLs", "Normal TE eQTLs", "Tumor gene eQTLs","Tumor TE eQTLs"), fill=c("grey",COL[1:4]), bty="n")
dev.off()

wilcox.test(t.gene$maf, t.te$maf)

#Wilcoxon rank sum test with continuity correction

#data:  t.gene$maf and t.te$maf
#W = 3648885, p-value = 0.1174
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(n.gene$maf, n.te$maf)

#Wilcoxon rank sum test with continuity correction

#data:  n.gene$maf and n.te$maf
#W = 27782999, p-value = 0.01913
#alternative hypothesis: true location shift is not equal to 0

median(n.gene$maf)
#[1] 0.2771739
median(n.te$maf)
#[1] 0.2844203
#> system("open .")
median(t.gene$maf)
#[1] 0.2880435
median(t.te$maf)
#[1] 0.2971014


