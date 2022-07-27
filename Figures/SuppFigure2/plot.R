#!/usr/bin/env Rscript 

setwd("~/Documents/PROJECTS/te_project/paper/rerun/Figures/SuppFigure2")

### SuppFigure2

ALL <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_TEannotation/TE_ensemblRegBuild_miRBase_overlap_ALL_TEs.bed.gz",header=FALSE,stringsAsFactors = FALSE,sep="\t"))
colnames(ALL) <- c("phe_chr","phe_from","phe_to", "phe_id","phe_gid","phe_strd","reg_chr","reg_from","reg_to","reg_id","reg_gid","reg_strd")

ALL$id <- paste(ALL$phe_id, ALL$phe_from+1, ALL$phe_to, ALL$phe_strd, sep="_")

tmp <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/stat/features_family.bed.gz", header=TRUE,stringsAsFactors = FALSE,sep="\t"))

ALL <- ALL[!ALL$id %in% tmp$id,]

p <- table(ALL$reg_id)
p <- p[-1]
p <- p / sum(p)

lbls <- paste(names(p), "\n" , round(as.numeric(p)*100,2),"%",sep="")
pdf("regulatory_element_overlap_TE_genomeWide.pdf",7,6.5)
pie(p,col=RColorBrewer::brewer.pal(8, "Dark2"),labels=lbls, main="TEs with overlapping regulatory regions in the genome\n(N=820,981)")
dev.off()

