#!/usr/bin/env Rscript 


tumor <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/mapping/tumor/TUMOR.conditional_bestPerRank.txt.gz",header=T,stringsAsFactors=F))

normal <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/mapping/normal/NORMAL.conditional_bestPerRank.txt.gz",header=T,stringsAsFactors=F,sep=" "))

COL <- RColorBrewer::brewer.pal(8, "Dark2")
col = list("TE"=COL[3], "GENE"=COL[1])


N.te <- normal[!grepl("ENSG",normal$phe_id),]
N.gene <- normal[grepl("ENSG",normal$phe_id),]



df<-merge(data.frame(table(N.gene$rank)),data.frame(table(N.te$rank)),by=1,all=T)
rownames(df) <- df$Var1
df<-df[-1]
names(df) <- c("genes","TEs")
df[,1] <- df[,1] / df[1,1]
df[,2] <- df[,2] / df[1,2]
df<-df[-1,]

pdf(file="conditional_new_colors.pdf",bg="white",8,4.5)
par(mfrow=c(1,2))
barplot(as.matrix(t(df)), main = "Independent eQTLs in normal", xlab = "Number of secondary eQTLs", ylab = "Proportion of genes/TEs with secondary eQTLs", beside = T,col=c(col$GENE, col$TE),legend.text =T)

T.te <- tumor[!grepl("ENSG",tumor$phe_id),]
T.gene <- tumor[grepl("ENSG",tumor$phe_id),]


df.t<-merge(data.frame(table(T.gene$rank)),data.frame(table(T.te$rank)),by=1,all=T)
rownames(df.t) <- df.t$Var1
df.t<-df.t[-1]
names(df.t) <- c("genes","TEs")
df.t[,1] <- df.t[,1] / df.t[1,1]
df.t[,2] <- df.t[,2] / df.t[1,2]
df.t<-df.t[-1,]

barplot(as.matrix(t(df.t)), main = "Independent eQTLs in tumor", xlab = "Number of secondary eQTLs", ylab = "Proportion of genes/TEs with secondary eQTLs", beside = T,col=c(col$GENE, col$TE),legend.text =T)
dev.off()

