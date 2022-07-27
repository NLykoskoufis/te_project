
suppressMessages(library(data.table))
suppressMessages(library(foreach))
suppressMessages(library(doMC))
suppressMessages(library(lmerTest))
suppressMessages(library(GenABEL))
suppressMessages(library(ggplot2))
suppressMessages(library(qvalue))
library(dplyr)
library(RColorBrewer)
COL = brewer.pal(9,"Set3")
library(bnlearn)





# LOAD FILES # 

BN.n <- as.data.frame(data.table::fread("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/normal/NORMAL.bnlearn_results.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
BN.t <- as.data.frame(data.table::fread("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/tumor/TUMOR.bnlearn_results.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))

UNION <- union(BN.n$triplet, BN.t$triplet)

bed.n <- as.data.frame(data.table::fread("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/normal/SYSCOL_normal_275_gene_TE.FM05.resid.cpm.bed.gz",header=T,stringsAsFactors=F,sep="\t"))
bed.t <- as.data.frame(data.table::fread("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/tumor/SYSCOL_tumor_276_gene_TE.FM05.resid.cpm.bed.gz",header=T,stringsAsFactors=F,sep="\t"))

QTL.n <- as.data.frame(data.table::fread("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/permute/normal/NORMAL_teqtls.chrALL.significant.txt",header=T,stringsAsFactors=F,sep=" "))
QTL.t <- as.data.frame(data.table::fread("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/permute/tumor/TUMOR_teqtls.chrALL.significant.txt",header=T,stringsAsFactors=F,sep=" "))
QTL <- rbind(QTL.n, QTL.t)

vcf.t <- "/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_genotypes/syscol_276_cancerNames_INFO_filtered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
vcf.n <- "/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_genotypes/syscol_275_normalNames_INFOfiltered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"



# LOAD FUNCTIONS # 

getBN <- function(df){
    df$V = (df$V - mean(df$V))/sd(df$V)
    #In statistics, @@standardization@@ is the process of putting different variables on the same scale. This process allows you to compare scores between different types of variables. Typically, to standardize variables, you calculate the mean and standard deviation for a variable. Then, for each observed value of the variable, you subtract the mean and divide by the standard deviation.

    # V = variant
    # T = Transposable ELEMENT
    # G = Gene

    #list of 4 models to test
    m=rep("", 3)
    m[1]="[V][T|V][G|T]"	#V -> T -> G # Scenario where the variants affects the TE that affects the G
    m[2]="[V][G|V][T|G]"	#V -> G -> T # Scenario where the variant affects the Gene and the Gene affects the TE
    m[3]="[V][T|V][G|V]"	#V -> G, V-> T # Scenario where the variants affects independently the TE and Gene.

    #calcuate P(D|G) where G=1,2,3
    loglik=rep(0, 3)
    for (i in 1:3) {
        net = model2network(m[i])
        loglik[i]=score(net, df, type="bge")
    }

    #assume constant prior over the 4 network configurations, i.e. P(G) = 0.25, and compute posterior
    prior=rep(1/3, 3)
    posteriors = exp(loglik - max(loglik)) * prior # scores are in log scale. substraction means division!!!
    posteriors = posteriors / sum(posteriors) # Division by the sum so that all of the probabilities add to 1!!

    #write output
    cat ("M1=", signif(posteriors[1],3), " M2=", signif(posteriors[2],3), "M3=", signif(posteriors[3],3), "\n")
    return(posteriors)
}



getdT <- function(var_id,gene,te){
    # Get Gene Expression 
    exp_gene <- bed.t[which(bed.t$id == gene),]
    exp_gene <- as.data.frame(t(exp_gene[,-c(1:6)]),stringsAsFactors=F)
    exp_gene$SampleID <- rownames(exp_gene)
    colnames(exp_gene) <- c("cpm_gene","SampleID")

    # Get TE expression 
    exp_te <- bed.t[which(bed.t$id == te),]
    exp_te <- as.data.frame(t(exp_te[,-c(1:6)]),stringsAsFactors=F)
    exp_te$SampleID <- rownames(exp_te)
    colnames(exp_te) <- c("cpm_te","SampleID")

    # Get Dosage 
    var <- QTL[which(QTL$var_id == var_id),]
    var <- as.data.frame(var[1,])
    var_pos <- paste(var$var_chr, ":" ,var$var_from, "-", var$var_to,sep="")
    var <- as.data.frame(data.table::fread(cmd=paste("bcftools view ", vcf.t ," -r ", var_pos, " | grep -v '##'", sep=""),header=T,stringsAsFactors=F,sep="\t",colClasses=c("character")))
    var <- var[which(var$ID == var_id) ,]
    if(var$ID != var_id){cat("There is an issue! NOt same variant extracted for triplet ",tri,"\n"); exit()}
    var <- var[,-c(1:9)]
    var <- as.data.frame(t(var),stringsAsFactors=F)
    colnames(var) <- c("genotype")
    var$SampleID <- rownames(var)
    var$DS <- sapply(strsplit(var$genotype,"\\:"),"[[",2)

    d <- merge(merge(var,exp_gene,by="SampleID"),exp_te,by="SampleID")
    d <- d[,c("DS","cpm_gene","cpm_te")]
    colnames(d) <- c("V","G","T")
    d$V <- as.numeric(as.character(d$V))
    d$G <- as.numeric(d$G)
    d$T <- as.numeric(d$T)
    return(d)
}


getdN <- function(var_id,gene,te){
    # Get Gene Expression 
    exp_gene <- bed.n[which(bed.n$id == gene),]
    exp_gene <- as.data.frame(t(exp_gene[,-c(1:6)]),stringsAsFactors=F)
    exp_gene$SampleID <- rownames(exp_gene)
    colnames(exp_gene) <- c("cpm_gene","SampleID")

    # Get TE expression 
    exp_te <- bed.n[which(bed.n$id == te),]
    exp_te <- as.data.frame(t(exp_te[,-c(1:6)]),stringsAsFactors=F)
    exp_te$SampleID <- rownames(exp_te)
    colnames(exp_te) <- c("cpm_te","SampleID")

    # Get Dosage 
    var <- QTL[which(QTL$var_id == var_id),]
    var <- as.data.frame(var[1,])
    var_pos <- paste(var$var_chr, ":" ,var$var_from, "-", var$var_to,sep="")
    var <- as.data.frame(data.table::fread(cmd=paste("bcftools view ", vcf.n ," -r ", var_pos, " | grep -v '##'", sep=""),header=T,stringsAsFactors=F,sep="\t",colClasses=c("character")))
    var <- var[which(var$ID == var_id) ,]
    if(var$ID != var_id){cat("There is an issue! NOt same variant extracted for triplet ",tri,"\n"); exit()}
    var <- var[,-c(1:9)]
    var <- as.data.frame(t(var),stringsAsFactors=F)
    colnames(var) <- c("genotype")
    var$SampleID <- rownames(var)
    var$DS <- sapply(strsplit(var$genotype,"\\:"),"[[",2)

    d <- merge(merge(var,exp_gene,by="SampleID"),exp_te,by="SampleID")
    d <- d[,c("DS","cpm_gene","cpm_te")]
    colnames(d) <- c("V","G","T")
    d$V <- as.numeric(as.character(d$V))
    d$G <- as.numeric(d$G)
    d$T <- as.numeric(d$T)
    return(d)
}



getBest <- function(posteriors){
    max_model <- max(posteriors)
    if(posteriors[1] == max_model){return("m1")}
    if(posteriors[2] == max_model){return("m2")}
    if(posteriors[3] == max_model){return("m3")}
}



res <- data.frame()
### START OF LOOP ### 
for(i in 1:length(UNION)){
    cat(i,"-",length(UNION),"\n")
    x <- UNION[i]
    
    var_id <- sapply(strsplit(x,"\\;"),"[[",1) 
    gene <- sapply(strsplit(x,"\\;"),"[[",2)
    te <- sapply(strsplit(x,"\\;"),"[[",3) 

    dn <- getdN(var_id, gene, te)
    dt <- getdT(var_id, gene, te)

    bnN <- getBN(dn)    
    bnT <- getBN(dt)
    r <- data.frame("var_id"=var_id, "gene"=gene,"te"=te,"M1_n"=bnN[1],"M2_n"=bnN[2],"M3_n"=bnN[3],"M1_t"=bnT[1],"M2_t"=bnT[2],"M3_t"=bnT[3])
    r$normal_best <- getBest(bnN)
    r$tumor_best <- getBest(bnT)
    res <- rbind(res,r)
}

write.table(res, file="UNION_normal_tumor_with_posteriors.txt",sep="\t",col.names=T,row.names=F,quote=F)


#UNION <- as.data.frame(data.table::fread("zcat < /Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/teqtls/shared/NORMAL.TUMOR_ALL_union.txt.gz",header=T,stringsAsFactors=F,sep="\t"))
#UNION$key <- paste(UNION$normal_best, UNION$tumor_best,sep="_")

#res <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/teqtls/shared/UNION_normal_tumor_with_posteriors.txt",header=T,stringsAsFactors=F,sep="\t"))
#res$key <- paste(res$normal_best, res$tumor_best,sep="_")

#res$triplet <- paste(res$var_id, res$gene, res$te,sep="_")
#res <- merge(res,UNION,by="triplet")




#N <- nrow(res)
#Cn <- sum(res$M1_n) / N 
#Rn <- sum(res$M2_n) / N
#In <- sum(res$M3_n) / N

#Ct <- sum(res$M1_t) / N
#Rt <- sum(res$M2_t) / N
#It <- sum(res$M3_t) / N


#TOTn= rbind(as.vector(Cn),as.vector(Rn),as.vector(In))
#TOTt= rbind(as.vector(Ct),as.vector(Rt),as.vector(It))

#TOTn <- data.frame("models"=c("m1","m2","m3"), "model_probability"=TOTn)
#TOTt <- data.frame("models"=c("m1","m2","m3"), "model_probability"=TOTt)
#TOTn$tissue <- "normal"
#TOTt$tissue <- "tumor"
#TOT <- rbind(TOTn,TOTt)


#ggunion<-ggplot(TOT,aes(x=tissue,y=model_probability,fill=factor(models),label=paste0(round(model_probability*100,0),"%")))+
#    geom_bar(stat="identity",size=1.25,color="black")+
#    geom_text(position=position_stack(vjust=0.5),size=5)+
#    scale_fill_manual(values=c(COL[3],COL[4],COL[6]),labels=c("m1" = "Causal (V->T->G)","m2"="Reactive (V->G->T)","m3"="Independent (V->T & V->G)"))+
#    #scale_x_discrete(labels=c("normal"=paste("normal\n(N=",sum(m_normal$Freq),")",sep=""),"tumor"=paste("tumor\n(N=",sum(m$Freq),")",sep="")))+
#    labs(y="Model probability",x="",title="Union of normal and tumor triplets (N=21647)\ngenome-wide mean probability")+
#    theme_classic(base_size=20)+
#    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
#        plot.title=element_text(hjust=0.5,face="bold"),
#        legend.position="top",
#        legend.direction="vertical",
#        legend.justification="left",
#        legend.title=element_blank())

#ggsave(ggunion, filename="UNION_genome_wide_probabilities.pdf",device="pdf",height=8,width=8)

#df <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/teqtls/shared/UNION_normal_tumor_with_posteriors.txt",header=T,stringsAsFactors=F,sep="\t"))

#df[which(df$normal_best == "m1"),]$normal_best <- "C"
#df[which(df$normal_best == "m2"),]$normal_best <- "R"
#df[which(df$normal_best == "m3"),]$normal_best <- "I"


#df[which(df$tumor_best == "m1"),]$tumor_best <- "C"
#df[which(df$tumor_best == "m2"),]$tumor_best <- "R"
#df[which(df$tumor_best == "m3"),]$tumor_best <- "I"

#df$key <- paste(df$normal_best, "_", df$tumor_best, sep="")
#dif <- as.data.frame(table(df$key),stringsAsFactors=F)

#dif$normal <- sapply(strsplit(dif$Var1,"\\_"),"[[",1)
#dif$tumor <- sapply(strsplit(dif$Var1,"\\_"),"[[",2)
#dif <- dif[,-1]
#colnames(dif) <- c("freq","key","model")
#dif$tissue <- "tumor"


#dn <- as.data.frame(table(df$normal_best),stringsAsFactors=F)
#colnames(dn) <- c("key","freq")
#dn <- dn[,c(2,1)]
#dn$model <- dn$key
#dn$tissue <- "normal"

#d <- rbind(dn,dif)
#d$tissue <- factor(d$tissue, levels=c("normal","tumor"))
#supp_labs <- c("C"="Causal model\nchanges","R"="Reactive model\nchanges","I"="Independent model\nchanges")
#ggbaby2 <- ggplot(data=d,aes(x=tissue, y=freq,fill=model,label=freq))+
#  geom_bar(stat="identity",colour="white",size=1,width = 0.9)+
#  geom_text(size=5,position=position_stack(vjust=0.5))+
#  facet_wrap(~key,scales="free",labeller = labeller(key = supp_labs))+
#  scale_fill_manual(values=c("C"=COL[3],"R"=COL[4],"I"=COL[6]),labels=c("I"="Independent (V->G & V->T)","R"="Reactive (V->G->T)","C"="Causal (V->T->G)"))+
#  labs(x="",y="frequency",title=paste("Model substitutions for the union between\nnormal and tumor triplets\n(N=",nrow(df),")",sep=""))+
#  theme_bw(base_size=20)+
#  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),strip.background = element_blank(),
#        panel.border = element_rect(colour = "black"),
#        plot.title=element_text(hjust=0.5,face="bold"),
#        legend.position="top",
#        legend.direction="vertical",
#        legend.justification="left",
#        legend.title=element_blank())

#ggsave(ggbaby2, filename="model_substitutions_normal_tumor_UNION.pdf",device="pdf",height=8,width=10)


