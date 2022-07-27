#!/usr/bin/env Rscript 

suppressMessages(library(ggplot2))
suppressMessages(library(qvalue))
library(dplyr)
library(RColorBrewer)
COL = brewer.pal(9,"Set3")
cpus <- 20
registerDoMC(cpus)

setwd("/Users/nikolaoslykoskoufis/Documents/PROJECTS/te_project/paper/rerun/Figures/Figure5")


getBest <- function(x){
  posteriors <- as.numeric(x[4:6])
  max_model <- max(posteriors)
  if(posteriors[1] == max_model){return("m1")}
  if(posteriors[2] == max_model){return("m2")}
  if(posteriors[3] == max_model){return("m3")}
}


results_tumor <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/tumor/TUMOR.bnlearn_results.txt",header=T,stringsAsFactors=F,sep="\t")
#results_tumor$max <- pmax(results_tumor$M1, results_tumor$M2,results_tumor$M3)
#results_tumor$best <- apply(results_tumor,1,function(x) getBest(x))

m <- as.data.frame(table(results_tumor$best),stringsAsFactors=F)
m$prop <- m$Freq / sum(m$Freq)


#### MODEL PROBABILITIES GENOME WIDE 

N = nrow(results_tumor)
C = sum(results_tumor$M1) / N 
R = sum(results_tumor$M2) / N 
I = sum(results_tumor$M3) / N 
TOT= rbind(as.vector(C),as.vector(R),as.vector(I))

TOT <- data.frame("models"=c("m1","m2","m3"), "model_probability"=TOT)


## NORMAL 
results_normal <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/normal/NORMAL.bnlearn_results.txt",header=T,stringsAsFactors=F,sep="\t")
#results_normal$max <- pmax(results_normal$M1, results_normal$M2,results_normal$M3)
#results_normal$best <- apply(results_normal,1,function(x) getBest(x))

m_normal <- as.data.frame(table(results_normal$best),stringsAsFactors=F)
m_normal$prop <- m_normal$Freq / sum(m_normal$Freq)

N = nrow(results_normal)
C = sum(results_normal$M1) / N 
R = sum(results_normal$M2) / N 
I = sum(results_normal$M3) / N 
TOT.N= rbind(as.vector(C),as.vector(R),as.vector(I))

TOT.N <- data.frame("models"=c("m1","m2","m3"), "model_probability"=TOT.N)
TOT.T <- TOT 
TOT.N$tissue <- "normal"
TOT.T$tissue <- "tumor"

df <- rbind(TOT.N, TOT.T)

# NORMAL TUMOR TOGETHER 

m$tissue <- "tumor"
m_normal$tissue <- "normal"
mnt <- rbind(m,m_normal)

dd <- merge(df,mnt, by.x=c("models","tissue"),by.y=c("Var1","tissue"))

ggNT <- ggplot(dd,aes(x=tissue,y=model_probability,fill=factor(models),label=paste0(round(model_probability*100,0),"%\n(",Freq,")")))+
  geom_bar(stat="identity",size=1.25,color="black")+
  geom_text(position=position_stack(vjust=0.5),size=7)+
  scale_fill_manual(values=c(COL[3],COL[4],COL[6]),labels=c("m1" = expression("Causal (V"%->%"T"%->%"G)"),"m2"=expression("Reactive (V"%->%"G"%->%"T)"),"m3"=expression("Independent (V"%->%"T & V"%->%"G)")))+
  labs(y="Mean probability",x="",title="Genome-wide mean probability")+
  scale_x_discrete(labels=c("tumor"=paste("Tumor\n(N=",nrow(results_tumor),")",sep=""),"normal"=paste("Normal\n(N=",nrow(results_normal),")",sep="")))+
  theme_classic(base_size=20)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="top",
        legend.direction="vertical",
        legend.justification="left",
        legend.title=element_blank())
ggsave(ggNT, filename="NORMAL.TUMOR_genome_wide_mean_probability_freq.pdf",height=8,width=8,device="pdf")



####### MODEL SIGNIFICANT DIFFERENCES BETWEEN THEM ####### 

t <- m 
n <- m_normal 

fisher <- list()
for(i in 1:3){
  m <- matrix(c(n[i,]$Freq,t[i,]$Freq,
                (sum(n$Freq) - n[i,]$Freq),(sum(t$Freq) - t[i,]$Freq)),
              ncol=2,nrow=2,byrow=T)
  
  fish <- fisher.test(m)
  fish_df <- data.frame("model" = paste("m",i,sep=""),
                        "p.value" = fish$p.value,
                        "oddsratio" = fish$estimate)
  fisher[[i]] <- fish_df
}
r <- Reduce(rbind,fisher)

#model       p.value oddsratio
#odds ratio     m1  0.000000e+00 0.3217265
#odds ratio1    m2  7.181141e-86 1.7198223
#odds ratio2    m3 1.344678e-126 2.6506424

#### MODEL SUBSTITUTIONS FOR THE TUMOR TRIPLETS ####

COL = RColorBrewer::brewer.pal(9,"Set3")


BN <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/tumor/TUMOR.BNlearn.ALL.extraINFO.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
#BN$triplet <- paste(BN$var_id, BN$gene, BN$te, sep="_")
BN.NORMAL <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/normal/NORMAL.BNlearn.ALL.extraINFO.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
#BN.NORMAL$triplet <- paste(BN.NORMAL$var_id, BN.NORMAL$gene, BN.NORMAL$te, sep="_")

UNION <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/union/UNION_normal_tumor_with_posteriors.txt.gz",header=T,stringsAsFactors=F,sep="\t"))
UNION$triplet <- paste(UNION$var_id, UNION$gene, UNION$te, sep=";")
UNION$switch <- "switch"
UNION[which(UNION$normal_best == UNION$tumor_best),]$switch <- "no_switch"
UNION[which(UNION$normal_best != UNION$tumor_best & UNION$tumor_best == "m1"),]$switch <- "switch_causal"


TUMOR <- UNION[UNION$triplet %in% BN$triplet,]
NORMAL <- UNION[UNION$triplet %in% BN.NORMAL$triplet,]
CAUSAL <- TUMOR[which(TUMOR$switch == "switch_causal"),]

df <- TUMOR

df[which(df$normal_best == "m1"),]$normal_best <- "C"
df[which(df$normal_best == "m2"),]$normal_best <- "R"
df[which(df$normal_best == "m3"),]$normal_best <- "I"


df[which(df$tumor_best == "m1"),]$tumor_best <- "C"
df[which(df$tumor_best == "m2"),]$tumor_best <- "R"
df[which(df$tumor_best == "m3"),]$tumor_best <- "I"

df$key <- paste(df$normal_best, "_", df$tumor_best, sep="")
dif <- as.data.frame(table(df$key),stringsAsFactors=F)

dif$normal <- sapply(strsplit(dif$Var1,"\\_"),"[[",1)
dif$tumor <- sapply(strsplit(dif$Var1,"\\_"),"[[",2)
dif <- dif[,-1]
colnames(dif) <- c("freq","key","model")
dif$tissue <- "tumor"

dn <- as.data.frame(table(df$normal_best),stringsAsFactors=F)
colnames(dn) <- c("key","freq")
dn <- dn[,c(2,1)]
dn$model <- dn$key
dn$tissue <- "normal"

d <- rbind(dn,dif)
d$tissue <- factor(d$tissue, levels=c("normal","tumor"))
supp_labs <- c("C"="Causal model\nchanges","R"="Reactive model\nchanges","I"="Independent model\nchanges")
ggbaby2 <- ggplot(data=d,aes(x=tissue, y=freq,fill=model,label=freq))+
  geom_bar(stat="identity",colour="white",size=1,width = 0.9)+
  geom_text(size=5,position=position_stack(vjust=0.5))+
  facet_wrap(~key,scales="free",labeller = labeller(key = supp_labs))+
  scale_fill_manual(values=c(COL[3],COL[4],COL[6]),labels=c("C" = expression("Causal (V"%->%"T"%->%"G)"),"R"=expression("Reactive (V"%->%"G"%->%"T)"),"I"=expression("Independent (V"%->%"T & V"%->%"G)")))+
  labs(x="",y="frequency",title=paste("Model substitutions for the significant tumor triplets\n(N=",nrow(df),")",sep=""))+
  theme_bw(base_size=20)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="top",
        legend.direction="vertical",
        legend.justification="left",
        legend.title=element_blank())
ggsave(ggbaby2, filename="model_substitutions_for_tumor_triplets_ggbaby.pdf",height=8,width=10)


#### TUMOR TRIPLETS SWITCHING TO CAUSAL #### 


ggswitch <- ggplot(data=as.data.frame(table(TUMOR$switch)), aes(x=Var1, y=Freq,fill=Var1,label=Freq))+
  geom_bar(stat="identity",color="black",size=1)+
  geom_text(position=position_stack(vjust=0.5),size=8)+
  theme_classic(base_size=20)+
  ggsci::scale_fill_npg()+
  scale_x_discrete(labels=c("no_switch"="No model\nswitch","switch"="Model switch\n(other than Causal)","switch_causal"="Switch to\nCausal model"))+
  labs(x="",y="Frequency",title="Model switch for the\ntumor triplets (N=9,528)")+
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=0)
ggsave(ggswitch, filename="TUMOR_triplets_model_switches.pdf",height=7,width=7.5,device="pdf")


#### TE effect on gene expression + model switches ####
## HOW MANY ARE CANCER DRIVER GENES ## 
oncogenes <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_CGC/CompiledOncogeneOncoExaptationList_NatGen2019.csv",header=TRUE,sep=",",stringsAsFactors = FALSE))
ensembl <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_CGC/CancerDriverGenes_list.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE))
oncogenes <- merge(ensembl[,c(1,2)], oncogenes, by.x="geneName", by.y="Gene Name")

cgc <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_CGC/cancer_gene_census.csv",header=TRUE,stringsAsFactors = FALSE,sep=","))

ensembl2GeneName <- read.table("/Users/srv/beegfs/scratch/groups/data/nikos/annotations/bed/ensemblID_to_genename.txt.gz",header=F,stringsAsFactors = FALSE,sep="\t")
colnames(ensembl2GeneName) <- c("EnsemblID", "geneName")
ensembl2GeneName$EnsemblID <- sapply(strsplit(ensembl2GeneName$EnsemblID, "\\."),"[[",1)
cgc <- merge(ensembl2GeneName, cgc, by.y="Gene Symbol", by.x="geneName")




intersect(CAUSAL$gene, cgc[grepl("colorectal",cgc$`Tumour Types(Somatic)`),]$EnsemblID)
# 5 genes "ENSG00000142208" "ENSG00000103126" "ENSG00000101115" "ENSG00000141510" "ENSG00000104408"
intersect(CAUSAL$gene, cgc$EnsemblID)
# 51 genes with CGC and 
intersect(CAUSAL$gene, oncogenes$EnsemblID)
#56 genes 

# 11 genes in common between CDG and oncogenes


#### ENRICHMENT FOR CGC SWITCH CAUSAL VS REST #####
REST <- TUMOR[which(TUMOR$switch != "switch_causal"),]
REST <- REST[!REST$gene %in% CAUSAL$gene,]

A <- length(intersect(CAUSAL$gene, cgc$EnsemblID))
B <- length(intersect(REST$gene, cgc$EnsemblID))
C <- length(setdiff(CAUSAL$gene, cgc$EnsemblID))
D <- length(setdiff(REST$gene, cgc$EnsemblID))

mat <- matrix(c(A,B,C,D), ncol=2, byrow=TRUE)
fisher.test(mat)
#### p-value 0.2903 with CDG thus no significant enrichment of switch to causal triplets for driver cancer genes. 

REST <- TUMOR[which(TUMOR$switch != "switch_causal"),]
REST <- REST[!REST$gene %in% CAUSAL$gene,]

A <- length(intersect(CAUSAL$gene, oncogenes$EnsemblID))
B <- length(intersect(REST$gene, oncogenes$EnsemblID))
C <- length(setdiff(CAUSAL$gene, oncogenes$EnsemblID))
D <- length(setdiff(REST$gene, oncogenes$EnsemblID))

mat <- matrix(c(A,B,C,D), ncol=2, byrow=TRUE)
fisher.test(mat)

##### ENRICHMENT FOR CAUSAL VS ALL THE REST ######## 
TUMOR_CAUSAL <- TUMOR[which(TUMOR$tumor_best == "m1"),]
TUMOR_REST <-   TUMOR[which(TUMOR$tumor_best != "m1"),]

common <- intersect(TUMOR_CAUSAL$gene, TUMOR_REST$gene)

TUMOR_CAUSAL <- TUMOR_CAUSAL[!TUMOR_CAUSAL$gene %in% common ,]
TUMOR_REST <- TUMOR_REST[!TUMOR_REST$gene %in% common ,]

A <- length(intersect(TUMOR_CAUSAL$gene, cgc$EnsemblID))
B <- length(intersect(TUMOR_REST$gene, cgc$EnsemblID))
C <- length(setdiff(TUMOR_CAUSAL$gene, cgc$EnsemblID))
D <- length(setdiff(TUMOR_REST$gene, cgc$EnsemblID))

mat <- matrix(c(A,B,C,D), ncol=2, byrow=TRUE)
fisher.test(mat)


#### FIGURE 5D

results <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/union/UNION_BN_ALL_correlations.txt.gz", header=TRUE,stringsAsFactors = FALSE,sep="\t"))


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
ggsave(gg, filename="CAUSAL_slopes_te_gene_withCDGs.pdf",device="pdf",height=5,width=6.5)

cgc[cgc$EnsemblID %in% TMP[TMP$gene %in% cgc$EnsemblID,]$gene ,]$geneName
#[1] "ACSL6"    "AKT1"     "ATIC"     "CASP8"    "CDKN1B"   "CRTC3"    "DDX5"     "DNM2"     "EIF3E"    "ELF3"     "ERC1"     "ERCC3"
#[13] "FANCA"    "FGFR1OP"  "FGFR4"    "FLCN"     "HLA-A"    "KNSTRN"   "LATS1"    "LPP"      "MAPK1"    "MYC"      "NCKIPSD"  "NCOR1"   
#[25] "NIN"      "NOTCH1"   "NTHL1"    "NUTM2D"   "PABPC1"   "PAX8"     "PBRM1"    "PDE4DIP"  "PMS2"     "POLG"     "PTCH1"    "PTEN"
#[37] "RAF1"     "RFWD3"    "RNF43"    "SBDS"     "SDHD"     "SETD1B"   "SLC45A3"  "SS18L1"   "TNFRSF14" "TNFRSF17" "TP53"     "TRAF7"   
#[49] "TRIP11"   "TSC2"     "USP8"     "VHL"      "XPC"      "ZMYM2"    "ZNF384"


cgc[cgc$EnsemblID %in% TMP[TMP$gene %in% cgc[grepl("colorectal",cgc$`Tumour Types(Somatic)`),]$EnsemblID,]$gene,]$geneName
#"AKT1"  "EIF3E" "TP53" 

###### CHECK FOR CDG ENRICHMENT IN CAUSAL TRIPLET VS OTHER TRIPLETS ######


SWC <- UNION[which(UNION$switch == "switch_causal" & UNION$triplet %in% TUMOR$triplet),]
OTH <- UNION[which(UNION$switch != "switch_causal" & UNION$triplet %in% TUMOR$triplet),]

cdg_in_swc <- length(unique(SWC[SWC$gene %in% cgc$EnsemblID,]$gene))
cdg_in_oth <- length(unique(OTH[OTH$gene %in% cgc$EnsemblID,]$gene))
total_swc_genes <- length(unique(SWC$gene))
total_oth_genes <- length(unique(OTH$gene))

mat <- matrix(c(cdg_in_swc,cdg_in_oth,(total_swc_genes - cdg_in_swc),(total_oth_genes - cdg_in_oth)),ncol=2, byrow=T)
mat
fisher.test(mat)
