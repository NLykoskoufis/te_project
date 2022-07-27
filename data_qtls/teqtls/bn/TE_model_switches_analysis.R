#!/usr/bin/env Rscript 

SCRIPT_COMMENT_ = "Find how TES switch models between tumor and normal"

library(ggplot2)
library(ggpubr)
COL = RColorBrewer::brewer.pal(9,"Set3")


setwd("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/te_project/Figures/OUTLINE/C")

BN <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/teqtls/bn/tumor/TUMOR.BNlearn.ALL.extraINFO.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
BN$triplet <- paste(BN$var_id, BN$gene, BN$te, sep="_")
BN.NORMAL <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/teqtls/bn/normal/NORMAL.BNlearn.ALL.extraINFO.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
BN.NORMAL$triplet <- paste(BN.NORMAL$var_id, BN.NORMAL$gene, BN.NORMAL$te, sep="_")

UNION <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/teqtls/shared/UNION_normal_tumor_with_posteriors.txt",header=T,stringsAsFactors=F,sep="\t"))
UNION$triplet <- paste(UNION$var_id, UNION$gene, UNION$te, sep="_")
UNION$switch <- "switch"
UNION[which(UNION$normal_best == UNION$tumor_best),]$switch <- "no_switch"
UNION[which(UNION$normal_best != UNION$tumor_best & UNION$tumor_best == "m1"),]$switch <- "switch_causal"


TUMOR <- UNION[UNION$triplet %in% BN$triplet,]
NORMAL <- UNION[UNION$triplet %in% BN.NORMAL$triplet,]
CAUSAL <- TUMOR[which(TUMOR$switch == "switch_causal"),]

POS <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_expression/stat/features_genomic_positions.bed.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))

toCheck <- POS[POS$id %in% TUMOR[!TUMOR$triplet %in% CAUSAL$triplet,]$te,]
write.table(toCheck, file="~/Desktop/TEs_not_switching_to_causal_check_TCGTs.bed",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)



N = nrow(TUMOR)
C = sum(TUMOR$M1_t) / N 
R = sum(TUMOR$M2_t) / N 
I = sum(TUMOR$M3_t) / N 
TOT= rbind(as.vector(C),as.vector(R),as.vector(I))

TOT <- data.frame("models"=c("m1","m2","m3"), "model_probability"=TOT)

N = nrow(TUMOR)
C = sum(TUMOR$M1_n) / N 
R = sum(TUMOR$M2_n) / N 
I = sum(TUMOR$M3_n) / N 
TOT.N= rbind(as.vector(C),as.vector(R),as.vector(I))

TOT.N <- data.frame("models"=c("m1","m2","m3"), "model_probability"=TOT.N)
TOT.T <- TOT 
TOT.N$tissue <- "normal"
TOT.T$tissue <- "tumor"

df <- rbind(TOT.N, TOT.T)

# NORMAL TUMOR TOGETHER 


m_normal <- as.data.frame(table(TUMOR$normal_best),stringsAsFactors=F)
m_normal$prop <- m_normal$Freq / sum(m_normal$Freq)
m_tumor <- as.data.frame(table(TUMOR$tumor_best),stringsAsFactors=F)
m_tumor$prop <- m_tumor$Freq / sum(m_tumor$Freq)


m_tumor$tissue <- "tumor"
m_normal$tissue <- "normal"
mnt <- rbind(m_tumor,m_normal)

dd <- merge(df,mnt, by.x=c("models","tissue"),by.y=c("Var1","tissue"))

ggNT <- ggplot(dd,aes(x=tissue,y=model_probability,fill=factor(models),label=paste0(round(model_probability*100,0),"%\n(",Freq,")")))+
  geom_bar(stat="identity",size=1.25,color="black")+
  geom_text(position=position_stack(vjust=0.5),size=7)+
  scale_fill_manual(values=c(COL[3],COL[4],COL[6]),labels=c("m1" = "Causal (V->T->G)","m2"="Reactive (V->G->T)","m3"="Independent (V->T & V->G)"))+
  labs(y="Mean probability",x="",title="Tumor triplets tested in normal")+
  scale_x_discrete(labels=c("tumor"=paste("Tumor\n(N=",nrow(TUMOR),")",sep=""),"normal"=paste("Normal\n(N=",nrow(TUMOR),")",sep="")))+
  theme_classic(base_size=20)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="top",
        legend.direction="vertical",
        legend.justification="left",
        legend.title=element_blank())
ggsave(ggNT, filename="NORMAL.TUMOR_genome_wide_mean_probability_freq.pdf",height=8,width=8,device="pdf")












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
  scale_fill_manual(values=c("C"=COL[3],"R"=COL[4],"I"=COL[6]),labels=c("I"="Independent (V->G & V->T)","R"="Reactive (V->G->T)","C"="Causal (V->T->G)"))+
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

ggswitch <- ggplot(data=as.data.frame(table(TUMOR$switch)), aes(x=Var1, y=Freq,fill=Var1,label=Freq))+
  geom_bar(stat="identity",color="black",size=1)+
  geom_text(position=position_stack(vjust=0.5),size=8)+
  theme_classic(base_size=20)+
  ggsci::scale_fill_npg()+
  scale_x_discrete(labels=c("no_switch"="No model\nswitch","switch"="Model switch\n(other than Causal)","switch_causal"="Switch to\nCausal model"))+
  labs(x="",y="Frequency",title="Model switch for the\ntumor triplets (N=9,714)")+
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=0)
ggsave(ggswitch, filename="TUMOR_triplets_model_switches.pdf",height=7,width=7.5,device="pdf")

TUMOR$sub <- paste(TUMOR$normal_best, "->", TUMOR$tumor_best,sep="")
ns <- as.data.frame(table(TUMOR[which(TUMOR$switch == "no_switch"),]$sub))
ns$switch <- paste0("No model\nswitch (N=",sum(ns$Freq),")")
s <- as.data.frame(table(TUMOR[which(TUMOR$switch == "switch"),]$sub))
s$switch <- paste0("Switch to\nReactive/Independent model (N=",sum(s$Freq),")")
sc <- as.data.frame(table(TUMOR[which(TUMOR$switch == "switch_causal"),]$sub))
sc$switch <- paste0("Switch to\ncausal model (N=",sum(sc$Freq),")")

bb <- rbind(rbind(ns,s),sc)

x_axis <- c("m1->m1"="Causal\nto\nCausal", "m2->m2" = "Reactive\nto\nReactive", "m2->m1"="Reactive\nto\nCausal", "m3->m1"="Independent\nto\nCausal", "m1->m2"="Causal\nto\nReactive", "m3->m2"="Independent\nto\nReactive","m1->m3"="Causal\nto\nIndependent","m3->m3"="Independent\nto\nIndependent","m2->m3"="Reactive\nto\nIndependent")


gg_tumor_switches <- ggplot(data=bb, aes(x= Var1, y = Freq, fill=switch, label=Freq))+
  geom_bar(stat="identity", color="black")+
  geom_text(position=position_stack(vjust=0.5))+
  #scale_fill_manual(values=c("No model\nswitch"=npg[1], "Switch to\nReactive/Independent model"=npg[2], "Switch to\ncausal model"=npg[3]))+
  scale_fill_manual(values=c(npg[1],npg[2],npg[3]))+
  scale_x_discrete(label = x_axis)+
  labs(x="", y="Frequency", title="Model switches between normal and tumor\nfor tumor triplets (N=9,714)")+
  facet_grid(~switch,scales="free")+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position=0,
        legend.direction="vertical",
        legend.justification="left",
        legend.title=element_blank(),
        strip.text.x=element_text(face="bold", size=12))
ggsave(gg_tumor_switches, filename="tumor_triplets_model_switches_new.pdf", height=4.5, width=10)


## HOW MANY ARE CANCER DRIVER GENES ## 
cgc <- read.csv("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_CGC/cancer_gene_census.csv",header=TRUE,stringsAsFactors=FALSE,na.strings="")
ensembl<- read.table("/Users/srv/beegfs/scratch/groups/data/nikos/annotations/bed/ensemblID_to_genename.txt.gz",header=FALSE,stringsAsFactors=F)
ensembl$V1 <- sapply(strsplit(ensembl$V1, "\\."),"[[",1)
cgc <- merge(ensembl,cgc,by.x="V2",by.y="Gene.Symbol")
colnames(cgc) <- c("geneName","EnsemblID",names(cgc[,-c(1,2)]))
cgc$EnsemblID <- gsub("\\..*","",cgc$EnsemblID)


intersect(CAUSAL$gene, cgc[grepl("colorectal",cgc$Tumour.Types.Somatic.),]$EnsemblID)
# 5 genes "ENSG00000142208" "ENSG00000103126" "ENSG00000101115" "ENSG00000141510" "ENSG00000104408"
intersect(CAUSAL$gene, cgc$EnsemblID)
# 51 genes

intersect(CAUSAL$gene, cgc$EnsemblID)
cgc[cgc$EnsemblID %in% intersect(CAUSAL$gene, cgc$EnsemblID),]$geneName
paste(cgc[cgc$EnsemblID %in% intersect(CAUSAL$gene, cgc$EnsemblID),]$geneName, collapse=" ")

switchRest <- TUMOR[which(TUMOR$switch != "switch_causal"),]

library(ggvenn)

lst <- list("cgc"=cgc$EnsemblID, "Switch Causal"=CAUSAL$gene, "other switch"=switchRest$gene)
ggvenn(lst)

A <- length(intersect(cgc$EnsemblID,CAUSAL$gene))
B <- length(intersect(cgc$EnsemblID, switchRest$gene))
C <- length(setdiff(CAUSAL$gene,cgc$EnsemblID))
D <- length(setdiff(switchRest$gene, cgc$EnsemblID))
mat <- matrix(c(A,B,C,D),ncol=2, byrow=TRUE)
fisher.test(mat)


x <- as.data.frame(data.table::fread("~/Downloads/41588_2019_373_MOESM4_ESM.csv",header=TRUE,stringsAsFactors = FALSE,sep=","))
colnames(x) <- x[1,]
x <- x[-1,]

oncogene <- read.csv("~/Desktop/CompiledOncogeneOncoExaptationEvents_NatureGenPaper.csv",header=TRUE,stringsAsFactors = FALSE)
onco <- oncogene$Gene.Name

causalGenes <- ensembl[ensembl$V1 %in% unique(CAUSAL$gene) ,]$V2
intersect(onco,causalGenes)



ensembl[which(ensembl$V1 == "ENSG00000158636"),]$V2 <- "EMSY"
ensembl[which(ensembl$V1 == "ENSG00000153292"),]$V2 <- "ADGRF1"
ensembl[which(ensembl$V1 == "ENSG00000143217"),]$V2 <- "NECTIN4"
ensembl[which(ensembl$V1 == "ENSG00000102794"),]$V2 <- "ACOD1"
ensembl[which(ensembl$V1 == "ENSG00000129204"),]$V2 <- "TRE17"
ensembl[which(ensembl$V1 == "ENSG00000062650"),]$V2 <- "WAPL"
ensembl[which(ensembl$V1 == "ENSG00000158669"),]$V2 <- "GPAT4"
data.table::fwrite(ensembl, file="/Users/srv/beegfs/scratch/groups/data/nikos/annotations/bed/ensemblID_to_genename.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

intersect(onco,ensembl[ensembl$V1 %in% CAUSAL$gene,]$V2)

A <- length(intersect(onco,ensembl[ensembl$V1 %in% CAUSAL$gene,]$V2))
B <- length(intersect(onco,ensembl[ensembl$V1 %in% switchRest$gene,]$V2))
C <- length(setdiff(ensembl[ensembl$V1 %in% CAUSAL$gene,]$V2,onco))
D <- length(setdiff(ensembl[ensembl$V1 %in% switchRest$gene,]$V2,onco))
mat <- matrix(c(A,B,C,D),ncol=2,byrow=T)

fisher.test(mat)







## PATHWAY ENRICHMENT ANALYSIS ## 
gostres <- gost(query = unique(CAUSAL$gene), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
gostplot(gostres, capped = TRUE, interactive = FALSE)

PATH_RES <- gostres$result
write.table(PATH_RES,file="SWITCH_CAUSAL_gprofiler2_pathway_enrichment.txt", sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)


##################################
ANALYSIS_ = "CHECKING FOR TCGT"  #
##################################

#TCGT <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/analysis_TcTs/TCGT_PROCESSED_ALL_with_Fisher.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))

findTCGT <- function(tcgts, file){
  ndata <- data.frame()
  for(i in 1:nrow(tcgts))
  {
    r <- tcgts[i,]
    bb <- file[which(file$te == r$id),]
    b <- bb[grepl(r$tcgt_gene, bb$gene),]
    if(nrow(b) != 0){
      ndata <- rbind(ndata,cbind(b,r))
    }
  }
  return(ndata)
}

CAUSAL.TCGT <- findTCGT(TCGT[which(TCGT$fisher_pvalue < 0.05),],CAUSAL)
# 50 triplets 
# 20 unique genes 
# 43 unique TEs 



#################################################################################
ANALYSIS_ = "CHECK FOR ENRICHMENT OF THE SWITCH TO CAUSAL VS OTHERS FOR TCGTs"  #
#################################################################################

TCGT <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/analysis_TcTs/TCGT_PROCESSED_ALL_withFisher.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t")
#cgc <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_CGC/CancerDriverGenes_list.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE))
BN.TCGT <- findTCGT(TCGT[which(TCGT$padjust_fdr <= 0.05),],TUMOR)
BN.TCGT$sub <- paste(BN.TCGT$normal_best, BN.TCGT$tumor_best, sep="->")
pdf("Causal_inference_triplets_with_tcGT.pdf",height=4.5,width=5)
bp <- barplot((table(BN.TCGT$tumor_best)/nrow(BN.TCGT)) *100,names.arg = c("Causal\n(V->T->G", "Reactive\n(V->G->T)"),col=c(COL[3],COL[4]),ylab = "Triplets (%)", main="Causal inference for the\ntriplets constituted of a TcGT (N=138)")
tt <- (table(BN.TCGT$tumor_best) / nrow(BN.TCGT))*100
text(bp, tt/2, labels = paste(round(tt, digits = 2),"%"))
dev.off()

pdf("Causal_inference_triplets_with_tcGT_actual_number.pdf",height=4.5,width=5)
bp <- barplot(table(BN.TCGT$tumor_best),names.arg = c("Causal\n(V->T->G)", "Reactive\n(V->G->T)"),col=c(COL[3],COL[4]),ylab = "Triplets", main="Causal inference for the\ntriplets constituted of a TcGT (N=138)")
tt <- table(BN.TCGT$tumor_best)
text(bp, tt/2, labels = round(tt, digits = 2))
dev.off()

pdf("Model_switches_triplets_with_tcGT_actual_number.pdf",height=5.5,width=6)

npg <- pal_npg()
npg <- npg(5)

bp <- barplot(table(BN.TCGT$switch),ylab = "Triplets", main="Model switch for the\ntriplets constituted of a TcGT (N=138)\nbetween normal and tumor", col = c(npg[1],npg[2],npg[3]), names.arg = c("No model\nswitch","Model switch\n(other than causal)", "Model switch\nto causal"))
tt <- table(BN.TCGT$switch)
text(bp, tt/2, labels = round(tt, digits = 2))

dev.off()

ns <- as.data.frame(table(BN.TCGT[which(BN.TCGT$switch == "no_switch"),]$sub))
ns$switch <- paste0("No model\nswitch (N=",sum(ns$Freq),")")
s <- as.data.frame(table(BN.TCGT[which(BN.TCGT$switch == "switch"),]$sub))
s$switch <- paste0("Switch to\nReactive/Independent model (N=",sum(s$Freq),")")
sc <- as.data.frame(table(BN.TCGT[which(BN.TCGT$switch == "switch_causal"),]$sub))
sc$switch <- paste0("Switch to\ncausal model (N=",sum(sc$Freq),")")

bb <- rbind(rbind(ns,s),sc)

x_axis <- c("m1->m1"="Causal\nto\nCausal", "m2->m2" = "Reactive\nto\nReactive", "m2->m1"="Reactive\nto\nCausal", "m3->m1"="Independent\nto\nCausal", "m1->m2"="Causal\nto\nReactive", "m3->m2"="Independent\nto\nReactive")

gg_tcgt_switches <- ggplot(data=bb, aes(x= Var1, y = Freq, fill=switch, label=Freq))+
  geom_bar(stat="identity", color="black")+
  geom_text(position=position_stack(vjust=0.5))+
  #scale_fill_manual(values=c("No model\nswitch"=npg[1], "Switch to\nReactive/Independent model"=npg[2], "Switch to\ncausal model"=npg[3]))+
  scale_fill_manual(values=c(npg[1],npg[2],npg[3]))+
  scale_x_discrete(label = x_axis)+
  labs(x="", y="Frequency", title="Model switches between normal and tumor\nfor tumor triplets with tcGTs")+
  facet_grid(~switch,scales="free")+
  theme_bw()+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position=0,
        legend.direction="vertical",
        legend.justification="left",
        legend.title=element_blank(),
        strip.text.x=element_text(face="bold"))
ggsave(gg_tcgt_switches, filename="tcgt_triplets_model_switches.pdf", height=4.5, width=7)

nrow(BN.TCGT[which(BN.TCGT$tumor_best == "m1"),])

SWITCH_CAUSAL <- BN.TCGT[which(BN.TCGT$switch == "switch_causal"),]
unique(SWITCH_CAUSAL[SWITCH_CAUSAL$gene %in% cgc$EnsemblID,]$gene)

ensembl <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_CGC/CancerDriverGenes_list.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE))

natGen <- read.csv("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_CGC/All_TEderived_isoforms.csv",header=TRUE,stringsAsFactors = FALSE)
natGen <- merge(natGen, ensembl[,c(1,2)], by.x="Oncogene", by.y="geneName")





SWITCH <- BN.TCGT[which(BN.TCGT$switch == "switch"),]
NO_SWITCH <- BN.TCGT[which(BN.TCGT$switch == "no_switch"),]

ToRemove <- intersect(SWITCH$identifier, NO_SWITCH$identifier)
ToRemove <- unique(c(ToRemove, intersect(SWITCH_CAUSAL$identifier, NO_SWITCH$identifier),intersect(SWITCH_CAUSAL$identifier, SWITCH$identifier)))
OTHER <- rbind(SWITCH,NO_SWITCH)
OTHER <- OTHER[!OTHER$identifier %in% ToRemove,]


SWITCH_CAUSAL <- SWITCH_CAUSAL[!SWITCH_CAUSAL$identifier %in% ToRemove,]

#OTHER <- rbind(SWITCH,NO_SWITCH[!NO_SWITCH$identifier %in% SWITCH$identifier,])
#OTHER <- OTHER[!OTHER$identifier %in% SWITCH_CAUSAL$identifier, ]
#OTHER_ALL <- TUMOR[which(TUMOR$switch != "switch_causal"),]

getEnrichment <- function(A,B,C,D){
  a <- nrow(A)
  b <- nrow(B)
  c <- nrow(C) - a 
  d <- nrow(D) - b 
  mat <- matrix(c(a,b,c,d), byrow=TRUE,ncol=2,dimnames = list(c("Triplet TCGT","Triplet w/o TCGT"),c("A","B")))
  f <- fisher.test(mat)
  cat(" * Matrix\n")
  print(mat)
  cat("\n")
  print(f)
}

getEnrichment(A=SWITCH_CAUSAL,
              B=OTHER,
              C=TUMOR[which(TUMOR$switch == "switch_causal"),],
              D=TUMOR[which(TUMOR$switch != "switch_causal"),])
#* Matrix
#A    B
#Triplet TCGT       38   82
#Triplet w/o TCGT 2613 6981


#Fisher's Exact Test for Count Data

#data:  mat
#p-value = 0.3022
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 0.8174917 1.8454699
#sample estimates:
#odds ratio 
#  1.238086

############################################################################
ANALYSIS_ = "Differential expression of TEs for Switch CAUSAL vs Others"   #
############################################################################
TEs <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/analysis_DEG/DGE_results_DESEq2_all_TE_genes.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
TEs <- TEs[which(TEs$padj <= 0.05),]
TES <- TEs[!grepl("ENSG", TEs$id),]

OTHER <- TUMOR[which(TUMOR$switch != "switch_causal"),]
OTHER_texp <- merge(OTHER, TEs, by.x="te", by.y="id")
CAUSAL_texp <- merge(CAUSAL, TEs, by.x="te", by.y="id")
common <- intersect(CAUSAL_texp$te, OTHER_texp$te)


OTHER_texp <- OTHER_texp[!OTHER_texp$te %in% common,]
CAUSAL_texp <- CAUSAL_texp[!CAUSAL_texp$te %in% common,]
OTHER_texp$type <- "Other"
CAUSAL_texp$type <- "Switch to causal"
TUMOR_texp <- rbind(CAUSAL_texp, OTHER_texp)
median(TUMOR_texp[which(TUMOR_texp$type == "Switch to causal"),]$log2FoldChange)
median(TUMOR_texp[which(TUMOR_texp$type == "Other"),]$log2FoldChange)

ggmedian <- ggplot(data=TUMOR_texp, aes(x= type, y=log2FoldChange))+
  geom_boxplot()+
  stat_compare_means(method="wilcox.test",label.x=1.4,label.y=8)+
  #annotate("text", x=1.4, y=6, label=paste("Switch to causal\nTEs median log2FC =", median(TUMOR_texp[which(TUMOR_texp$type == "Switch to causal"),]$ratio),"\nOther TEs\nmedian log2FC =",median(TUMOR_texp[which(TUMOR_texp$type == "Other"),]$ratio), sep=""))+
  labs(x="", y="log2(tumor expression / normal expression)", title="TE expression log2 fold change difference\nbetween triplets switching to causal and other triplets")+
  theme_bw()+
  theme(panel.grid=element_blank(), 
        plot.title = element_text(hjust=0.5, face="bold"))
ggsave(ggmedian, filename="te_median_expression_difference_switch_causal_others.pdf",height=5,width=5.5,device="pdf")










###############################################################
ANALYSIS_ = "CHECK FOR TUMOR-SPECIFIC EQTLS INSIDE THE TEs"   #
###############################################################
library(lmerTest)
library(foreach)
library(doMC)
cpus <- 10
registerDoMC(cpus)

TMP <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/teqtls/enrichment/analysis_QTLenrichment/TE_variant_intersection.txt.gz",header=FALSE,stringsAsFactors = FALSE,sep="\t"))

IN_CAUSAL <- TMP[TMP$V4 %in% unique(SWITCH_CAUSAL$te),]
IN_OTHER <- TMP[TMP$V4 %in% unique(OTHER$te),]

IN_CAUSAL <- merge(SWITCH_CAUSAL, IN_CAUSAL,by.x="te",by.y="V4")
IN_OTHER <- merge(OTHER, IN_OTHER,by.x="te",by.y="V4")



BED.t <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_expression/tumor/raw/CPM_trimmed_Genes_proteinCoding_lincRNA_TEs_filtered_50_percent_and_more_276_tumor.raw.bed.gz",header=T,stringsAsFactors=F,sep="\t"))
BED.n <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_expression/normal/raw/CPM_trimmed_Genes_proteinCoding_lincRNA_TEs_filtered_50_percent_and_more_275_normal.raw.bed.gz",header=T,stringsAsFactors=F,sep="\t"))


VCF.N <- "/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V2/GENOTYPES/syscol_275_normalNames_INFOfiltered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
VCF.T <- "/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V2/GENOTYPES/syscol_276_cancerNames_INFO_filtered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"

cov.t <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_pca/tumor/PC/TUMOR.PC30_PC3_geno.pca",header=T,stringsAsFactors=F, sep="\t"))
rownames(cov.t) <- cov.t$SampleID
cov.t <- cov.t[,-c(1)] 
cov.n <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_pca/normal/PC/NORMAL.PC30_PC3_geno.pca",header=T,stringsAsFactors=F,sep="\t"))
cov.n <- cov.n[,-c(1)] 
cov.nt <- cbind(cov.t,cov.n)
COV.NT <- as.data.frame(t(cov.nt))
COV.NT <- scale(COV.NT,center=T,scale=F)

getExp <- function(feature, fbed){
  exp <- as.data.frame(t(fbed[which(fbed$gene == feature),-c(1:6)]),stringsAsFactors = FALSE)
  exp$samples <- rownames(exp)
  colnames(exp) <- c("exp","samples")
  return(exp)
}

getDS <- function(rsid, pos, fvcf)
{
  df <- as.data.frame(data.table::fread(cmd=paste("bcftools view ", fvcf ," -r ", pos, " | grep -v '##'", sep=""),header=T,stringsAsFactors=F,sep="\t",colClasses=c("character")))
  df <- df[which(df$ID == rsid) ,]
  df <- df[,-c(1:9)]
  df <- as.data.frame(t(df),stringsAsFactors=F)
  df$samples <- rownames(df)
  colnames(df) <- c("genotype","samples")
  #df$GT <- sapply(strsplit(df$genotype,"\\:"),"[[",1)
  df$DS <- sapply(strsplit(df$genotype,"\\:"),"[[",2)
  #df$GT <- sub("1\\|0","0\\|1",df$GT)
  return(df[,c(2,3)])
}

getVarPos <- function(tmp)
{
  return (paste(tmp$V7, ":",tmp$V8, "-",tmp$V9,sep=""))
}

IN_CAUSAL_RES <- foreach(rrr=1:nrow(IN_CAUSAL), .combine = rbind, .multicombine = T) %dopar% {
  
  x <- IN_CAUSAL[rrr,]
  
  n.ds <- getDS(rsid = x$V10, pos = getVarPos(x),fvcf = VCF.N)
  n.exp <- getExp(feature = x$gene, fbed = BED.n)
  
  t.ds <-  getDS(rsid = x$V10, pos = getVarPos(x),fvcf = VCF.T)
  t.exp <- getExp(feature = x$gene, fbed = BED.t)
  
  df.n <- merge(n.ds, n.exp, by = "samples")
  df.n$id <- gsub("[NRTP].*","",df.n$samples)
  df.n$tissue <- "normal"
  
  df.t <- merge(t.ds, t.exp, by = "samples")
  df.t$id <- gsub("[NRTP].*","",df.t$samples)
  df.t$tissue <- "tumor"
  
  df <- rbind(df.n, df.t)
  df$cpm_nt <- GenABEL::rntransform(df$exp)
  c.nt <- as.matrix(COV.NT[df$samples,])
  
  # Preparing variable for linear mix model 
  
  Y <- as.numeric(df$cpm_nt)
  X <- as.numeric(df$DS)
  TISSUE <- as.factor(df$tissue)
  ids <- as.factor(df$id)
  
  # Linear mix model (lmm)
  
  lmm <- try(lmer(formula = Y ~ X + TISSUE + (1|ids) + X * TISSUE + c.nt[,c(1:33)]), silent=T)
  
  if(class(lmm)=='try-error'){
    pval_eqtl_int <- NA
    beta_eqtl_int <- NA
  
  }else{
    pval_eqtl_int <- try(summary(lmm)$coefficients[37,5],silent = T)
    beta_eqtl_int <- try(summary(lmm)$coefficients[37,1],silent = T)
  }
  
  # RESIDUALIZE DATA FOR EQTL ANALYSIS # 
  
  c.n <- as.matrix(c.nt[df.n$samples,])
  df.n$resid <- GenABEL::rntransform(resid(lm(df.n$exp~ c.n[,c(1:33)])))
  
  c.t <- as.matrix(c.nt[df.t$samples,])
  df.t$resid <- GenABEL::rntransform(resid(lm(df.t$exp~ c.t[,c(1:33)])))
  
  # EQTL ANALYSIS # 
  lm.n.res <- summary(lm(df.n$resid ~ as.numeric(as.character(df.n$DS))))
  lm.t.res <- summary(lm(df.t$resid ~ as.numeric(as.character(df.t$DS))))
  
  
  dd <- cbind(x,data.frame("pval_eqtl_int"=pval_eqtl_int,"beta_eqtl_int"=beta_eqtl_int,"nom_pvalNormal"=lm.n.res$coefficients[2,4],"slopeNormal"=lm.n.res$coefficients[2,1],"nom_pvalTumor"=lm.t.res$coefficients[2,4],"slopeTumor"=lm.t.res$coefficients[2,1]))
  
}

IN_CAUSAL$tspe_hit <- IN_CAUSAL_RES$pval_eqtl_int < 0.05 & IN_CAUSAL_RES$nom_pvalTumor < 0.05 & IN_CAUSAL_RES$nom_pvalNormal > 0.05


fnomList <- Sys.glob("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/eqtls/nominal/normal/data_split/*.gz")


  
  
  
  



###############################################################
ANALYSIS_ = "PREPARING FOR CYTOSCAPE PLOTTING ALL TCGTs"      #
###############################################################
cyto <- data.frame()
SS <- merge(SWITCH_CAUSAL, ensembl, by.x="gene",by.y="V1")
for(i in 1:nrow(SS))
{
  x <- SS[i,]
  if(x$tumor_best == "m1")
  {
    df <- rbind(data.frame("A"=x$var_id,"B"=x$te,"directed"="true","id"=x$triplet,"model"="Causal"),data.frame("A"=x$te,"B"=x$V2,"directed"="true","id"=x$triplet,"model"="Causal"))
    cyto <- rbind(cyto,df)
  }
}
write.table(cyto, file="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/analysis_TcTs/SWITCH_CAUSAL_cytoscape.tab",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

SWITCH_CAUSAL[SWITCH_CAUSAL$gene %in% cgc$EnsemblID,]
TCGA <- read.table("~/Downloads/gene_set_Cancer_driver_genes_TCGA.2021-02-23.tsv",header=TRUE,stringsAsFactors = FALSE,sep="\t")
TCGA <- TCGA$id



###############################################################
ANALYSIS_ = "CHECK FOR DGE between samples with / without TCGT" #
###############################################################

BED.t <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_expression/tumor/corrected/CPM.GENES.TES_tumor.resid_notTransformed.bed.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))

getExp <- function(feature, fbed){
    exp <- as.data.frame(t(fbed[which(fbed$gene == feature),-c(1:6)]),stringsAsFactors = FALSE)
    exp$samples <- rownames(exp)
    colnames(exp) <- c("exp","samples")
    return(exp)
}
plotList <- list()
wilcox_res <- data.frame()
all_tcgt <- data.frame()
for(rrr in 1:nrow(BN.TCGT))
  {
  x <- BN.TCGT[rrr,]
  te <- x$te
  gene <- x$gene 
  sam <- unlist(strsplit(x$samples,"\\,"))
  n.samples <- sam[grepl("N",sam)]
  t.samples <- sam[!grepl("N",sam)]
  
  #te_exp <- getExp(te, BED.t)
  #colnames(te_exp) <- gsub("exp","te_exp",colnames(te_exp))
  gene_exp <- getExp(gene, BED.t)
  colnames(gene_exp) <- gsub("exp","gene_exp",colnames(gene_exp))
  
  gene_exp$tcgt <- FALSE 
  gene_exp$tcgt[gene_exp$samples %in% t.samples] <- TRUE
  gene_exp$id <- replicate(nrow(gene_exp),paste0(te,";",gene))
  all_tcgt <- rbind(all_tcgt, gene_exp)
  gg <- ggplot(data=gene_exp, aes(x=tcgt,y=gene_exp,fill=tcgt))+
    geom_boxplot()+
    geom_jitter(width=0.1)+
    stat_compare_means(method="wilcox.test")+
    scale_fill_manual(values=c("TRUE"="red4","FALSE"="lightgrey"))+
    scale_x_discrete(label=c("TRUE"="Samples\nwith TCGT","FALSE"="Samples\nwithout TCGT"))+
    labs(x="",y="normalized gene expression",title=paste(te,gene, sep="\n"))+
    theme_classic()+
    theme(legend.position=0, plot.title=element_text(hjust=0.5, face="bold"))
  
  plotList[[rrr]] <- gg
  
  test <- wilcox.test(gene_exp[gene_exp$tcgt,]$gene_exp, gene_exp[!gene_exp$tcgt,]$gene_exp)
  tcgt_median <- median(gene_exp[gene_exp$tcgt,]$gene_exp)
  notcgt_median <- median(gene_exp[!gene_exp$tcgt,]$gene_exp)
  wilcox_res <- rbind(wilcox_res, data.frame("te"=te,"gene"=gene,"triplet"=x$triplet,"tcgt_median"=tcgt_median,"notcgt_median"=notcgt_median,"wilcox_pval"=as.numeric(test$p.value)))
  
}
ml <- gridExtra::marrangeGrob(plotList, ncol=2, nrow=2)
ggsave(ml, filename="BN.TCGT_gene_dge_with_without_tcgt_tumor.pdf",device="pdf",height=8,width=8.5)


wilcox_res$fdr <- p.adjust(wilcox_res$wilcox_pval, method="fdr")
nrow(wilcox_res[which(wilcox_res$fdr < 0.05),])

toplot <- wilcox_res[which(wilcox_res$fdr < 0.05),]
toplot <- paste0(toplot$te, ";", toplot$gene)


toplot_data <- all_tcgt[all_tcgt$id %in% toplot ,]
ggplot(toplot_data, aes(x = id, y = log2(gene_exp), fill=tcgt))+
  geom_boxplot(outlier.shape=NA)+
  scale_fill_manual(values=c("FALSE"= "grey", "TRUE"="red4"))+
  theme_bw()+
  coord_flip()+
  theme(panel.grid = element_blank())

cgc_tcgt <- c("ENSG00000108375","ENSG00000157557")
cgc_tcgt_data <- all_tcgt[all_tcgt$id %in% with(wilcox_res[wilcox_res$gene %in% cgc_tcgt,], paste(te, gene, sep=";")),]

cgc_tcgt_data$id <- gsub("ENSG00000108375","RNF43",cgc_tcgt_data$id)
cgc_tcgt_data$id <- gsub("ENSG00000157557","ETS2",cgc_tcgt_data$id)
cgc_tcgt_data$id <- gsub("_.*;", "-", cgc_tcgt_data$id)
gg_cgc_tcgt <- ggplot(cgc_tcgt_data, aes(x = id, y = log2(gene_exp),fill=tcgt))+
  geom_boxplot(outlier.shape=NA)+
  scale_fill_manual(values=c("FALSE"= "grey", "TRUE"="red4"), labels=c("FALSE"="samples without tcGTs", "TRUE"="samples with tcGTs"))+
  labs(y="log2 normalized gene expression", x = "", title="tcGTs with cancer driver genes")+
  theme_bw()+
  ylim(11, 15)+
  coord_flip()+
  theme(panel.grid = element_blank(),
        plot.title=element_text(face="bold", hjust=0.5), 
        legend.position="top", legend.title = element_blank())
ggsave(gg_cgc_tcgt, filename="BN.TCGT_cgcs_expression_boxplot.pdf", height=8,width=6)





TMP <- CAUSAL.TCGT[CAUSAL.TCGT$triplet %in% wilcox_res[which(wilcox_res$fdr < 0.05),]$triplet ,]
plotList <- list()
for(rrr in 1:nrow(TMP) )
{
  x <- TMP[rrr,]
  te <- x$te
  gene <- x$gene 
  sam <- unlist(strsplit(x$samples,"\\,"))
  n.samples <- sam[grepl("N",sam)]
  t.samples <- sam[!grepl("N",sam)]
  
  #te_exp <- getExp(te, BED.t)
  #colnames(te_exp) <- gsub("exp","te_exp",colnames(te_exp))
  gene_exp <- getExp(gene, BED.t)
  colnames(gene_exp) <- gsub("exp","gene_exp",colnames(gene_exp))
  
  gene_exp$tcgt <- FALSE 
  gene_exp$tcgt[gene_exp$samples %in% t.samples] <- TRUE
  
  gg <- ggplot(data=gene_exp, aes(x=tcgt,y=gene_exp,fill=tcgt))+
    geom_boxplot()+
    geom_jitter(width=0.1)+
    stat_compare_means(method="wilcox.test",size=8,)+
    scale_fill_manual(values=c("TRUE"="red4","FALSE"="lightgrey"))+
    scale_x_discrete(label=c("TRUE"="Samples\nwith TCGT","FALSE"="Samples\nwithout TCGT"))+
    labs(x="",y="normalized gene expression",title=paste(te,gene, sep="\n"))+
    theme_classic(base_size=20)+
    theme(legend.position=0, plot.title=element_text(hjust=0.5, face="bold"))
  
  plotList[[rrr]] <- gg
}
ml <- gridExtra::marrangeGrob(plotList, ncol=2, nrow=2)
ggsave(ml, filename="CAUSAL.TCGT_gene_dge_with_without_tcgt_tumor_dgeSignificant.pdf",device="pdf",height=14,width=15)




##########################################################
ANALYSIS_ = "Get all correlations between eQTL-TE-Gene"  #
##########################################################
library(foreach)
library(doMC)
cpus <- 4
registerDoMC(cpus)

BED.t <- as.data.frame(data.table::fread("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/normal/SYSCOL_normal_275_gene_TE.FM05.resid.cpm.bed.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
BED.n <- as.data.frame(data.table::fread("/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/tumor/SYSCOL_tumor_276_gene_TE.FM05.resid.cpm.bed.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))

EMOD.T <- as.data.frame(data.table::fread("",header=TRUE,stringsAsFactors = FALSE,sep=" "))
EMOD.N <- as.data.frame(data.table::fread("",header=TRUE,stringsAsFactors = FALSE,sep=" "))

QTL.N <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/teqtls/permute/normal/NORMAL_teqtls_chrALL.significant.txt",header=TRUE,stringsAsFactors = FALSE, sep=" "))

VCF.N <- "/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V2/GENOTYPES/syscol_275_normalNames_INFOfiltered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
VCF.T <- "/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V2/GENOTYPES/syscol_276_cancerNames_INFO_filtered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"



getExp <- function(feature, fbed){
  exp <- as.data.frame(t(fbed[which(fbed$id == feature),-c(1:6)]),stringsAsFactors = FALSE)
  exp$samples <- rownames(exp)
  colnames(exp) <- c("exp","samples")
  return(exp)
}

getDS <- function(rsid, pos, fvcf)
{
  df <- as.data.frame(data.table::fread(cmd=paste("bcftools view ", fvcf ," -r ", pos, " | grep -v '##'", sep=""),header=T,stringsAsFactors=F,sep="\t",colClasses=c("character")))
  df <- df[which(df$ID == rsid) ,]
  df <- df[,-c(1:9)]
  df <- as.data.frame(t(df),stringsAsFactors=F)
  df$samples <- rownames(df)
  colnames(df) <- c("genotype","samples")
  #df$GT <- sapply(strsplit(df$genotype,"\\:"),"[[",1)
  df$DS <- sapply(strsplit(df$genotype,"\\:"),"[[",2)
  #df$GT <- sub("1\\|0","0\\|1",df$GT)
  return(df[,c(2,3)])
}

getVarPos <- function(rsid, fqtl)
{
  tmp <- fqtl[which(fqtl$var_id == rsid),]
  tmp <- tmp[1,]
  return (paste(tmp$var_chr, ":",tmp$var_from, "-",tmp$var_to,sep=""))
}

INTERSECTION <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/teqtls/shared/NORMAL_TUMOR_shared_BNlearn.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
CAUSAL <- INTERSECTION[which(INTERSECTION$normal_best != "m1" & INTERSECTION$tumor_best == "m1"),]


results <- foreach(rrr=1:nrow(TUMOR), .combine = rbind, .multicombine = T) %dopar% {
  cat(rrr,"-",nrow(TUMOR),"\n")
  x <- TUMOR[rrr,]
  
  # TE-gene correlations
  emod.n.tmp <- EMOD.N[which(EMOD.N$phe_id == x$gene & EMOD.N$var_id == x$te),]
  emod.t.tmp <- EMOD.T[which(EMOD.T$phe_id == x$gene & EMOD.T$var_id == x$te),]
  
  #gene 
  n.ds <- getDS(rsid = x$var_id, pos = getVarPos(x$var_id, QTL.N),fvcf = VCF.N)
  n.exp <- getExp(feature = x$gene, fbed = BED.n)
  t.ds <-  getDS(rsid = x$var_id, pos = getVarPos(x$var_id, QTL.N),fvcf = VCF.T)
  t.exp <- getExp(feature = x$gene, fbed = BED.t)
  
  df.n <- merge(n.ds, n.exp, by = "samples")
  df.t <- merge(t.ds, t.exp, by = "samples")
  
  lm.n.res <- summary(lm(df.n$exp~ as.numeric(as.character(df.n$DS))))
  lm.t.res <- summary(lm(df.t$exp~ as.numeric(as.character(df.t$DS))))
  
  # te
  n.exp <- getExp(feature = x$te, fbed = BED.n)
  t.exp <- getExp(feature = x$te, fbed = BED.t)
  
  df.n <- merge(n.ds, n.exp, by = "samples")
  df.t <- merge(t.ds, t.exp, by = "samples")
  
  lm.n.res.te <- summary(lm(df.n$exp~ as.numeric(as.character(df.n$DS))))
  lm.t.res.te <- summary(lm(df.t$exp~ as.numeric(as.character(df.t$DS))))
  
  d <- data.frame("teGeneSlope_n"=emod.n.tmp$slope, "teGenePval_n"=emod.n.tmp$nom_pval,
                  "teGeneSlope_t"=emod.t.tmp$slope, "teGenePval_t"=emod.t.tmp$nom_pval,
                  "geneQTL_slope_n"=lm.n.res$coefficients[2,1],"geneQTL_Pval_n"=lm.n.res$coefficients[2,4],
                  "geneQTL_slope_t"=lm.t.res$coefficients[2,1],"geneQTL_Pval_t"=lm.t.res$coefficients[2,4],
                  "teQTL_slope_n"=lm.n.res.te$coefficients[2,1],"teQTL_Pval_n"=lm.n.res.te$coefficients[2,4],
                  "teQTL_slope_t"=lm.t.res.te$coefficients[2,1],"teQTL_Pval_t"=lm.t.res.te$coefficients[2,4])
  return(cbind(x,d))
}
write.table(results, file="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/teqtls/shared/NORMAL_BN_ALL_correlations.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

nrow(results[which(results$teQTL_Pval_n >= 0.05),]) # 1076
nrow(results[which(results$geneQTL_Pval_n >= 0.05),]) # 420

TUMOR <- as.data.frame(data.table::fread("",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
NORMAL <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_qtls/teqtls/shared/NORMAL_BN_ALL_correlations.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))

un <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/union/UNION_BN_ALL_correlations.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t"))



lab.n <- c(-0.5,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1)
lab.t <- c(-0.4,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1)

results <- tumor_results[which(tumor_results$switch == "switch_causal"),]
res_cgc <- results[results$gene %in% cgc$EnsemblID,]
unique(res_cgc[which(res_cgc$teGenePval_n >= 0.05),]$gene)
unique(res_cgc[which(res_cgc$teGenePval_n < 0.05),]$gene)
intersect(res_cgc[which(res_cgc$teGenePval_n >= 0.05),]$gene,res_cgc[which(res_cgc$teGenePval_n < 0.05),]$gene)

gg <- ggplot()+
  geom_point(data=results[which(results$teGenePval_n >= 0.05),], aes(x=teGeneSlope_n,y=teGeneSlope_t,colour="p-value>=0.05 in Normal",shape=normal_best))+
  geom_point(data=results[which(results$teGenePval_n < 0.05),], aes(x=teGeneSlope_n,y=teGeneSlope_t,colour="p-value<0.05 in Normal",shape=normal_best))+
  geom_point(data=results[results$gene %in% ensembl[ensembl$V2 %in% onco,]$V1,], aes(x=teGeneSlope_n,y=teGeneSlope_t,colour="Cancer driver genes",shape=normal_best))+
  geom_point(data=results[results$gene %in% cgc[grepl("colorectal",cgc$Tumour.Types.Somatic.),]$EnsemblID,],aes(x=teGeneSlope_n,y=teGeneSlope_t,colour="CRC genes",shape=normal_best))+
  geom_hline(yintercept=0,lty=2)+
  geom_vline(xintercept=0,lty=2)+
  labs(x="Normal TE-gene effect size (regression slope)",y="Tumor TE-gene effect size (regression slope)",title="Tumor triplets switching\nto causal in tumor")+
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
ggsave(gg, filename="~/Documents/PROJECTS/te_project/Figures/OUTLINE/C/CAUSAL_slopes_te_gene_withCDGs.pdf",device="pdf",height=5,width=5.5)


df <- merge(tumor_results, BN, by=c("var_id","gene","te"))

res_cdg <- results[results$gene %in% cgc$EnsemblID,]
#table(res_cdg$normal_best)
#m2 m3 
#17 62 



###### TE Enrichment for marks ########### 

ens <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_expression/stat/TE_ensemblRegulatoryBuild_overlap.bed.gz",header=FALSE,stringsAsFactors = FALSE,sep="\t"))

df <- merge(TUMOR,ens, by.x="te", by.y="V4")
table(df[which(df$switch == "switch_causal"),]$V10) / nrow(df[which(df$switch == "switch_causal"),])
table(df[which(df$switch == "switch"),]$V10) / nrow(df[which(df$switch == "switch"),])
table(df[which(df$switch == "no_switch"),]$V10) / nrow(df[which(df$switch == "no"),])


marks <- unique(df$V10)
sc <- TUMOR[which(TUMOR$switch == "switch_causal"),]
other <- TUMOR[which(TUMOR$switch != "switch_causal"),]













