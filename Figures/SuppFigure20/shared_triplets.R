#!/usr/bin/env Rscript 


library(ggplot2)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
COL = brewer.pal(9,"Set3")
setwd("/Users/nikolaoslykoskoufis/Documents/PROJECTS/te_project/paper/rerun/Figures/SuppFigure20")

df <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/shared/NORMAL_TUMOR_shared_BNlearn.txt.gz",header=T,stringsAsFactors=F,sep="\t")


df[which(df$normal_best == "m1"),]$normal_best <- "C"
df[which(df$normal_best == "m2"),]$normal_best <- "R"
df[which(df$normal_best == "m3"),]$normal_best <- "I"


df[which(df$tumor_best == "m1"),]$tumor_best <- "C"
df[which(df$tumor_best == "m2"),]$tumor_best <- "R"
df[which(df$tumor_best == "m3"),]$tumor_best <- "I"

df$key <- paste(df$normal_best, " to ", df$tumor_best, sep="")

### GENOME-WIDE MODEL PERCENTAGES ### 


N = nrow(df)
CN = sum(df$normal_m1) / N 
RN = sum(df$normal_m2) / N 
IN= sum(df$normal_m3) / N 

CT= sum(df$tumor_m1) / N 
RT= sum(df$tumor_m2) / N 
IT= sum(df$tumor_m3) / N 
TOTN= rbind(as.vector(CN),as.vector(RN),as.vector(IN))
TOTT=rbind(as.vector(CT),as.vector(RT),as.vector(IT))
TOTN <- data.frame("models"=c("m1","m2","m3"), "model_probability"=TOTN,"tissue"="normal")
TOTT <- data.frame("models"=c("m1","m2","m3"), "model_probability"=TOTT,"tissue"="tumor")

dd <- rbind(TOTN,TOTT)


gg_dd <- ggplot(dd,aes(x=tissue,y=model_probability,fill=factor(models),label=paste0(round(model_probability*100,0),"%")))+
    geom_bar(stat="identity",size=1.25,color="black")+
    geom_text(position=position_stack(vjust=0.5),size=5)+
    scale_fill_manual(values=c(COL[3],COL[4],COL[6]),labels=c("m1" = expression("Causal (V"%->%"T"%->%"G)"),"m2"=expression("Reactive (V"%->%"G"%->%"T)"),"m3"=expression("Independent (V"%->%"T & V"%->%"G)")))+    #scale_x_discrete(labels=c("normal"=paste("normal\n(N=",nrow(df),")",sep=""),"tumor"=paste("tumor\n(N=",nrow(results_tumor),")",sep="")))+
    labs(y="Model probability",x="",title=paste("Shared triplets (N=",nrow(df),")\ngenome-wide mean probability",sep=""))+
    theme_classic(base_size=20)+
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="top",
        legend.direction="vertical",
        legend.justification="left",
        legend.title=element_blank())
ggsave(gg_dd, filename="SHARED.model_percentages_genomeWide.pdf",device="pdf",height=7,width=7)


########## Testing new plot ################ 


df <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/shared/NORMAL_TUMOR_shared_BNlearn.txt.gz",header=T,stringsAsFactors=F,sep="\t")


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

supp_labs <- c("C"="Causal model\nchanges","R"="Reactive model\nchanges","I"="Independent model\nchanges")
gg_baby <- ggplot(data=d,aes(x=tissue, y=freq,fill=model,label=freq))+
     geom_bar(stat="identity",colour="white",size=1,width = 0.9)+
     geom_text(size=5,position=position_stack(vjust=0.5))+
     facet_wrap(~key,scales="free",labeller = labeller(key = supp_labs))+
    scale_fill_manual(values=c(COL[3],COL[4],COL[6]),labels=c("C" = expression("Causal (V"%->%"T"%->%"G)"),"R"=expression("Reactive (V"%->%"G"%->%"T)"),"I"=expression("Independent (V"%->%"T & V"%->%"G)")))+    #scale_x_discrete(labels=c("normal"=paste("normal\n(N=",nrow(df),")",sep=""),"tumor"=paste("tumor\n(N=",nrow(results_tumor),")",sep="")))+
    labs(x="",y="frequency",title="Model substitutions between\nnormal and tumor\nfor shared triplets (N=1571)")+
     theme_bw(base_size=20)+
     theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="top",
        legend.direction="vertical",
        legend.justification="left",
        legend.title=element_blank())

ggsave(gg_baby, filename="model_substitutions_shared_triplets.pdf",device="pdf",height=8,width=10)


##### MODEL SIGNIFICANT DIFFERENCES BETWEEN NORMAL AND TUMOR ###### 


############## FISHER EXACT TEST COMPARISON BETWEEN EACH MODEL #################

# Causal model (m1)

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

# model    p.value      oddsratio
#  m1     2.632708e-37    0.3333315
# m2      4.364760e-06    0.7113500
# m3      3.204723e-66    4.1381766