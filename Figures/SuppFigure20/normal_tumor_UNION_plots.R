#!/usr/bin/env Rscript 


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
cpus <- 20
registerDoMC(cpus)


results <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/union/UNION_normal_tumor_with_posteriors.txt.gz",sep="\t",header=T,stringsAsFactors=F)

NO <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/normal/NORMAL.BNlearn.ALL.extraINFO.txt.gz",sep="\t",header=T,stringsAsFactors = FALSE)

TU <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/tumor/TUMOR.BNlearn.ALL.extraINFO.txt.gz",sep="\t",header=T,stringsAsFactors = FALSE)

df <- results

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
    scale_fill_manual(values=c(COL[3],COL[4],COL[6]),labels=c("C" = expression("Causal (V"%->%"T"%->%"G)"),"R"=expression("Reactive (V"%->%"G"%->%"T)"),"I"=expression("Independent (V"%->%"T & V"%->%"G)")))+    #scale_x_discrete(labels=c("normal"=paste("normal\n(N=",nrow(df),")",sep=""),"tumor"=paste("tumor\n(N=",nrow(results_tumor),")",sep="")))+
    labs(x="",y="frequency",title=paste("Model substitutions for the union between\nnormal and tumor triplets\n(N=",nrow(df),")",sep=""))+
     theme_bw(base_size=20)+
     theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="top",
        legend.direction="vertical",
        legend.justification="left",
        legend.title=element_blank())

ggsave(ggbaby2, filename="model_substitutions_normal_tumor_UNION.pdf",device="pdf",height=8,width=10)



n_model_perc <- as.data.frame(table(results$tumor_best))
n_model_perc$perc <- n_model_perc$Freq / sum(n_model_perc$Freq)
n_model_perc$type <- "tumor"
n_model_perc$Var1 <- as.character(n_model_perc$Var1)
n_model_perc <- n_model_perc[order(n_model_perc$Var1),]

n_model_tperc <- as.data.frame(table(results$normal_best))
n_model_tperc$perc <- n_model_tperc$Freq / sum(n_model_tperc$Freq)
n_model_tperc$type <- "normal"
n_model_tperc$Var1 <- as.character(n_model_tperc$Var1)
n_model_tperc <- n_model_tperc[order(n_model_tperc$Var1),]

dn <- rbind(n_model_perc,n_model_tperc)


############## FISHER EXACT TEST COMPARISON BETWEEN EACH MODEL #################

# Causal model (m1)

fisher <- list()
for(i in 1:3){
    m <- matrix(c(n_model_perc[i,]$Freq,n_model_tperc[i,]$Freq,
                (sum(n_model_perc$Freq) - n_model_perc[i,]$Freq),(sum(n_model_tperc$Freq) - n_model_tperc[i,]$Freq)),
                ncol=2,nrow=2,byrow=T)

    fish <- fisher.test(m)
    fish_df <- data.frame("model" = paste("m",i,sep=""),
                        "p.value" = fish$p.value,
                        "oddsratio" = fish$estimate)
    fisher[[i]] <- fish_df
}
r <- Reduce(rbind,fisher)

# model    p.value      oddsratio
# m1    8.810239e-167   2.2200365
# m2    1.181087e-37    0.7113492
# m3    5.325404e-64    0.5452098


#################################################################################



COL = brewer.pal(9,"Set3")
gg2 <- ggplot(dn,aes(x=factor(type),y=perc*100,fill=factor(Var1),label=paste0(round(perc*100,0),"%\n(",Freq,")")))+
    geom_bar(stat="identity",size=1.25,color="black")+
    geom_text(position=position_stack(vjust=0.5),size=6)+
    scale_fill_manual(values=c("m1"=COL[3],"m2"=COL[4],"m3"=COL[6]),labels=c("m1" = "Causal (V->T->G)","m2"="Reactive (V->G->T)","m3"="Independent (V->G & V->T)"))+
    labs(y="Model percentage",x="",title=paste("Union of normal and tumor triplets\n(N=",nrow(results),")",sep=""))+
    theme_classic(base_size=20)+
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="top",
        legend.direction="vertical",
        legend.justification="left",
        legend.title=element_blank())
ggsave(gg2, filename="normal_tumor_UNION_model_percentages.pdf",device="pdf",height=8,width=8)

#gg_nt <- ggarrange(gg2,ggbaby2, ncol=2,nrow=1,common.legend=TRUE)
ggsave(gg2, filename="normal_tumor_UNION_model_percentages.pdf",device='pdf',height=7,width=7)
