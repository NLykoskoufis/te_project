#!/usr/bin/env Rscript 

library(RColorBrewer)
COL <- brewer.pal(9,"Set3")
library(ggplot2)
library(ggthemes)
library(ggpubr)


setwd("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/analysis_qtlReplication/teqtl")

trans <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/analysis_qtlReplication/teqtl/bn/normal/SYSCOL.NORMAL.in.GTEx.ColonTransverse.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t")) 
syscol <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/normal/NORMAL.BNlearn.ALL.extraINFO.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))

trans_N <- nrow(trans)
trans_C <- sum(trans$M1) / trans_N
trans_R <- sum(trans$M2) / trans_N
trans_I <- sum(trans$M3) / trans_N
trans_T <- rbind(as.vector(trans_C),as.vector(trans_R),as.vector(trans_I))
trans_df <- data.frame("models" = c("Causal","Reactive","Independent"), "model_prob"=trans_T)

syscol_N <- nrow(syscol)
syscol_C <- sum(syscol$M1) / syscol_N
syscol_R <- sum(syscol$M2) / syscol_N
syscol_I <- sum(syscol$M3) / syscol_N
syscol_T <- rbind(as.vector(syscol_C),as.vector(syscol_R),as.vector(syscol_I))
syscol_df <- data.frame("models" = c("Causal","Reactive","Independent"), "model_prob"=syscol_T)



trans_df$tissue <- "GTEx Colon Transverse"
syscol_df$tissue <- "SYSCOL normal"

toPlot <- rbind(trans_df,syscol_df)
ggReplication <- ggplot(data=toPlot, aes(x=models, y=model_prob,fill=models,label=paste0(round(model_prob*100,0),"%")))+
  geom_bar(stat="identity",color="black",size=1.15)+
  geom_text(position=position_stack(vjust=0.5),size=5)+
  theme_base(base_size=15)+
  facet_wrap(~tissue,)+
  scale_fill_manual(values=c("Causal"=COL[3],"Reactive"=COL[4],"Independent"=COL[6]))+
  labs(y="Model probability",x="")+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="top",
        legend.title=element_blank(),
        strip.text.x = element_text(face = "bold"))+
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave(ggReplication, filename="bnReplications_normal.pdf",height=5,width=15)
  



#### COLON TRANSVERSE #### 
cts <- merge(trans, syscol,by.y="triplet", by.x="syscol_triplet")
cts$colonTransverse_syscol_models <- paste(cts$best.y, cts$best.x,sep="->")
cts$colonTransverse_syscol_models <- gsub("m1","C",cts$colonTransverse_syscol_models)
cts$colonTransverse_syscol_models <- gsub("m2","R",cts$colonTransverse_syscol_models)
cts$colonTransverse_syscol_models <- gsub("m3","I",cts$colonTransverse_syscol_models)
cts$dfs <- "Same model\npredicted"
cts$dfs[which(cts$best.x != cts$best.y)] <- "Different model\npredicted" 



m <- as.data.frame(table(cts$colonTransverse_syscol_models))

plotXlabels <- c("I->R" = expression("I" * symbol('\256') * "R"),
                 "I->C" = expression("I" * symbol('\256') * "C"),
                 "I->I" = expression("I" * symbol('\256') * "I"),
                 "R->I" = expression("R" * symbol('\256') * "I"),
                 "C->I" = expression("C" * symbol('\256') * "I"),
                 "C->R" = expression("C" * symbol('\256') * "R"),
                 "R->C" = expression("R" * symbol('\256') * "C"),
                 "R->R" = expression("R" * symbol('\256') * "R"),
                 "C->C" = expression("C" * symbol('\256') * "C")) 

ggCS <- ggplot(data=m, aes(x=reorder(Var1,-Freq), y = Freq))+
  geom_bar(stat="identity",colour="black", fill="lightgrey")+
  labs(x=expression("SYSCOL " * symbol('\256') * " GTEx Colon Transvserse models"),y="Frequency")+
  scale_x_discrete(labels = plotXlabels)+
  coord_flip()+
  theme_base()+
  theme(plot.title=element_text(hjust=0.5, face="bold"))

d <- as.data.frame(prop.table(table(cts$dfs)))
ggColTran <- ggplot(d, aes(x=Var1, y=Freq*100, label=paste0(round(Freq * 100), "%")))+
  geom_bar(stat="identity", color="black", fill="lightgrey")+
  geom_text(position=position_stack(vjust=0.5), size=5, fontface="bold")+
  labs(x="", y="Percentage of same models\nbetween SYSCOL and GTEx Colon Transverse")+
  theme_base()
gg <- ggpubr::ggarrange(ggReplication,ggarrange(ggColTran,ggCS, ncol=2, nrow=1),nrow=2)  

ggsave(gg, filename="ggColTrans_syscol_GTEx_ColonTransverse_modelDifferences.pdf",height=11,width=11)


##### TUMOR ######


trans <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/analysis_qtlReplication/teqtl/bn/tumor/SYSCOL.TUMOR.in.TCGA_COAD.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t")) 
syscol <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/tumor/TUMOR.BNlearn.ALL.extraINFO.txt.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))

trans_N <- nrow(trans)
trans_C <- sum(trans$M1) / trans_N
trans_R <- sum(trans$M2) / trans_N
trans_I <- sum(trans$M3) / trans_N
trans_T <- rbind(as.vector(trans_C),as.vector(trans_R),as.vector(trans_I))
trans_df <- data.frame("models" = c("Causal","Reactive","Independent"), "model_prob"=trans_T)

syscol_N <- nrow(syscol)
syscol_C <- sum(syscol$M1) / syscol_N
syscol_R <- sum(syscol$M2) / syscol_N
syscol_I <- sum(syscol$M3) / syscol_N
syscol_T <- rbind(as.vector(syscol_C),as.vector(syscol_R),as.vector(syscol_I))
syscol_df <- data.frame("models" = c("Causal","Reactive","Independent"), "model_prob"=syscol_T)



trans_df$tissue <- "TCGA-COAD"
syscol_df$tissue <- "SYSCOL tumor"

toPlot <- rbind(trans_df,syscol_df)
ggReplication <- ggplot(data=toPlot, aes(x=models, y=model_prob,fill=models,label=paste0(round(model_prob*100,0),"%")))+
  geom_bar(stat="identity",color="black",size=1.15)+
  geom_text(position=position_stack(vjust=0.5),size=5)+
  theme_base(base_size=15)+
  facet_wrap(~tissue,)+
  scale_fill_manual(values=c("Causal"=COL[3],"Reactive"=COL[4],"Independent"=COL[6]))+
  labs(y="Model probability",x="")+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="top",
        legend.title=element_blank(),
        strip.text.x = element_text(face = "bold"))+
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave(ggReplication, filename="bnReplications_tumor.pdf",height=5,width=15)




#### COLON TRANSVERSE #### 
cts <- merge(trans, syscol,by.y="triplet", by.x="syscol_triplet")
cts$colonTransverse_syscol_models <- paste(cts$best.y, cts$best.x,sep="->")
cts$colonTransverse_syscol_models <- gsub("m1","C",cts$colonTransverse_syscol_models)
cts$colonTransverse_syscol_models <- gsub("m2","R",cts$colonTransverse_syscol_models)
cts$colonTransverse_syscol_models <- gsub("m3","I",cts$colonTransverse_syscol_models)
cts$dfs <- "Same model\npredicted"
cts$dfs[which(cts$best.x != cts$best.y)] <- "Different model\npredicted" 



m <- as.data.frame(table(cts$colonTransverse_syscol_models))

plotXlabels <- c("I->R" = expression("I" * symbol('\256') * "R"),
                 "I->C" = expression("I" * symbol('\256') * "C"),
                 "I->I" = expression("I" * symbol('\256') * "I"),
                 "R->I" = expression("R" * symbol('\256') * "I"),
                 "C->I" = expression("C" * symbol('\256') * "I"),
                 "C->R" = expression("C" * symbol('\256') * "R"),
                 "R->C" = expression("R" * symbol('\256') * "C"),
                 "R->R" = expression("R" * symbol('\256') * "R"),
                 "C->C" = expression("C" * symbol('\256') * "C")) 

ggCS <- ggplot(data=m, aes(x=reorder(Var1,-Freq), y = Freq))+
  geom_bar(stat="identity",colour="black", fill="lightgrey")+
  labs(x=expression("SYSCOL " * symbol('\256') * " TCGA-COAD models"),y="Frequency")+
  scale_x_discrete(labels = plotXlabels)+
  coord_flip()+
  theme_base()+
  theme(plot.title=element_text(hjust=0.5, face="bold"))

d <- as.data.frame(prop.table(table(cts$dfs)))
ggColTran <- ggplot(d, aes(x=Var1, y=Freq*100, label=paste0(round(Freq * 100), "%")))+
  geom_bar(stat="identity", color="black", fill="lightgrey")+
  geom_text(position=position_stack(vjust=0.5), size=5, fontface="bold")+
  labs(x="", y="Percentage of same models\nbetween SYSCOL and TCGA-COAD")+
  theme_base()
gg <- ggpubr::ggarrange(ggReplication,ggarrange(ggColTran,ggCS, ncol=2, nrow=1),nrow=2)  

ggsave(gg, filename="ggColTrans_syscol_TCGA_COAD_modelDifferences.pdf",height=11,width=11)
