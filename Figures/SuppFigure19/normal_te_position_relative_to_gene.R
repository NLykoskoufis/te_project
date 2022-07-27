
suppressMessages(library(ggplot2))
suppressMessages(library(qvalue))
library(dplyr)
library(RColorBrewer)
COL = brewer.pal(9,"Set3")
library(reshape2)


# USEFUL FUNCTIONS # 

getBest <- function(x){
    posteriors <- as.numeric(x[4:6])
    max_model <- max(posteriors)
    if(posteriors[1] == max_model){return("m1")}
    if(posteriors[2] == max_model){return("m2")}
    if(posteriors[3] == max_model){return("m3")}
}



###################### POSITION OF TE COMPARED TO GENE ########################

TPM <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/normal/NORMAL.BNlearn.ALL.extraINFO.txt.gz", header = TRUE,stringsAsFactors = FALSE,sep="\t"))

N = rep(0,7)
C = rep(0,7)
R = rep(0,7)
I = rep(0,7)

N[1] = nrow(TPM[TPM$ting,])
N[2] = nrow(TPM[!TPM$ting,])
N[3] = nrow(TPM[TPM$updown == "up",])
N[4] = nrow(TPM[TPM$updown == "down",])
N[5] = nrow(TPM[TPM$vint,])
N[6] = nrow(TPM[!TPM$vint,])
N[7] = nrow(TPM)

C[1] = sum(TPM[TPM$ting,]$M1) / N[1]
C[2] = sum(TPM[!TPM$ting,]$M1) / N[2]
C[3] = sum(TPM[TPM$updown == "up",]$M1) / N[3]
C[4] = sum(TPM[TPM$updown == "down",]$M1) / N[4]
C[5] = sum(TPM[TPM$vint,]$M1) / N[5]
C[6] = sum(TPM[!TPM$vint,]$M1) / N[6]
C[7] = sum(TPM$M1) / N[7]

R[1] = sum(TPM[TPM$ting,]$M2) / N[1]
R[2] = sum(TPM[!TPM$ting,]$M2) / N[2]
R[3] = sum(TPM[TPM$updown == "up",]$M2) / N[3]
R[4] = sum(TPM[TPM$updown == "down",]$M2) / N[4]
R[5] = sum(TPM[TPM$vint,]$M2) / N[5]
R[6] = sum(TPM[!TPM$vint,]$M2) / N[6]
R[7] = sum(TPM$M2) / N[7]

I[1] = sum(TPM[TPM$ting,]$M3) / N[1]
I[2] = sum(TPM[!TPM$ting,]$M3) / N[2]
I[3] = sum(TPM[TPM$updown == "up",]$M3) / N[3]
I[4] = sum(TPM[TPM$updown == "down",]$M3) / N[4]
I[5] = sum(TPM[TPM$vint,]$M3) / N[5]
I[6] = sum(TPM[!TPM$vint,]$M3) / N[6]
I[7] = sum(TPM$M3) / N[7]
TOT=rbind(C, R, I)


TOT <- as.data.frame(TOT,stringsAsFactors=F)
colnames(TOT) <- c("TinG","ToutG","TEupGene","TEdownGene","VinT","VoutT","All")
rownames(TOT) <- c("C","R","I")
TOT$model <- rownames(TOT)

tot <- melt(TOT, id.vars="model")
gg_normal <- ggplot(data=tot, aes(x=variable, y=value*100,fill=model,label=paste0(round(value*100,0),"%")))+
    geom_bar(stat="identity",color="black",size=1.25)+
    geom_text(position=position_stack(vjust=0.5),size=5)+
    labs(y="Mean probability",x="",title=paste("Normal (N=",nrow(TPM),")",sep=""))+
    scale_x_discrete(labels=c("TinG"=paste("TE in Gene\n(N=",N[1],")",sep=""),
                              "TouG"=paste("TE out of Gene\n(N=",N[2],")",sep=""),
                              "TEupGene"=paste("TE upstream\nof gene\n(N=",N[3],")",sep=""),
                              "TEdownGene"=paste("TE downstream\nof gene\n(N=",N[4],")",sep=""),
                              "VinT"=paste("Variant inside\nTE\n(N=",N[5],")",sep=""),
                              "VoutT"=paste("Variant outside\nTE\n(N=",N[6],")",sep=""),
                              "All" = paste("All cases\nTogether\n(N=",N[7],")",sep="")))+
    theme_classic(base_size=20)+
    scale_fill_manual(values=c("C"=COL[3],"R"=COL[4],"I"=COL[6]),labels=c("I"="Independent (V->G & V->T)","R"="Reactive (V->G->T)","C" = "Causal (V->T->G)"))+
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          plot.title=element_text(hjust=0.5,face="bold"),
          legend.position="top",
          legend.direction="vertical",
          legend.justification="left",
          legend.title=element_blank(),
          axis.text.x=element_text())


ggsave(gg_normal, filename="normal_te_position_relative_to_gene_ggplot.pdf",device="pdf",height=7,width=12.5)



#### SIGNIFICANT TESTS ##### 

# TE in gene vs out of gene

N_in = nrow(TPM[TPM$ting,])
N_out = nrow(TPM[!TPM$ting,])
C_in = nrow(CTPM[CTPM$ting,])
C_out = nrow(CTPM[!CTPM$ting,])
R_in = nrow(RTPM[RTPM$ting,])
R_out = nrow(RTPM[!RTPM$ting,])
I_in = nrow(ITPM[ITPM$ting,])
I_out = nrow(ITPM[!ITPM$ting,])

c <- matrix(c(C_in,C_out,(N_in - C_in),(N_out-C_out)),nrow=2,ncol=2,byrow=TRUE)
r <- matrix(c(R_in,R_out,(N_in - R_in),(N_out-R_out)),nrow=2,ncol=2,byrow=TRUE)
i <- matrix(c(I_in,I_out,(N_in - I_in),(N_out-I_out)),nrow=2,ncol=2,byrow=TRUE)

c.f = fisher.test(c) # OR = 0.6380234 Pvalue  < 2.2e-16
r.f = fisher.test(r) # OR = 2.165599 Pvalue < 2.2e-16
i.f = fisher.test(i) # OR = 0.4983539 Pvalue = < 2.2e-16


# TE upstream vs downstream of gene

N_up = nrow(TPM[TPM$updown == "up",])
N_down = nrow(TPM[TPM$updown == "down",])
C_up = nrow(CTPM[CTPM$updown == "up",]) 
C_down = nrow(CTPM[CTPM$updown == "down",])
R_up = nrow(RTPM[RTPM$updown == "up",]) 
R_down = nrow(RTPM[RTPM$updown == "down",]) 
I_up = nrow(ITPM[ITPM$updown == "up",]) 
I_down = nrow(ITPM[ITPM$updown == "down",]) 


c1 <- matrix(c(C_up,C_down,(N_up - C_up),(N_down-C_down)),nrow=2,ncol=2,byrow=TRUE)
r1 <- matrix(c(R_up,R_down,(N_up - R_up),(N_down-R_down)),nrow=2,ncol=2,byrow=TRUE)
i1 <- matrix(c(I_up,I_down,(N_up - I_up),(N_down-I_down)),nrow=2,ncol=2,byrow=TRUE)

c1.f <- fisher.test(c1) # OR = 1.128016 Pval 0.03327
r1.f <- fisher.test(r1) # OR = 0.6987501 Pvalue 4.05e-13
i1.f <- fisher.test(i1) # OR = 1.413631 Pvalue  7.619e-10



