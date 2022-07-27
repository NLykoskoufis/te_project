#!/usr/bin/env Rscript 


library(ggplot2)
library(ggpubr)
library(reshape2)
source("/home/users/l/lykoskou/bin/tools.R")
library(gridExtra)


#### eQTL plotting

bed <- as.data.frame(data.table::fread("zcat /srv/beegfs/scratch/users/l/lykoskou/TE/V3/EXPRESSION/normal/corrected/CPM_trimmed_Genes_proteinCoding_lincRNA_TEs_filtered_50_percent_and_more_275_normal.resid.bed.gz",header=T,stringsAsFactors=F,sep="\t"))
vcf.n <- "/srv/beegfs/scratch/users/l/lykoskou/TE/V3/GENOTYPES/syscol_275_normalNames_INFOfiltered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"



col = list("TE"="#ef8a62", "GENE"="#67a9cf")


phen <- as.data.frame(data.table::fread("/srv/beegfs/scratch/users/l/lykoskou/TE/V2/EXPRESSION/phen_types.txt.gz",header=F,stringsAsFactors = F,sep="\t"))
colnames(phen) <- c("phe_id","type","fam")
phen$FAM <- "GENES"
phen$FAM[grepl("DNA",phen$fam)] <- "DNA"
phen$FAM[grepl("LINE",phen$fam)] <- "LINE"
phen$FAM[grepl("LTR",phen$fam)] <- "LTR"
phen$FAM[grepl("SINE",phen$fam)] <- "SINE"


n.c <- as.data.frame(data.table::fread("zcat /srv/beegfs/scratch/users/l/lykoskou/TE/V3/QTL/eqtls/mapping/normal/NORMAL.conditional_bestPerRank.txt.gz",header=T,stringsAsFactors = F,sep=" "))

plot_list <- list()
for(r in 1:nrow(n.c)){
      cat(r,"-",nrow(n.c))
      x <- n.c[r,]
      phe_id <- x$phe_id 
      var_id <- x$var_id 
      var_pos <- paste(x$var_chr,":",x$var_from,"-",x$var_to,sep="")

      gen.n <- as.data.frame(data.table::fread(cmd=paste("bcftools view ", vcf.n ," -r ", var_pos, " | grep -v '##'", sep=""),header=T,stringsAsFactors=F,sep="\t",colClasses=c("character")))
      gen.n <- gen.n[which(gen.n$ID == var_id) ,]
      ref <- gen.n$REF 
      alt <- gen.n$ALT
      
      gen.n <- gen.n[,-c(1:9)]
      gen.n <- as.data.frame(t(gen.n),stringsAsFactors=F)
      gen.n$samples <- rownames(gen.n)
      colnames(gen.n) <- c("genotype","samples")
      gen.n$GT <- sapply(strsplit(gen.n$genotype,"\\:"),"[[",1)
      gen.n$DS <- sapply(strsplit(gen.n$genotype,"\\:"),"[[",2)
      gen.n$GT <- sub("1\\|0","0\\|1",gen.n$GT)

      exp <- bed[which(bed$id == phe_id),]
      gene_name <- gsub(".*N=","",exp$gid)
      exp <- exp[,-c(1:6)]
      exp <- as.data.frame(t(exp),stringsAsFactors=F)
      exp$samples <- rownames(exp)
      colnames(exp) <- c("exp","samples")

      df <- merge(gen.n,exp,by="samples")

      gg1 <- ggplot(data=df,aes(x=as.factor(GT),y=exp))+
            geom_boxplot(size=1.25,color="black",fill="lightgrey")+
            geom_jitter(width=0.2)+
            labs(x="genotypes",y="normalized CPM")+
            scale_x_discrete(labels=c("0|0"=paste(ref,"|",ref,sep=""),"0|1"=paste(ref,"|",alt,sep=""),"1|1"=paste(alt,"|",alt,sep="")))+
            theme_bw(base_size=24)+
            theme(panel.grid.minor=element_blank(),
            panel.grid.major=element_blank(),
            plot.title=element_text(hjust=0.5))

      gg2 <- ggplot(data=df,aes(x=as.numeric(DS),y=exp))+
            geom_jitter(width=0.2)+
            geom_smooth(method="lm",se=F)+
            labs(x="Dosage",y="normalized CPM")+
            theme_bw(base_size=24)+
            theme(panel.grid.minor=element_blank(),
            panel.grid.major=element_blank(),
            plot.title=element_text(hjust=0.5))

      gg <- ggarrange(gg1,gg2,ncol=2,nrow=1)
      g <- annotate_figure(gg,top=text_grob(paste("var_id: ",var_id," ; phe_id: ",phe_id,"\nphe_id2: ",gene_name,"\np-value= ",x$bwd_pval," ;slope= ",x$bwd_slope,sep=""),size=24))
      plot_list[[r]] <- g
}
ml <- marrangeGrob(plot_list, nrow=1, ncol=1)
ggsave(ml,filename="NORMAL_conditional_eQTL_All.pdf",device="pdf",height=8,width=16)



n.t <- as.data.frame(data.table::fread("zcat /srv/beegfs/scratch/users/l/lykoskou/TE/V3/QTL/eqtls/mapping/tumor/TUMOR.conditional_bestPerRank.txt.gz",header=T,stringsAsFactors = F,sep=" "))
vcf.t <- "/srv/beegfs/scratch/users/l/lykoskou/TE/V2/GENOTYPES/syscol_276_cancerNames_INFO_filtered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"


plot_list <- list()
for(r in 1:nrow(n.t)){
      cat(r,"-",nrow(n.t))
      x <- n.t[r,]
      phe_id <- x$phe_id 
      var_id <- x$var_id 
      var_pos <- paste(x$var_chr,":",x$var_from,"-",x$var_to,sep="")

      gen.n <- as.data.frame(data.table::fread(cmd=paste("bcftools view ", vcf.t ," -r ", var_pos, " | grep -v '##'", sep=""),header=T,stringsAsFactors=F,sep="\t",colClasses=c("character")))
      gen.n <- gen.n[which(gen.n$ID == var_id) ,]
      ref <- gen.n$REF 
      alt <- gen.n$ALT
      
      gen.n <- gen.n[,-c(1:9)]
      gen.n <- as.data.frame(t(gen.n),stringsAsFactors=F)
      gen.n$samples <- rownames(gen.n)
      colnames(gen.n) <- c("genotype","samples")
      gen.n$GT <- sapply(strsplit(gen.n$genotype,"\\:"),"[[",1)
      gen.n$DS <- sapply(strsplit(gen.n$genotype,"\\:"),"[[",2)
      gen.n$GT <- sub("1\\|0","0\\|1",gen.n$GT)

      exp <- bed[which(bed$id == phe_id),]
      gene_name <- gsub(".*N=","",exp$gid)
      exp <- exp[,-c(1:6)]
      exp <- as.data.frame(t(exp),stringsAsFactors=F)
      exp$samples <- rownames(exp)
      colnames(exp) <- c("exp","samples")

      df <- merge(gen.n,exp,by="samples")

      gg1 <- ggplot(data=df,aes(x=as.factor(GT),y=exp))+
            geom_boxplot(size=1.25,color="black",fill="lightgrey")+
            geom_jitter(width=0.2)+
            labs(x="genotypes",y="normalized CPM")+
            scale_x_discrete(labels=c("0|0"=paste(ref,"|",ref,sep=""),"0|1"=paste(ref,"|",alt,sep=""),"1|1"=paste(alt,"|",alt,sep="")))+
            theme_bw(base_size=24)+
            theme(panel.grid.minor=element_blank(),
            panel.grid.major=element_blank(),
            plot.title=element_text(hjust=0.5))

      gg2 <- ggplot(data=df,aes(x=as.numeric(DS),y=exp))+
            geom_jitter(width=0.2)+
            geom_smooth(method="lm",se=F)+
            labs(x="Dosage",y="normalized CPM")+
            theme_bw(base_size=24)+
            theme(panel.grid.minor=element_blank(),
            panel.grid.major=element_blank(),
            plot.title=element_text(hjust=0.5))

      gg <- ggarrange(gg1,gg2,ncol=2,nrow=1)
      g <- annotate_figure(gg,top=text_grob(paste("var_id: ",var_id," ; phe_id: ",phe_id,"\nphe_id2: ",gene_name,"\np-value= ",x$bwd_pval," ;slope= ",x$bwd_slope,sep=""),size=24))
      plot_list[[r]] <- g
}
ml <- marrangeGrob(plot_list, nrow=1, ncol=1)
ggsave(ml,filename="TUMOR_conditional_eQTLs_All.pdf",device="pdf",height=8,width=16)
