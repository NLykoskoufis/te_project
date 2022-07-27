#!/usr/bin/env Rscript 

library(ggplot2)
library(ggpubr)

if (! file.exists(c("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/specific/tissueSpecific_eQTLs.txt","/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/specific/shared_eQTLs.txt")))
{
    n <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/specific/NORMAL.specificity_All.txt.gz",header=T,stringsAsFactors=F,sep="\t"))
    n$key <- paste(n$phe_id, n$var_id, sep=";")
    
    cond.n <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/mapping/normal/NORMAL.conditional_bestPerRank.txt.gz",header=T,stringsAsFactors=F,sep=" "))
    cond.n$key <- paste(cond.n$phe_id, cond.n$var_id, sep=";")
    
    cond.t <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/mapping/tumor/TUMOR.conditional_bestPerRank.txt.gz",header=T,stringsAsFactors=F,sep=" "))
    cond.t$key <- paste(cond.t$phe_id, cond.t$var_id, sep=";")
    
    t <-  as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/specific/TUMOR.specificity_All.txt.gz",header=T,stringsAsFactors=F,sep="\t"))
    t$key <- paste(t$phe_id, t$var_id, sep=";")
    
    n_spe <- n[which(n$spe_hit == 1),]
    t_spe <- t[which(t$spe_hit == 1),]
    
    n_spe$tissue <- "normal"
    t_spe$tissue <- "tumor"
    
    tissueSpe <- rbind(n_spe, t_spe)
    
    cond.n$tissue <- "normal"
    cond.t$tissue <- "tumor"
    
    dd <- rbind(cond.n,cond.t)
    dd <- unique(dd[,c(7,23)])
    
    tissueSpe <- merge(tissueSpe,dd, by="key")
    tissueSpe <- tissueSpe[,-1]
    write.table(tissueSpe,file="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/specific/tissueSpecific_eQTLs.txt",col.names=T,row.names=F,sep=" ",quote=F)
    
    ### Shared eQTLs 
    
    n_shared <- n[which(n$eqtl_int_pval >= 0.05),]
    t_shared <- t[which(t$eqtl_int_pval >= 0.05),]
    
    pot_shared <- unique(c(n_shared$key, t_shared$key))
    count <- 1
    list_counter <- 1
    true_shared <- vector()
    shared <- list()
    for(r in pot_shared){
        cat(count, "-",length(pot_shared),"\r")
        count <- count + 1
        x.n <- cond.n[which(cond.n$key == r),]
        x.t <- cond.t[which(cond.t$key == r),]
        if(nrow(x.n)!=0  && nrow(x.t)!=0){
            if(x.n$bwd_pval < 0.05 && x.t$bwd_pval < 0.05){
                true_shared <- c(true_shared, r)
                phe_id = sapply(strsplit(r, "\\;"),"[[",1)
                var_id = sapply(strsplit(r, "\\;"),"[[",2)
                phe_chr = x.n$phe_chr 
                phe_start = x.n$phe_from 
                phe_end = x.n$phe_to 
                phe_strand = x.n$phe_strd
                var_chr = x.n$var_chr 
                var_from = x.n$var_from 
                var_to = x.n$var_to            
                n_eqtl_bwd_pval = x.n$bwd_pval 
                n_eqtl_bwd_slope = x.n$bwd_slope 
                t_eqtl_bwd_pval = x.t$bwd_pval 
                t_eqtl_bwd_slope = x.t$bwd_slope 
                dist_phe_var = x.n$dist_phe_var
    
    
                d <- data.frame("phe_id"=phe_id,
                                "phe_chr"=phe_chr, 
                                "phe_from"=phe_start,
                                "phe_to"=phe_end,
                                "phe_strd"=phe_strand,
                                "n_var_in_cis"=".",
                                "dist_phe_var"=dist_phe_var,
                                "var_id"=var_id,
                                "var_chr"=var_chr,
                                "var_from"=var_from,
                                "var_to"=var_to,
                                "n_eqtl_bwd_pval"=n_eqtl_bwd_pval,
                                "n_eqtl_bwd_slope"=n_eqtl_bwd_slope,
                                "t_eqtl_bwd_pval"=t_eqtl_bwd_pval,
                                "t_eqtl_bwd_slope"=t_eqtl_bwd_slope)
                shared[[list_counter]] <- d
                list_counter <- list_counter + 1
            }
        }
    
    }
    shared <- Reduce(rbind, shared)
    write.table(shared, file="/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/specific/shared_eQTLs.txt",quote=F,sep=" ",row.names=F,col.names=T)
}else{
    tissueSpe <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/specific/tissueSpecific_eQTLs.txt",header=T,stringsAsFactors=F,sep=" ")
    shared <- read.table("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/specific/shared_eQTLs.txt",header=T,stringsAsFactors=F,sep=" ")  
}
    


gg_te <- ggplot()+
    geom_point(data=tissueSpe[!grepl("ENSG",tissueSpe$phe_id),], aes(x=dist_phe_var.y,y=-log10(bwd_pval),colour="tissue-specific"))+
    geom_point(data=shared[!grepl("ENSG",shared$phe_id),], aes(x=dist_phe_var,y=-log10(t_eqtl_bwd_pval),colour="shared"))+
    #geom_point(data=tissueSpe[grepl("ENSG",tissueSpe$phe_id),], aes(x=dist_phe_var,y=log10(bwd_pval),colour="tissue-specific",shape="gene"))+
    #geom_point(data=shared[grepl("ENSG",shared$phe_id),], aes(x=dist_phe_var,y=log10(t_eqtl_bwd_pval),colour="shared",shape="gene"))+
    labs(x="distance (bp)",y="-log10(P-value)",title="TE-eQTLs")+
    theme_classic(base_size=20)+
    scale_colour_manual(values=c("tissue-specific"="red","shared"="black"))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.position=c(0.8,0.9),
          legend.title=element_blank(),
          legend.box.background = element_rect(color="black", size=0.5),
          legend.box.margin = margin(2, 4, 4, 4),
          plot.title=element_text(hjust=0.5,face="bold"))
ggsave(gg_te, filename="TE_tissue_specific_shared_eQTLs_dist.pdf",device="pdf",height=8,width=8)

gg_gene <- ggplot()+
    geom_point(data=tissueSpe[grepl("ENSG",tissueSpe$phe_id),], aes(x=dist_phe_var.y,y=-log10(bwd_pval),colour="tissue-specific"))+
    geom_point(data=shared[grepl("ENSG",shared$phe_id),], aes(x=dist_phe_var,y=-log10(t_eqtl_bwd_pval),colour="shared"))+
    labs(x="distance (bp)",y="-log10(P-value)",title="Gene-eQTLs")+
    theme_classic(base_size=20)+
    scale_colour_manual(values=c("tissue-specific"="red","shared"="black"))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.position=c(0.8,0.9),
          legend.title=element_blank(),
          legend.box.background = element_rect(color="black", size=0.5),
          legend.box.margin = margin(2, 4, 4, 4),
          plot.title=element_text(hjust=0.5,face="bold"))
ggsave(gg_gene, filename="gene_tissue_specific_shared_eQTLs_dist.pdf",device="pdf",height=8,width=8)


gg_all <-  ggplot()+
    geom_point(data=tissueSpe, aes(x=dist_phe_var.y,y=-log10(bwd_pval),colour="tissue-specific"))+
    geom_point(data=shared, aes(x=dist_phe_var,y=-log10(t_eqtl_bwd_pval),colour="shared"))+
    labs(x="distance (bp)",y="-log10(P-value)",title="TE and Gene-eQTLs")+
    theme_classic(base_size=20)+
    scale_colour_manual(values=c("tissue-specific"="red","shared"="black"))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.position=c(0.8,0.9),
          legend.title=element_blank(),
          legend.box.background = element_rect(color="black", size=0.5),
          legend.box.margin = margin(2, 4, 4, 4),
          plot.title=element_text(hjust=0.5,face="bold"))
ggsave(gg_all, filename="all_tissue_specific_shared_eQTLs.pdf",device="pdf",width=8,height=8)


### Diff in distance between shared and tissueSpe?
wilcox.test(abs(shared$dist_phe_var),abs(tissueSpe$dist_phe_var.x)) 
#Wilcox rank sum test Pvalue < 2.2e-16
# median abs shared distance 6'557 bp
# median abs tissueSpe distance 51'868 bp
 
### Diff in distance between shared and tissueSPe for TE-eQTLs
wilcox.test(abs(shared[!grepl("ENSG",shared$phe_id),]$dist_phe_var),abs(tissueSpe[!grepl("ENSG",tissueSpe$phe_id),]$dist_phe_var.y)) 
# Wilcox rank sum test P<2.2e-16
# median abs shared distance TE: 5'580.5
# median abs tissueSPe distance TE: 60'870

## mosaic plot 

data1 <- as.table(matrix(c(nrow(t_spe),nrow(n_spe),nrow(shared),nrow(shared)),byrow=T,2,2))
colnames(data1) <- c("Tumor","Normal")
rownames(data1) <- c("Unique to cell type","shared")
pdf("eqtl_level_sharing_mosaic_plot.pdf",8,8)
mosaicplot(t(data1),cex.axis=1.5,main="eQTL level sharing")
dev.off()
data1
#                    Tumor Normal
#Unique to cell type   546   2599
#shared                700    700

### mosaic plot for TE-eQTLs 

t_te <- tissueSpe[which(tissueSpe$tissue == "tumor"),]
t_te <- t_te[!grepl("ENSG",t_te$phe_id),]

n_te <- tissueSpe[which(tissueSpe$tissue == "normal"),]
n_te <- n_te[!grepl("ENSG",n_te$phe_id),]


data1 <- matrix(c(nrow(t_te),nrow(n_te),nrow(shared[!grepl("ENSG",shared$phe_id),]),nrow(shared[!grepl("ENSG",shared$phe_id),])),byrow=T,2,2)
colnames(data1) <- c("Tumor","Normal")
rownames(data1) <- c("Unique to cell type","shared")
pdf("TE_eqtl_level_sharing_mosaic_plot.pdf",8,8)
mosaicplot(t(data1),cex.axis=1.5,main="TE-eQTL level sharing")
dev.off()

#Tumor Normal
#Unique to cell type   429   1697
#shared                525    525

t_gene <- tissueSpe[which(tissueSpe$tissue == "tumor"),]
t_gene <- t_gene[grepl("ENSG",t_gene$phe_id),]
n_gene <- tissueSpe[which(tissueSpe$tissue == "normal"),]
n_gene <- n_gene[grepl("ENSG",n_gene$phe_id),]

data2 <- matrix(c(nrow(t_gene),nrow(n_gene),nrow(shared[grepl("ENSG",shared$phe_id),]),nrow(shared[grepl("ENSG",shared$phe_id),])),byrow=T,2,2)
colnames(data2) <- c("Tumor","Normal")
rownames(data2) <- c("Unique to cell type","shared")
pdf("Gene_eqtl_level_sharing_mosaic_plot.pdf",8,8)
mosaicplot(t(data2),cex.axis=1.5,main="Gene-eQTL level sharing")
dev.off()
data2
#                    Tumor Normal
#Unique to cell type   117    902
#shared                175    175



#### effect size plot #### 
gg_effectSize <- ggplot()+
    geom_point(data=shared, aes(x=n_eqtl_bwd_slope,y=t_eqtl_bwd_slope,colour="Shared"))+
    geom_point(data=tissueSpe[which(tissueSpe$tissue == "tumor"),], aes(x=nominal_slope, y=bwd_slope,colour="Tumor specific"))+
    geom_point(data=tissueSpe[which(tissueSpe$tissue == "normal"),], aes(y=nominal_slope, x=bwd_slope,colour="Normal specific"))+
    labs(x="Normal eQTL slope",y="Tumor eQTL slope")+
    theme_classic()+
    scale_colour_manual(values=c("Shared"="black","Tumor specific"="#e41a1c","Normal specific"="#377eb8"))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.position="top",
    legend.direction="vertical",
    legend.justification="left",
    legend.title=element_blank())

ggsave(gg_effectSize,filename="tissueSpe_effectSize_plot.pdf",device="pdf",height=5,width=5.5)

gg_effectSize_te <- ggplot()+
    geom_point(data=shared[!grepl("ENSG",shared$phe_id),], aes(x=n_eqtl_bwd_slope,y=t_eqtl_bwd_slope,colour="Shared"))+
    geom_point(data=t_te, aes(x=nominal_slope, y=bwd_slope, colour="Tumor specific"))+
    geom_point(data=n_te, aes(x=bwd_slope, y=nominal_slope, colour="Normal specific"))+
    labs(x="Normal TE-eQTL slope",y="Tumor TE-eQTL slope")+
    theme_classic()+
    scale_colour_manual(values=c("Shared"="black","Tumor specific"="#e41a1c","Normal specific"="#377eb8"))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.position="top",
    legend.direction="vertical",
    legend.justification="left",
    legend.title=element_blank())
ggsave(gg_effectSize_te, filename="TE_tissu e_specific_shared_slopes.pdf",device="pdf",height=5, width=5.5)
##### TE family enrichment ########
# tissueSPe 

gg_effectSize_gene <- ggplot()+
    geom_point(data=shared[grepl("ENSG",shared$phe_id),], aes(x=n_eqtl_bwd_slope,y=t_eqtl_bwd_slope,colour="Shared"))+
    geom_point(data=t_gene, aes(x=nominal_slope, y=bwd_slope, colour="Tumor specific"))+
    geom_point(data=n_gene, aes(x=bwd_slope, y=nominal_slope, colour="Normal specific"))+
    labs(x="Normal Gene-eQTL slope",y="Tumor Gene-eQTL slope")+
    theme_classic()+
    scale_colour_manual(values=c("Shared"="black","Tumor specific"="#e41a1c","Normal specific"="#377eb8"))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    legend.position="top",
    legend.direction="vertical",
    legend.justification="left",
    legend.title=element_blank())
ggsave(gg_effectSize_gene, filename="Gene_tissue_specific_shared_slopes.pdf",device="pdf",height=5, width=5.5)


