#!/usr/bin/env Rscript 

ALL <- as.data.frame(data.table::fread("/Users/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_TEannotation/hg19_TE_repmask_LTRm_s_20140131.maptable.bed.gz",header=TRUE,stringsAsFactors = FALSE,sep="\t"))
ALL$id <- paste(ALL$id, ALL$start+1, ALL$end, ALL$strand, sep="_")

tab <- as.data.frame(prop.table(table(ALL$gid)))
tab$Var1 <- as.character(tab$Var1)
tab$family <- sapply(strsplit(tab$Var1, "\\/"),"[[",1)
tab <- tab[!grepl("\\?",tab$family),]

gg1 <- ggplot(data=tab, aes(x=Var1, y=Freq,fill=family))+
  geom_bar(stat="identity",color="black")+
  labs(x="TE subfamilies",y="Proportion",title="proportion of TE subfamilies in the genome", fill="TE family")+
  scale_fill_npg()+
  theme_base(base_size=12)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),legend.title=element_text(face="bold"),
        plot.title=element_text(hjust=0.5), axis.title = element_text(face="bold"),
        plot.background=element_blank())
ggsave(gg1, filename="suppFigure1.pdf",height=6,width=8,device="pdf")

