#!/bin/bash 

module load subread

RAW_DIR="/scratch/lykoskou/te_project/data_syscol"
GTF="/scratch/lykoskou/te_project/hg19_genes_TE_1709.gtf"

for bam in ${RAW_DIR}/*/*.bam; 
do 
    SAMPLE_ID=$(basename ${bam} | cut -d"." -f1)

    #### UNIQUELY MAPPED 
    wsbatch --partition=serial --mem=40G --time=03:00:00 -N1 -n1 -c4 --wrap="featureCounts -p -T 4 -t exon -g gene_id -a ${GTF} -Q 10 -o ${SAMPLE_ID}.ufc.txt ${bam}"

    #### MULTI_MAPPED 
    wsbatch --partition=serial --mem=40G --time=03:00:00 -N1 -n1 -c4 --wrap="featureCounts -M --fraction -p -T 4 -t exon -g gene_id -a ${GTF} -Q 0 -o ${SAMPLE_ID}.mfc.txt ${bam}"

done 