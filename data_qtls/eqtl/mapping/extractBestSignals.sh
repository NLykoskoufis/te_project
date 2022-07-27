#!/bin/bash
set -e

NORMAL="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/mapping/normal"
TUMOR="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/mapping/tumor"

# NORMAL
zcat ${NORMAL}/*.txt.gz > ${NORMAL}/NORMAL.conditional_All.txt
awk 'NR==1; FNR>1 {if(($21 == 1 && $22 == 1)) {print}}' ${NORMAL}/NORMAL.conditional_All.txt > ${NORMAL}/NORMAL.conditional_bestPerRank.txt 
bgzip ${NORMAL}/NORMAL.conditional_All.txt
bgzip ${NORMAL}/NORMAL.conditional_bestPerRank.txt


# TUMOR 

zcat ${TUMOR}/*.txt.gz > ${TUMOR}/TUMOR.conditional_All.txt 
awk 'NR==1; FNR>1 {if(($21 == 1 && $22 == 1)) {print}}' ${TUMOR}/TUMOR.conditional_All.txt > ${TUMOR}/TUMOR.conditional_bestPerRank.txt
bgzip ${TUMOR}/TUMOR.conditional_All.txt 
bgzip ${TUMOR}/TUMOR.conditional_bestPerRank.txt
