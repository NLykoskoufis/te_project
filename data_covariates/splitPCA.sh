#!/bin/bash 
set -e 

INPUT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_covariates/SYSCOL_normal_3PC_geno_merged.pca"

for i in 0 1 2 5 10 20 30 40 50 60 70 80 90 100; do 
    OUT="NORMAL.PC${i}_PC3_geno.pca"
    head -$((4+i)) ${INPUT} > /srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_covariates/normal/${OUT}
done

INPUT2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_covariates/SYSCOL_tumor_3PC_geno_merged.pca"

for i in 0 1 2 5 10 20 30 40 50 60 70 80 90 100; do
    OUT="TUMOR.PC${i}_PC3_geno.pca"
    head -$((4+i)) ${INPUT2} > /srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_covariates/tumor/${OUT}
done




