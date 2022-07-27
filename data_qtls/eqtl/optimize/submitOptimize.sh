#!/bin/bash 
set -e 

BIN="/srv/beegfs/scratch/groups/funpopgen/bin/QTLtools"
VCF="/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_genotypes/syscol_275_normalNames_INFOfiltered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
BED="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/normal/SYSCOL_normal_275_gene_TE.FM05.raw.cpm.bed.gz"
PC_DIR="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_covariates/normal"
OUT_DIR="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/optimize/normal"

# NORMAL 
for PC in 0 1 2 5 10 20 30 40 50 60 70 80 90 100; do 
    echo ${PC}
    for i in {0..50}; do 
        wsbatch --partition=shared-cpu -o ${OUT_DIR}/slurm-%j.out --mem=10G --time=12:00:00 --wrap="${BIN} cis --vcf ${VCF} --bed ${BED} --cov ${PC_DIR}/NORMAL.PC${PC}_PC3_geno.pca --normal --permute 1000 --out ${OUT_DIR}/NORMAL.PC${PC}.${i}_50.txt.gz --chunk ${i} 50"; 
    done
done  

VCF2="/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_genotypes/syscol_276_cancerNames_INFO_filtered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
BED2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/tumor/SYSCOL_tumor_276_gene_TE.FM05.raw.cpm.bed.gz"
PC_DIR2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_covariates/tumor"
OUT_DIR2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/optimize/tumor"


# TUMOR 
for PC in 0 1 2 5 10 20 30 40 50 60 70 80 90 100; do 
    echo ${PC}
    for i in {0..50}; do 
        wsbatch --partition=shared-cpu -o ${OUT_DIR2}/slurm-%j.out --mem=10G --time=12:00:00 --wrap="${BIN} cis --vcf ${VCF2} --bed ${BED2} --cov ${PC_DIR2}/TUMOR.PC${PC}_PC3_geno.pca --normal --permute 1000 --out ${OUT_DIR2}/TUMOR.PC${PC}.${i}_50.txt.gz --chunk ${i} 50"; 
    done
done  
