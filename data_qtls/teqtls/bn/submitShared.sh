#!/bin/bash 
set -e

BIN="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/binaries/data_qtls/teqtls/bn/find_shared_triplets.py"
ROOT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/shared"

TUMOR="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/tumor/TUMOR.bnlearn_results.txt.gz"
NORMAL="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/normal/NORMAL.bnlearn_results.txt.gz"
TEQTLT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/permute/tumor/TUMOR_teqtls.chrALL.significant.txt"
TEQTLN="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/permute/normal/NORMAL_teqtls.chrALL.significant.txt"
VCF="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_genotypes/syscol_276_cancerNames_INFO_filtered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
R2=0.90
OUT="${ROOT}/NORMAL_TUMOR_shared_BNlearn.txt"


wsbatch -p shared-cpu -n 1 -t 12:00:00 --mem=20G \
    --wrap="python3 ${BIN} \
            -t ${TUMOR} \
            -n ${NORMAL} \
            -teqtl ${TEQTLT} \
            -neqtl ${TEQTLN} \
            -vcf ${VCF} \
            -r2 ${R2} \
            -out ${OUT}"
