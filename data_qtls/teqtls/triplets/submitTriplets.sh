#!/bin/bash 

BIN="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/binaries/data_qtls/teqtls/triplets/create_triplets.py"


# NORMAL 
TEQTLN="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/permute/normal/NORMAL_teqtls.chrALL.significant.txt"
VCFN="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_genotypes/syscol_275_normalNames_INFOfiltered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
BEDN="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/normal/SYSCOL_normal_275_gene_TE.FM05.resid.cpm.bed.gz"
OUT_DIR="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/triplet/normal"

cd ${OUT_DIR}
wsbatch --partition=shared-cpu --time=12:00:00 --mem=10G --wrap="python3 ${BIN} ${TEQTLN} ${VCFN} ${BEDN}"


#TUMOR
TEQTLT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/permute/tumor/TUMOR_teqtls.chrALL.significant.txt"
VCFT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_genotypes/syscol_276_cancerNames_INFO_filtered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
BEDT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/tumor/SYSCOL_tumor_276_gene_TE.FM05.resid.cpm.bed.gz"
OUT_DIR2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/triplet/tumor"

cd ${OUT_DIR2}
wsbatch --partition=shared-cpu --time=12:00:00 --mem=10G --wrap="python3 ${BIN} ${TEQTLT} ${VCFT} ${BEDT}"

cd /srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/triplet
