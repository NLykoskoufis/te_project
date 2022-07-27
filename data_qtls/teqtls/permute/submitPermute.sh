#!/bin/bash 

BIN="/srv/beegfs/scratch/shares/brauns_lab/bin/QTLtools"

# NORMAL 
VCFN="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_genotypes/syscol_275_normalNames_INFOfiltered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
BEDN="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/emods/quantifications/normal/NORMAL_TE_GENE_pc1_chrALL.bed.gz"
OUT_DIR="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/permute/normal"

#TUMOR
VCFT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_genotypes/syscol_276_cancerNames_INFO_filtered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
BEDT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/emods/quantifications/tumor/TUMOR_TE_GENE_pc1_chrALL.bed.gz"
OUT_DIR2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/permute/tumor"


for c in {0..100}; do 
    # NORMAL
    wsbatch --partition=shared-cpu --time=12:00:00 --mem=10G -o ${OUT_DIR}/slurm-%j_${c}_100.out --wrap="${BIN} cis --vcf ${VCFN} --bed ${BEDN} --permute 1000 --normal --chunk ${c} 100 --out ${OUT_DIR}/NORMAL_teqtls_${c}_100.txt.gz";
    # TUMOR
    wsbatch --partition=shared-cpu --time=12:00:00 --mem=10G -o ${OUT_DIR2}/slurm-%j_${c}_100.out --wrap="${BIN} cis --vcf ${VCFT} --bed ${BEDT} --permute 1000 --normal --chunk ${c} 100 --out ${OUT_DIR2}/TUMOR_teqtls_${c}_100.txt.gz";
done 
