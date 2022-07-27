#!/bin/bash 
set -e

BIN="/srv/beegfs/scratch/groups/funpopgen/bin/QTLtools"

VCF="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_genotypes/syscol_275_normalNames_INFOfiltered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
BED="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/normal/SYSCOL_normal_275_gene_TE.FM05.resid.cpm.bed.gz"
OUT_DIR="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/mapping/normal"
THRESHOLD="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/permute/normal/NORMAL.PC30.All.thresholds.txt"

for i in {0..100}; do 
    wsbatch --partition=shared-cpu --mem=10G --time=12:00:00 -o ${OUT_DIR}/slurm-%j_${i}_100.out --wrap="${BIN} cis --vcf ${VCF} --bed ${BED} --normal --mapping ${THRESHOLD} --chunk ${i} 100 --out ${OUT_DIR}/NORMAL.conditional_${i}_100.txt.gz"; 
done 

VCF2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_genotypes/syscol_276_cancerNames_INFO_filtered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
BED2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/tumor/SYSCOL_tumor_276_gene_TE.FM05.resid.cpm.bed.gz" 
OUT_DIR2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/mapping/tumor"
THRESHOLD2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/permute/tumor/TUMOR.PC30.All.thresholds.txt"

for i in {0..100}; do 
    wsbatch --partition=shared-cpu --mem=10G --time=12:00:00 -o ${OUT_DIR2}/slurm-%j_${i}_100.out --wrap="${BIN} cis --vcf ${VCF2} --bed ${BED2} --normal --mapping ${THRESHOLD2} --chunk ${i} 100 --out ${OUT_DIR2}/TUMOR.conditional_${i}_100.txt.gz"; 
done

