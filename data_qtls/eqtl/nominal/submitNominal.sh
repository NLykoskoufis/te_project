#!/bin/bash 
set -e

BIN="/srv/beegfs/scratch/shares/brauns_lab/bin/QTLtools"

VCF="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_genotypes/syscol_275_normalNames_INFOfiltered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
BED="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/normal/SYSCOL_normal_275_gene_TE.FM05.resid.cpm.bed.gz"
OUT_DIR="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/nominal/NORMAL.All.txt.gz"

wsbatch --partition=shared-cpu --mem=20G --time=05:00:00 --wrap="${BIN} cis --vcf ${VCF} --bed ${BED} --nominal 1 --out ${OUT_DIR} --normal"

VCF2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_genotypes/syscol_276_cancerNames_INFO_filtered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
BED2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/tumor/SYSCOL_tumor_276_gene_TE.FM05.resid.cpm.bed.gz" 
OUT_DIR2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/nominal/TUMOR.All.txt.gz"


wsbatch --partition=shared-cpu --mem=20G --time=05:00:00 --wrap="${BIN} cis --vcf ${VCF2} --bed ${BED2} --nominal 1 --out ${OUT_DIR2} --normal"
