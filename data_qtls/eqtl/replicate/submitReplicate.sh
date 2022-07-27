#!/bin/bash 
set -e 

BIN="/srv/beegfs/scratch/shares/brauns_lab/bin/QTLtools"

#### NORMAL ####
VCF="/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_genotypes/syscol_275_normalNames_INFOfiltered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
BED="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/normal/SYSCOL_normal_275_gene_TE.FM05.resid.cpm.bed.gz"
QTL="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/replicate/normal/TUMOR_eqtls.txt"
OUT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/replicate/normal/TUMOR_eqtls_replicated_in_NORMAL.txt.gz"

#wsbatch --partition=shared-cpu --mem=30G --time=02:00:00 -o /srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/replicate/normal/slurm-%j.out  --wrap="${BIN} rep --vcf ${VCF} --bed ${BED} --qtl ${QTL} --out ${OUT} --normal"




VCF2="/srv/beegfs/scratch/users/l/lykoskou/TE/V3/data_genotypes/syscol_276_cancerNames_INFO_filtered_NEW_maf5_callrate95_threshold7_hwe6.vcf.gz"
BED2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/tumor/SYSCOL_tumor_276_gene_TE.FM05.resid.cpm.bed.gz"
QTL2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/replicate/tumor/NORMAL_eqtls.txt"
OUT2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/replicate/tumor/NORMAL_eqtls_replicated_in_TUMOR.txt.gz"

wsbatch --partition=shared-cpu --mem=30G --time=02:00:00 -o /srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/replicate/tumor/slurm-%j.out  --wrap="${BIN} rep --vcf ${VCF2} --bed ${BED2} --qtl ${QTL2} --out ${OUT2} --normal"
