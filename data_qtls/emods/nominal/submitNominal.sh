#!/bin/bash 

BIN="/srv/beegfs/scratch/shares/brauns_lab/bin/QTLtools"

# NORMAL 
OUT_DIR="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/emods/nominal/normal"
VCFT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/normal/SYSCOL_normal_275_TEs.FM05.resid.cpm.bed.gz"
BEDT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/normal/SYSCOL_normal_275_genes.FM05.resid.cpm.bed.gz"
THRESHOLD="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/emods/permute/normal/NORMAL.TE.GENE.permute.chrAll.thresholds.txt"

wsbatch --partition=shared-cpu --mem=20G --time=02:00:00 --wrap="${BIN} cis --vcf ${VCFT} --bed ${BEDT} --nominal ${THRESHOLD} --normal --out ${OUT_DIR}/NORMAL.TE.GENE.nominal.chrAll.txt.gz"


# TUMOR 

OUT_DIR2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/emods/nominal/tumor"
VCFT2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/tumor/SYSCOL_tumor_276_TEs.FM05.resid.cpm.bed.gz"
BEDT2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/tumor/SYSCOL_tumor_276_genes.FM05.resid.cpm.bed.gz"
THRESHOLD2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/emods/permute/tumor/TUMOR.TE.GENE.permute.chrAll.thresholds.txt"


wsbatch --partition=shared-cpu --mem=20G --time=02:00:00 --wrap="${BIN} cis --vcf ${VCFT2} --bed ${BEDT2} --nominal ${THRESHOLD2} --normal --out ${OUT_DIR2}/TUMOR.TE.GENE.nominal.chrAll.txt.gz"
