#!/bin/bash 

BIN="/srv/beegfs/scratch/shares/brauns_lab/bin/QTLtools"

# NORMAL 
OUT_DIR="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/emods/permute/normal"
VCFT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/normal/SYSCOL_normal_275_TEs.FM05.resid.cpm.bed.gz"
BEDT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/normal/SYSCOL_normal_275_genes.FM05.resid.cpm.bed.gz"

#for c in {0..100}; do
#    wsbatch --partition=shared-cpu --time=05:00:00 -o ${OUT_DIR}/slurm-%j_${c}_100.out --mem=10G --wrap="${BIN} cis --vcf ${VCFT} --bed ${BEDT} --normal --permute 1000 --out ${OUT_DIR}/NORMAL.TE.GENE.permute_${c}_100.txt.gz --chunk ${c} 100"; 
#done

# TUMOR 

OUT_DIR2="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/emods/permute/tumor"
VCFT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/tumor/SYSCOL_tumor_276_TEs.FM05.resid.cpm.bed.gz"
BEDT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_expression/tumor/SYSCOL_tumor_276_genes.FM05.resid.cpm.bed.gz"

for c in {0..100}; do
    wsbatch --partition=shared-cpu --time=05:00:00 -o ${OUT_DIR2}/slurm-%j_${c}_100.out --mem=10G --wrap="${BIN} cis --vcf ${VCFT} --bed ${BEDT} --normal --permute 1000 --out ${OUT_DIR2}/TUMOR.TE.GENE.permute_${c}_100.txt.gz --chunk ${c} 100"; 
done


