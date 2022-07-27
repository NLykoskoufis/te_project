#!/bin/bash 

BIN="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/binaries/data_qtls/teqtls/bn/5_bnlearn.call3.R"
NORMAL_TRIPLETS="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/triplet/normal"
TUMOR_TRIPLETS="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/triplet/tumor"

NORMAL_OUTPUT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/normal"
TUMOR_OUTPUT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/bn/tumor"


wsbatch --partition=shared-cpu -N1 -n1 -c1 --mem=10G --time=01:00:00 --wrap="Rscript ${BIN} ${NORMAL_TRIPLETS} ${NORMAL_OUTPUT}/NORMAL.bnlearn_results.txt"

wsbatch --partition=shared-cpu -N1 -n1 -c1 --mem=10G --time=01:00:00 --wrap="Rscript ${BIN} ${TUMOR_TRIPLETS} ${TUMOR_OUTPUT}/TUMOR.bnlearn_results.txt"

