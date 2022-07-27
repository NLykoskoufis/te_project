#!/bin/bash 
set -e 

DIR="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/teqtl/permute"
BIN="/srv/beegfs/scratch/users/l/lykoskou/TE/V3/te_project/common/runFDR_cis.R"

NORMAL_DIR="${DIR}/normal"
TUMOR_DIR="${DIR}/tumor"

## NORMAL 
zcat ${NORMAL_DIR}/*.txt.gz  | bgzip -c > ${NORMAL_DIR}/NORMAL_teqtls.chrALL.txt.gz
Rscript ${BIN} ${NORMAL_DIR}/NORMAL_teqtls.chrALL.txt.gz 0.05 ${NORMAL_DIR}/NORMAL_teqtls.chrALL >> ${NORMAL_DIR}/NORMAL_teqtls.chrALL.log

## TUMOR 
zcat ${TUMOR_DIR}/*.txt.gz | bgzip -c > ${TUMOR_DIR}/TUMOR_teqtls.chrALL.txt.gz 
Rscript ${BIN} ${TUMOR_DIR}/TUMOR_teqtls.chrALL.txt.gz 0.05 ${TUMOR_DIR}/TUMOR_teqtls.chrALL >> ${TUMOR_DIR}/TUMOR_teqtls.chrALL.log


echo "Done"
