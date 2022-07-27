#!/bin/bash 

NORMAL="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/optimize/normal"
TUMOR="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/optimize/tumor"

BIN="/srv/beegfs/scratch/users/l/lykoskou/TE/V3/te_project/common/runFDR_cis.R"

for i in 0 1 2 5 10 20 30 40 50 60 70 80 90 100; do
    # NORMAL
    echo "${i}"
    zcat ${NORMAL}/NORMAL.PC${i}.*.txt.gz > ${NORMAL}/NORMAL.PC${i}.All.txt 
    Rscript ${BIN} ${NORMAL}/NORMAL.PC${i}.All.txt 0.05 ${NORMAL}/NORMAL.PC${i}.All >> ${NORMAL}/NORMAL.PC${i}.log
    # TUMOR
    zcat ${TUMOR}/TUMOR.PC${i}.*.txt.gz > ${TUMOR}/TUMOR.PC${i}.All.txt
    Rscript ${BIN} ${TUMOR}/TUMOR.PC${i}.All.txt 0.05 ${TUMOR}/TUMOR.PC${i}.All >> ${TUMOR}/TUMOR.PC${i}.log
done 

wc -l ${NORMAL}/NORMAL.PC*.All.significant.txt | tr -s " " "\t" | cut -f2,3 | sed 's/NORMAL.PC\|.All.significant.txt//g' | awk '{print $2"\t"$1"\tnormal"}' | sort -V -k1,1 > ${NORMAL}/NORMAL.Neqtls.PC.txt

wc -l ${TUMOR}/*.significant.txt | tr -s " " "\t" | cut -f2,3 | sed 's/TUMOR.PC\|.All.significant.txt//g' | awk '{print $2"\t"$1"\ttumor"}' | sort -V -k1,1 > ${TUMOR}/TUMOR.Neqtls.PC.txt
