#!/bin/bash

BIN="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/binaries/data_qtls/eqtl/specific/Check_opposite_cellType_for_nonSignificance.py"

DIR="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/specific"
EQTLN="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/interaction/normal/normal_specific_eQTL_analysis_results_ALL_conditional.txt.gz"
EQTLT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/interaction/tumor/tumor_specific_eQTL_analysis_results_ALL_conditional.txt.gz"

NOMN="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/nominal/NORMAL.All.txt.gz"
NOMT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/nominal/TUMOR.All.txt.gz"

CONDN="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/mapping/normal/NORMAL.conditional_bestPerRank.txt.gz"
CONDT="/srv/beegfs/scratch/users/l/lykoskou/TE/V4/data_qtl/eqtl/mapping/tumor/TUMOR.conditional_bestPerRank.txt.gz"

# TUMOR
wsbatch --partition=shared-cpu --time=12:00:00 --mem=10G --wrap="python3 ${BIN} ${EQTLT} ${NOMN} ${CONDT} ${DIR}/TUMOR.specificity_All.txt" 

# NORMAL 
wsbatch --partition=shared-cpu --time=12:00:00 --mem=10G --wrap="python3 ${BIN} ${EQTLN} ${NOMT} ${CONDN} ${DIR}/NORMAL.specificity_All.txt"

 
