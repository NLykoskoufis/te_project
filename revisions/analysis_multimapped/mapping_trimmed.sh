#!/bin/bash 

module load hisat2/2.1.0

RAW_DIR="/scratch/lykoskou/te_project/data_syscol"

TRIMMED_DIR=($(dirname ${RAW_DIR}/*/*.fastq.gz | sort | uniq))

INDEX="/scratch/lykoskou/te_project/hg19_hs2/hg19"

for dir in ${TRIMMED_DIR[@]}; 
do 
    SAMPLE_ID=$(basename ${dir} | cut -d"_" -f1 )
    FASTQ1=${dir}/${SAMPLE_ID}_1.fastq.gz 
    FASTQ2=${dir}/${SAMPLE_ID}_2.fastq.gz 
    OUTPUT=${dir}/${SAMPLE_ID}.accepted_hits.bam

    #mv ${dir}/accepted_hits.bam ${dir}/accepted_hits.bam.tmp
    wsbatch --partition=serial --mem=50G --time=03:00:00 -N1 -n1 -c 6 --wrap="hisat2 -k 5 --seed 42 -p 6 -x ${INDEX} -1 ${FASTQ1} -2 ${FASTQ2} | samtools view -bSh - > ${OUTPUT}"
    #echo "hisat2 -k 5 --seed 42 -p 6 -x ${INDEX} -1 ${FASTQ1} -2 ${FASTQ2} | sambools view -bSh - -@ 6> ${OUTPUT}"
done 