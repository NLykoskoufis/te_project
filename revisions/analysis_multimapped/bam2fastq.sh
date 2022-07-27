#!/bin/bash 


BIN="/home/lykoskou/bin/picard.jar"

RAW_DIR="/scratch/lykoskou/te_project/data_syscol"
BAM_FILE_OLD_NAME="accepted_hits.bam" # all bam files are named like this...

DESIRED_READ_LENGTH=49

for DIR in ${RAW_DIR}/*;
do 
    SAMPLE=$(basename ${DIR})
    
    
    
    SAMPLE_ID=$(basename ${SAMPLE} | cut -d"_" -f1)
    BAM_FILE=$(ls ${RAW_DIR}/${SAMPLE}/*.bam)
    READ_LENGTH=$(( $(samtools view ${BAM_FILE} | head -1 | cut -f10 | wc -m) -1 ))
    
    if [[ ${READ_LENGTH} != 49 ]]; 
    then 
        echo $SAMPLE
        READ_TRIM_COUNT=$((${READ_LENGTH} - ${DESIRED_READ_LENGTH}))
        #echo ${READ_TRIM_COUNT}
        #wsbatch --partition=serial --mem=20G --time=03:00:00 -N1 -n1 -c1 --wrap="java -jar ${BIN} SamToFastq --I ${BAM_FILE} --FASTQ ${SAMPLE_ID}_1.fastq.gz --SECOND_END_FASTQ ${SAMPLE_ID}_2.fastq.gz --READ1_TRIM ${READ_TRIM_COUNT} --READ2_TRIM ${READ_TRIM_COUNT}"
        wsbatch --partition=serial --mem=20G --time=03:00:00 -N1 -n1 -c1 --wrap="java -jar ${BIN} SamToFastq --I ${BAM_FILE} --FASTQ ${DIR}/${SAMPLE_ID}_1.fastq.gz --SECOND_END_FASTQ ${DIR}/${SAMPLE_ID}_2.fastq.gz --READ1_TRIM ${READ_TRIM_COUNT} --READ2_TRIM ${READ_TRIM_COUNT}"
    fi
done

