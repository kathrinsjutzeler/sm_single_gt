echo """
#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N SAMPLE
#$ -o /master/kbailey/sm_single_gt/results/logs/trim.SAMPLE.log
#$ -j y
#$ -q all.q
#$ -pe smp 10

""" >header.qsub.sh

mkdir /master/kbailey/sm_single_gt/results/trimmed_reports

REF_DIR="/master/kbailey/sm_single_gt/data/seq_data"
RESULTS_DIR="/master/kbailey/sm_single_gt/results"

for sample in $(cat /master/kbailey/sm_single_gt/data/seq_data/sample_list.txt);
  do
  SAM=$(echo ${sample} | cut -f1 -d"_")
  
  CMD="trim_galore -q 28 --fastqc_args \"--outdir $RESULTS_DIR/trimmed_reports\" --illumina --max_n 1 
--trim-n -o $RESULTS_DIR/trimmed_fastq --clip_R1 9 --clip_R2 9 --paired $REF_DIR/${sample}_R1_001.fastq.gz 
$REF_DIR/${sample}_R2_001.fastq.gz"

    cat header.qsub.sh >Sm_${SAM}.sh
    sed -i "s/SAMPLE/Sm_${SAM}/g" Sm_${SAM}.sh

    echo $CMD >>Sm_${SAM}.sh

    qsub Sm_${SAM}.sh
  done
