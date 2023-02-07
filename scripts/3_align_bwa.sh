echo """
#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N SAMPLE
#$ -o /master/kbailey/sm_single_gt/results/logs/align.SAMPLE.log
#$ -j y
#$ -q all.q
#$ -pe smp 4

""" >header.qsub.sh

mkdir /master/kbailey/sm_single_gt/results/aligned
mkdir /master/kbailey/sm_single_gt/results/nodups

REF_DIR="/master/kbailey/sm_single_gt/data/reference"
RESULTS_DIR="/master/kbailey/sm_single_gt/results"

for sample in $(cat /master/kbailey/sm_single_gt/data/seq_data/test_list.txt);
  do
  SAM=$(echo ${sample} | cut -f1 -d"_")

  CMD="bwa mem -M $REF_DIR/GCA_000237925.5_SM_V9_genomic.fna.gz $RESULTS_DIR/trimmed_fastq/${sample}_R1_001_val_1.fq.gz $RESULTS_DIR/trimmed_fastq/${sample}_R2_001_val_2.fq.gz | samtools sort -o \"$RESULTS_DIR/aligned/${sample}.bam\"
  picard MarkDuplicates I=\"$RESULTS_DIR/aligned/${sample}.bam\" O=\"$RESULTS_DIR/nodups/${sample}.bam\" M=\"${sample}.txt\" REMOVE_DUPLICATES=true ASSUME_SORTED=true
  samtools flagstat $RESULTS_DIR/nodups/${sample}.bam > ${sample}.flagstat;
  samtools index $RESULTS_DIR/nodups/${sample}.bam"

    cat header.qsub.sh >Sm_${SAM}.sh
    sed -i "s/SAMPLE/Sm_${SAM}/g" Sm_${SAM}.sh

    echo "$CMD" >>Sm_${SAM}.sh

    #qsub -V -cwd -S /bin/bash -N haplo_${sample}_${interval} -o /master/kbailey/ascaris_test/results/logs/haplo.o_${sample}_${interval}.log  ${sample}_${interval}.sh
    qsub Sm_${SAM}.sh
  done
