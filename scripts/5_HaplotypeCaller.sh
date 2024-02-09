echo """
#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N As_SAMPLE
#$ -o /master/kbailey/sm_single_gt/results/logs/haplo.SAMPLE_INTERVAL.log
#$ -j y
#$ -q all.q
#$ -pe smp 1

""" >header.qsub.sh

REF_DIR="/master/kbailey/sm_single_gt/data/reference"
RESULTS_DIR="/master/kbailey/sm_single_gt/results"

mkdir /master/kbailey/sm_single_gt/results/haplo
mkdir /master/kbailey/temp/haplo

for interval in $(cat $REF_DIR/intervals/intervals.list); do
  for sample in $(cat /master/kbailey/sm_single_gt/data/part_samples.txt); do
    NUM=$(echo ${interval} | cut -f1 -d"-")
    SAM=$(echo ${sample} | cut -f1 -d"-")
	
    CMD="gatk --java-options \"-Xmx12g\" HaplotypeCaller \
    -R $REF_DIR/GCA_013433145.1_ASM1343314v1_genomic.fna \
    -I \"$RESULTS_DIR/finalbam/${sample}.bam\" \
    -O \"$RESULTS_DIR/haplo/${sample}.${NUM}.g.vcf.gz\" \
    -ERC GVCF \
    -L $REF_DIR/intervals/${interval} \
    --tmp-dir /master/kbailey/temp/haplo"

    cat header.qsub.sh >${SAM}_$NUM.sh
    sed -i "s/As_SAMPLE/As_${SAM}/g" ${SAM}_${NUM}.sh
    sed -i "s/SAMPLE/${SAM}/g" ${SAM}_${NUM}.sh
    sed -i "s/INTERVAL/$NUM/g" ${SAM}_${NUM}.sh

    echo $CMD >>${SAM}_${NUM}.sh

    #qsub -V -cwd -S /bin/bash -N haplo_${sample}_${interval} -o /master/kbailey/sm_single_gt/results/logs/haplo.o_${sample}_${interval}.log  ${sample}_${interval}.sh
    qsub ${SAM}_${NUM}.sh
  done
done
