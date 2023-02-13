# This part will set up the qsub properties and save them in a new script

echo """
#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N NAME
#$ -o /master/kbailey/sm_single_gt/results/logs/gvcf.NAME.log
#$ -j y
#$ -q all.q
#$ -pe smp 10

""" >header.qsub.sh

# This is where we set up new folders and variables
mkdir /master/kbailey/sm_single_gt/results/Genotype

ref=/master/kbailey/sm_single_gt/data/reference

# This is where we write the code that will be run. Note that the main command is stored as a variable.
for interval in $(cat $ref/intervals/intervals.list); do
  NAME=$(echo ${interval} | cut -f1 -d"-")

  CMD="gatk --java-options \"-Xmx12g\" GenotypeGVCFs \
        -R $ref/GCA_000237925.5_SM_V9_genomic.fna \
        -V gendb://$master/master/kbailey/sm_single_gt/results/gdbsingle_gt/${NAME} \
        -O /master/kbailey/sm_single_gt/results/Genotype/${NAME}.vcf.gz \
        --tmp-dir /master/kbailey/sm_single_gt/temp"

# Save the header as a script for each individual loop
      cat header.qsub.sh >Sm_${NAME}gvcf.sh
# Modify the header to change the name of the job and log file
      sed -i "s/NAME/Sm_${NAME}/g" Sm_${NAME}gvcf.sh
# Append the command to the new script
      echo $CMD >>Sm_${NAME}gvcf.sh
# Execute the script on the server
      qsub Sm_${NAME}gvcf.sh
done
