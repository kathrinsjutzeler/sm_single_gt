# This part will set up the qsub properties and save them in a new script

echo """
#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N NAME
#$ -o /master/kbailey/sm_single_gt/results/logs/gdbi.NAME.log
#$ -j y
#$ -q all.q
#$ -pe smp 10

""" >header.qsub.sh

# This is where we set up new folders and variables
mkdir /master/kbailey/sm_single_gt/results/gdbsingle_gt

REF=/master/kbailey/sm_single_gt/data/reference/

# This is where we write the code that will be run. Note that the main command is stored as a variable.
for interval in $(cat $REF/intervals/intervals.list); do
  NAME=$(echo ${interval} | cut -f1 -d"-")

  CMD="gatk --java-options \"-Xmx12g -Xms12g\" GenomicsDBImport \
      --genomicsdb-workspace-path /master/kbailey/sm_single_gt/results/gdbsingle_gt/${NAME} \
      --batch-size 50 \
      --reader-threads 5 \
      --sample-name-map /master/kbailey/sm_single_gt/results/haplo_merged/all_samples.sample_map \
      --tmp-dir /master/kbailey/sm_single_gt/temp \
      -L $REF/intervals/${interval}"


# Save the header as a script for each individual loop
      cat header.qsub.sh >${NAME}.sh
# Modify the header to change the name of the job and log file
      sed -i "s/NAME/S${NAME}/g" ${NAME}.sh
# Append the command to the new script
      echo $CMD >>${NAME}.sh
# Execute the script on the server
      qsub ${NAME}.sh
done
