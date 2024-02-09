gatk --java-options \"-Xmx1024g\" GenotypeGVCFs \
        -R $ref/GCA_000237925.5_SM_V9_genomic.fna \
        -V gendb://$master/master/kbailey/sm_single_gt/previous_results2/gdbimport/${NAME} \
        --include-non-variant-sites \
        -O /master/kbailey/sm_single_gt/results/allsites/${NAME}.vcf.gz \
        --tmp-dir /master/kbailey/sm_single_gt/temp -L $ref/intervals/${interval}

gatk --java-options "-Xmx1024g" LiftoverVcf \
     -I /master/kbailey/sm_single_gt/results/allsites/all_sites_merged_v9.vcf \
     -O /master/kbailey/sm_single_gt/results/allsites/all_sites_merged_v10.vcf \
     -CHAIN /master/kbailey/smv9tosmv10.chain \
     -REJECT rejected_variants.vcf \
      --TMP_DIR /master/kbailey/sm_single_gt/temp \
     -R /master/kbailey/egg_RNA/data/reference/schistosoma_mansoni.PRJEA36577.WBPS18.genomic.fa

bedtools intersect -a all_sites_merged_v10.vcf -b ~/sm_single_gt/data/reference/CDS_v10.gff -u -wa -header > all_sites_CDS_v10.vcf

gatk IndexFeatureFile all_sites_CDS_v10.vcf

gatk SelectVariants -V /master/kbailey/sm_single_gt/results/allsites/all_sites_CDS_v10.vcf \
--sample-name samples.args \
--exclude-intervals SM_V10_WSR \
--exclude-intervals SM_V10_Z \
--exclude-intervals SM_V10_MITO \
-O /master/kbailey/sm_single_gt/results/allsites/all_sites_CDS_v10_for_pixy.vcf \
-R /master/kbailey/egg_RNA/data/reference/schistosoma_mansoni.PRJEA36577.WBPS18.genomic.fa

bgzip all_sites_CDS_v10_for_pixy.vcf && tabix all_sites_CDS_v10_for_pixy.vcf.gz
conda activate pixy
    pixy --stats pi \
         --vcf /master/kbailey/sm_single_gt/results/allsites/all_sites_CDS_v10_for_pixy.vcf.gz \
         --populations /master/kbailey/sm_single_gt/results/allsites/pops.args \
         --window_size 25000 \
         --n_cores 192 \
         --output_prefix single_gt \
         --output_folder /master/kbailey/sm_single_gt/results/allsites/pi
    conda deactivate
