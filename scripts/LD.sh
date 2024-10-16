#!/bin/bash

# Calculate LD stats in windows 
for pop in EG NMRI BRE LE OR; do 

gatk SelectVariants -V /master/kbailey/sm_single_gt/results/filtered/annotated_snps_${pop}.vcf.gz \
    --exclude-intervals SM_V10_WSR \
    --exclude-intervals SM_V10_Z \
    --exclude-intervals SM_V10_MITO \
    --exclude-non-variants \
    -O /master/kbailey/sm_single_gt/results/filtered/annotated_snps_autosomal_${pop}.vcf.gz;
done

#MAF filter
for pop in EG NMRI BRE LE OR; do 
vcftools \
  --gzvcf /master/kbailey/sm_single_gt/results/filtered_lab/annotated_snps_autosomal_${pop}.vcf.gz \
  --maf 0.05 `# Rare variants` \
  --recode \
  --recode-INFO-all \
  --stdout \
  >/master/kbailey/sm_single_gt/results/filtered_lab/maf05_autosomal_${pop}.vcf ;
done

for pop in EG NMRI BRE LE OR; do 
plink \
            --threads 6 \
            --r2 \
            --vcf master/kbailey/sm_single_gt/results/sfs/maf05_CDS_autosome_${pop}_v10.vcf \
            --out ~/sm_single_gt/results/ld/maf05_${pop}.ld \
            --double-id \
            --allow-extra-chr \
            --ld-window-r2 0.0 \
            --ld-window 1000000 \
            --ld-window-kb 1000
done
