#!/bin/bash

mkdir /master/kbailey/sm_single_gt/results/filtered

RESULTS_DIR="/master/kbailey/sm_single_gt/results"

# This tool is designed for hard-filtering variant calls based on certain criteria. https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
gatk VariantFiltration \
  -V $RESULTS_DIR/merged_v10.vcf \
  -O $RESULTS_DIR/filtered_v10/merged_soft.vcf.gz \
  -R /master/kbailey/egg_RNA/data/reference/schistosoma_mansoni.PRJEA36577.WBPS18.genomic.fa \
  --filter-expression "FS > 60.0" `# FisherStrand: Phred-scaled probability, no strand bias = close to 0` \
  --filter-name "fs_lt_60" \
  --filter-expression "SOR > 3.0" `# StrandOddsRatio: This is another way to estimate strand bias using a test similar to the symmetric odds ratio test` \
  --filter-name "sor_lt_3" \
  --filter-expression "MQ < 40.0" `# RMSMapping Quality: This is the root mean square mapping quality over all the reads at the site` \
  --filter-name "mq_gt_40" \
  --filter-expression "MQRankSum < -12.5" `# This is the u-based z-approximation from the Rank Sum Test for mapping qualities` \
  --filter-name "MQRankSum_lt_-12.5" \
  --filter-expression "ReadPosRankSum < -8.0" `# It compares whether the positions of the reference and alternate alleles are different within the reads` \
  --filter-name "ReadPosRankSum_lt_-8" \
  --filter-expression "QD < 2.0" `#  QualByDepth: This annotation is intended to normalize the variant quality in order to avoid inflation caused when there is deep coverage.` \
  --filter-name "qd_lt_2"

# This tool makes it possible to select a subset of variants based on various criteria in order to facilitate certain analyses
# SNPs
gatk SelectVariants \
  -V $RESULTS_DIR/filtered_v10/merged_soft.vcf.gz --exclude-filtered -O $RESULTS_DIR/filtered_v10/merged_hard.vcf.gz -R /master/kbailey/egg_RNA/data/reference/schistosoma_mansoni.PRJEA36577.WBPS18.genomic.fa

# Use VCFtools for more rigourous filtering
#get bialleleic snps with min gq and depth
vcftools \
	--gzvcf $RESULTS_DIR/filtered_v10/merged_hard.vcf.gz \
	--minGQ 15 `# Exclude all genotypes with a quality below the threshold specified` \
	--minDP 10 `# DP is the filtered depth, at the sample level` \
  --min-alleles 2 `# include only biallelic sites` \
	--max-alleles 2 \
	--remove-indels `# remove sites that contain an indel` \
	--recode `# generates a new file ... ` \
	--recode-INFO-all `#... and keeps the INFO field` \
	--stdout `# pipe into command or new filename `\
	>$RESULTS_DIR/filtered_v10/merged_bisnp_v10.vcf

#remove sites where missingness is > 0.58 (3 step process)
vcftools --vcf $RESULTS_DIR/filtered_v10/merged_bisnp_v10.vcf \
  --missing-site `#stats table to count missing sites` \
  --stdout \
	>$RESULTS_DIR/filtered_v10/missing_per_site.tbl

# Filter missing sites 
  cat $RESULTS_DIR/filtered_v10/missing_per_site.tbl \
  	| sed 1d `# removed header line` \
  	| awk '{if ($6<=0.58) print $1"\t"$2}' `# keep lines with <=58% missingness and print out columns 1 (CHR) and 2 (POS)`\
  	>$RESULTS_DIR/filtered_v10/high_freq_gt_sites.list

  # Generate new working file excluding sites with >=58% missingness
  vcftools \
	   --vcf $RESULTS_DIR/filtered_v10/merged_bisnp_v10.vcf \
	   --positions $RESULTS_DIR/filtered_v10/high_freq_gt_sites.list `# only keep sites with <50% missingness`\
	   --recode \
	   --recode-INFO-all \
	   --stdout \
	   >$RESULTS_DIR/filtered_v10/high_freq_gt_sites_v10.vcf 

#remove individuals genotyped at less than 50% of sites (3 step process as before)
vcftools \
  --vcf $RESULTS_DIR/filtered_v10/high_freq_gt_sites_v10.vcf \
  --missing-indv \
  --stdout \
  >$RESULTS_DIR/filtered_v10/indv_gt_freq.tbl

  cat $RESULTS_DIR/filtered_v10/indv_gt_freq.tbl \
    | sed 1d \
    | awk '{if ($5<=0.5) print $1}' \
    >$RESULTS_DIR/filtered_v10/indvs_to_keep.list

  vcftools \
    --vcf $RESULTS_DIR/filtered_v10/high_freq_gt_sites_v10.vcf \
    --keep $RESULTS_DIR/filtered_v10/indvs_to_keep.list   \
    --recode \
    --recode-INFO-all \
    --stdout \
    >$RESULTS_DIR/filtered_v10/indv_and_site_filt_v10.vcf

#annotate remaining snps
bcftools annotate \
  --set-id +'%CHROM\:%POS' \
  $RESULTS_DIR/filtered_v10/indv_and_site_filt_v10.vcf \
  >$RESULTS_DIR/filtered_v10/annotated_snps_v10.vcf

#get maf05s (Minor Allele Frequency)
vcftools \
  --vcf $RESULTS_DIR/filtered_v10/annotated_snps_v10.vcf \
  --maf 0.05 `# Rare variants` \
  --recode \
  --recode-INFO-all \
  --stdout \
  >$RESULTS_DIR/filtered_v10/maf05_v10.vcf
