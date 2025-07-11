# Filtering for high quality variants

gatk IndexFeatureFile -I probe_filtered_vars.vcf

# Exclude mito an sex chromosomes
pop=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ng" "sn")

for pop in ${pop[@]}; do

gatk SelectVariants -V /master/kbailey/sm_single_gt/results/post_review/probe_filtered_vars.vcf \
--sample-name /master/kbailey/sm_single_gt/results/lists/${pop}.args \
--exclude-intervals SM_V10_Z \
--exclude-intervals SM_V10_WSR \
--exclude-intervals SM_V10_MITO \
--allow-nonoverlapping-command-line-samples \
-O /master/kbailey/sm_single_gt/results/post_review/probe_filtered_invars_autosome_${pop}.vcf \
-R /master/kbailey/egg_RNA/data/reference/schistosoma_mansoni.PRJEA36577.WBPS18.genomic.fa
done

# Only retain high quality gt sites
pop=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ng" "sn")

for pop in ${pop[@]}; do

vcftools \
	   --vcf /master/kbailey/sm_single_gt/results/post_review/probe_filtered_invars_autosome_${pop}.vcf \
	   --max-missing 0.8 \
           --recode \
	   --recode-INFO-all \
	   --stdout \
	   >/master/kbailey/sm_single_gt/results/post_review/high_qual_sites${pop}.vcf
 done


#Sanity check

#!/bin/bash

# Define population codes
populations=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ng" "sn")

# Base paths
VCF_DIR="/master/kbailey/sm_single_gt/results/post_review"
OUTPUT_DIR="/master/kbailey/sm_single_gt/results/post_review"

# Loop through each population
for pop in "${populations[@]}"; do
  echo "Processing population: $pop"

  input_vcf="${VCF_DIR}/high_qual_sites${pop}.vcf"
  site_missing_out="${OUTPUT_DIR}/missing_per_site_80${pop}.tbl"
  indv_missing_out="${OUTPUT_DIR}/gt_freq_indv_${pop}.tbl"

  # Generate missing site statistics
  vcftools --vcf "$input_vcf" \
    --missing-site \
    --stdout \
    > "$site_missing_out"

  echo "Finished processing $pop"
done

echo "All populations processed."
#Looks good.

#Next filter out individuals
#!/bin/bash

# Define population codes
populations=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ng" "sn")

# Base paths
VCF_DIR="/master/kbailey/sm_single_gt/results/post_review"
OUTPUT_DIR="/master/kbailey/sm_single_gt/results/post_review"

# Loop through each population
for pop in "${populations[@]}"; do
  echo "Processing population: $pop"

  input_vcf="${VCF_DIR}/high_qual_sites${pop}.vcf"
  site_missing_out="${OUTPUT_DIR}/missing_per_site_80${pop}.tbl"
  indv_missing_out="${OUTPUT_DIR}/gt_freq_indv_80${pop}.tbl"

  # Generate missing site statistics
  vcftools --vcf "$input_vcf" \
    --missing-indv \
    --stdout \
    > "$indv_missing_out"

  echo "Finished processing $pop"
done

echo "All populations processed."

#Get a list with indivdual GT frequency > 80%
populations=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ng" "sn")

VCF_DIR="/master/kbailey/sm_single_gt/results/post_review"
OUTPUT_DIR="/master/kbailey/sm_single_gt/results/post_review"

for pop in "${populations[@]}"; do
cat ${OUTPUT_DIR}/gt_freq_indv_80${pop}.tbl \
    | sed 1d \
    | awk '{if ($5<=0.2) print $1}' \
    >${OUTPUT_DIR}/indvs_to_keep_${pop}.list
 done
Filter high quality sites to only keep these individuals:
populations=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ng" "sn")

VCF_DIR="/master/kbailey/sm_single_gt/results/post_review"
OUTPUT_DIR="/master/kbailey/sm_single_gt/results/post_review"

for pop in "${populations[@]}"; do
input_vcf="${VCF_DIR}/high_qual_sites${pop}.vcf"

vcftools \
    --vcf "$input_vcf" \
    --keep ${OUTPUT_DIR}/indvs_to_keep_${pop}.list \
    --recode \
    --recode-INFO-all \
    --stdout \
    >${OUTPUT_DIR}/indv_and_site_filt_${pop}.vcf
done

populations=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ng" "sn")
for pop in "${populations[@]}"; do
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' indv_and_site_filt_${pop}.vcf > GT_${pop}.tbl
done


#Make a list in R -> 131,207 positions
populations=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ng" "sn")
for pop in "${populations[@]}"; do
bgzip indv_and_site_filt_${pop}.vcf && tabix indv_and_site_filt_${pop}.vcf.gz
done
populations=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ng" "sn")
for pop in "${populations[@]}"; do
bcftools view -R loci.tsv indv_and_site_filt_${pop}.vcf.gz -o high_qual_filter/filtered_${pop}.vcf
done

#gatk IndexFeatureFile -I probe_filtered_vars.vcf

pop=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ng" "sn")

for pop in ${pop[@]}; do

gatk SelectVariants -V /master/kbailey/sm_single_gt/results/post_review/high_qual_filter/filtered_${pop}.vcf \
--exclude-non-variants \
-O /master/kbailey/sm_single_gt/results/post_review/high_qual_filter/filtered_vars_${pop}.vcf \
-R /master/kbailey/egg_RNA/data/reference/schistosoma_mansoni.PRJEA36577.WBPS18.genomic.fa
done
