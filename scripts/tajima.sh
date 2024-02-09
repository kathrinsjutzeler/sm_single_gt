#!/bin/bash

pop=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ne" "sn")

mkdir ~/sm_single_gt/results/tajima_v10

for pop in ${pop[@]}; do

vcftools --vcf annotated_snps_nuclear_CDS_v10.vcf --TajimaD 25000 --keep /master/kbailey/sm_single_gt/results/lists/${pop}.args --out /master/kbailey/sm_single_gt/results/tajima_v10/annotated_snps_nuclear_CDS_tajima_${pop};

done 
