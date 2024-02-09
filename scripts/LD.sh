# Calculate LD stats in windows 
pop=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ne" "sn")

mkdir ~/sm_single_gt/results/ld_v10

for pop in ${pop[@]}; do

gatk SelectVariants -V /master/kbailey/sm_single_gt/results/filtered_v10/maf05_nuclear_CDS_v10.vcf \
--exclude-sample-name /master/kbailey/sm_single_gt/results/lists/outgroup.args \
--exclude-intervals SM_V10_Z \
--exclude-intervals SM_V10_WSR \
--allow-nonoverlapping-command-line-samples \
--sample-name /master/kbailey/sm_single_gt/results/lists/${pop}.args \
-O /master/kbailey/sm_single_gt/results/ld_v10/maf05_CDS_nosex_${pop}.vcf \
-R /master/kbailey/egg_RNA/data/reference/schistosoma_mansoni.PRJEA36577.WBPS18.genomic.fa ;
done

pop=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ne" "sn")
for pop in ${pop[@]}; do
plink \
            --threads 6 \
            --r2 \
            --vcf /master/kbailey/sm_single_gt/results/ld_v10/maf05_CDS_nosex_${pop}.vcf \
            --out ~/sm_single_gt/results/ld_v10/maf05_CDS_nosex_${pop}.ld \
            --double-id \
            --allow-extra-chr \
            --ld-window-r2 0.0 \
            --ld-window 1000000 \
            --ld-window-kb 1000
done
