# Code used to calculate unlinked LD

pop=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ng" "sn")

for pop in ${pop[@]}; do
# Step 1: Extract header lines
zcat filtered_vars_${pop}.vcf.gz | grep '^#' > header_${pop}.vcf

# Step 2: Randomly sample 1000 data rows
zcat filtered_vars_${pop}.vcf.gz | grep -v '^#'  | shuf -n 1000 > sampled_variants_${pop}.vcf

# Step 3: Combine header with sampled data
cat header_${pop}.vcf > sampled_${pop}.vcf
cat sampled_variants_${pop}.vcf >> sampled_${pop}.vcf
done

pop=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ng" "sn")

for pop in ${pop[@]}; do
gatk SortVcf I=sampled_${pop}.vcf O=final_${pop}.vcf
done

# Prepare a chromosome size file (e.g. from .fai index)
cut -f1,2 /master/kbailey/egg_RNA/data/reference/schistosoma_mansoni.PRJEA36577.WBPS18.genomic.fa.fai > chrom_sizes.txt

# Compute cumulative offsets
awk 'BEGIN {OFS="\t"} {
  chr = $1
  len = $2
  offset[chr] = total
  total += len
}
END {
  for (chr in offset)
    print chr, offset[chr]
}' chrom_sizes.txt | sort -k1,1 > chrom_offsets.txt

pop=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ng" "sn")
for pop in ${pop[@]}; do

awk 'BEGIN {OFS="\t"} FNR==NR {
  offset[$1] = $2; next
}
/^#/ { print; next }
{
  chr = $1
  pos = $2
  gpos = offset[chr] + pos
  $2 = gpos  # Replace POS with genomic position
  print
}' chrom_offsets.txt final_${pop}.vcf > pos_${pop}.vcf
done

pop=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ng" "sn")
for pop in ${pop[@]}; do

awk '
# Skip header lines
/^#/ { next }

{
  chrom = $1
  pos = $2

  if (!(chrom in min) || pos < min[chrom]) min[chrom] = pos
  if (!(chrom in max) || pos > max[chrom]) max[chrom] = pos
}

END {
  for (c in min) {
    print c, min[c], max[c]
  }
}
' OFS='\t' pos_${pop}.vcf > chrom_ranges_${pop}.bed
done

pop=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ng" "sn")

for pop in ${pop[@]}; do
awk '
BEGIN { OFS = "\t" }
# Keep header lines
/^#/ { print; next }

# Replace CHROM column with "1"
{ $1 = "SM_V10_1"; print }

' chrom_ranges_${pop}.bed > chrom_pos_ranges_${pop}.bed
done

pop=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ng" "sn")

for pop in ${pop[@]}; do
awk '
BEGIN { OFS = "\t" }
# Keep header lines
/^#/ { print; next }

# Replace CHROM column with "1"
{ $1 = "SM_V10_1"; print }

' pos_${pop}.vcf > chr_pos_${pop}.vcf
done

# Calculate LD with plink
pop=("EG" "BRE" "LE" "OR" "NMRI" "br" "tz" "ng" "sn")
for pop in ${pop[@]}; do

plink \
--vcf chr_pos_${pop}.vcf \
            --out ${pop} \
            --double-id \
            --allow-extra-chr \
            --ld-window 999999 \
            --ld-window-kb 10000 \
            --ld-window-r2 0 \
            --r2
done
# The rest was done in R
