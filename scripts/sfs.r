# Allele frequency distributions Sm single genotype project
# Author: Kathrin Jutzeler
# Date: May 15, 2023
# Last updated: February 8, 2024


# Run with R version 4.2.3
library(vcfR)
library(pegas)
library(ape)

library(fields)
library(tidyverse)


# Approach with pegas and vcf file 
# Define colors and order ####
setwd('./sfs')

lab_colors = c("OR" = "#B22222", "BRE" = "#00A08A", "LE" = "#F98400", "EG" = "#5BBCD6", "NMRI" = '#8A2BE2')
field_colors <- c('Brazil' = "#9986A5", 'Senegal' = "#79402E", 'Niger'="#CCBA72" ,'Tanzania'='#0F0D0E')

all_colors <- c(lab_colors, field_colors)

order <- c('BRE', 'EG', 'LE', 'OR', 'NMRI', 'Brazil', 'Niger', 'Senegal', 'Tanzania')

# Import files ####
# Note these files were generated from annotated_snps_nuclear_CDS_nosex_v10.vcf
files <- list.files()

# Write functions  ####
f_sfs <- function(x){
  vcf.in <- read.vcfR(x)  
  #myID <- getID(vcf.in)
  #vcf.in <- vcf.in[!duplicated(myID, incomparables = NA), ]
  
  vcf.bin <- vcfR2DNAbin(vcf.in)
}

# Process data ####

out <- lapply(files[2:10], f_sfs)

out2 <- lapply(out, function(x) site.spectrum(x))

names(out2) <-  c('Brazil', 'BRE', 'EG', 'LE', 'Niger', 'NMRI', 'OR', 'Senegal', 'Tanzania') 

out3 <- lapply(out2, function(x) as.vector(x))

out4 <- lapply(out3, function(x) data.frame(x))

df <- bind_rows(out4, .id = 'pop')

df <- df %>% mutate(origin = ifelse(pop %in% c('BRE', 'EG', 'LE', 'OR', 'NMRI'), 'Lab', 'Field'))

df$pop  <- factor(df$pop, levels = order)

df <- df %>%
  group_by(pop) %>%
  mutate(sample = seq(x))

df %>% group_by(pop) %>% summarize(max = max(x))
# A tibble: 9 × 2
#pop        max
#<fct>    <int>
#  1 BRE       2142
#2 EG        2207
#3 LE        1683
#4 OR        2158
#5 NMRI      1324
#6 Brazil    3885
#7 Niger     3433
#8 Senegal   4118
#9 Tanzania 11026

# Plot folded SFS ####
p_sfs_folded <- ggplot(df, aes(x = sample, y = x, fill = pop)) +
  geom_col() +
  facet_wrap(pop ~ ., scales = 'free') +
  scale_fill_manual(values = all_colors) +
  labs(x = 'Sample size', y = 'Number of sites') +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 18), axis.line = element_line(), legend.position = 'none')

ggsave('Plots/SFS_folded.png', width = 8, height = 6)

# Plot ECDF ####
p_ecdf <- ggplot(df, aes(x, color = origin)) + stat_ecdf() +
  #scale_color_manual(values = all_colors) + 
  labs(x = "Allele frequency", y = "ECDF") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 18), axis.line = element_line())

ggsave('ECDF_origin.png', width = 8, height = 6)

# STATS ####
combos <- combn(order, 2, simplify = F)

f_ks <- function(x) {
  a <- x[1]
  b <- x[2]
  
  a <- df %>% filter(pop == a)
  b <- df %>% filter(pop == b)
  
  return(ks.test(a$x, b$x))
}

field <- df %>% filter(pop != 'Tanzania', origin == 'Field')
lab <- df %>% filter(origin == 'Lab')

ks.test(field$x, lab$x)

ks_out <- lapply(combos, f_ks)

names(ks_out) <- combos

results <- data.frame(
  comparison = names(ks_out), 
  D = sapply(ks_out, "[[", "statistic"), 
  pvalue = sapply(ks_out, "[[", "p.value")
)

results <- arrange(results, pvalue)
write_delim(results, 'SFS_stats.txt')


# Allele frequency appraoch (generated with bcftools) ####
setwd('./AF')

freq_files <- list.files()

# Import data
freq <- lapply(freq_files, read.table, stringsAsFactors = F)

# Assign names
names(freq) <- c('Brazil', 'BRE', 'EG', 'LE', 'Niger', 'NMRI', 'OR', 'Senegal', 'Tanzania')

df_af <- bind_rows(freq, .id = 'pop')
df_af$pop <- factor(df_af$pop, levels = order)

# Plot AF distribution ####
p_af <- ggplot(df_af, aes(V3, fill = pop)) +
  geom_histogram() +
  facet_wrap(pop ~ ., scales = 'free') +
  scale_fill_manual(values = all_colors) +
  labs(x = 'Allele frequency', y = 'Number of sites') +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 18), axis.line = element_line(), legend.position = 'none')

setwd('..')
ggsave('Plots/AF_distribution.png')


# Plot ECDF of AF ####
df_af <- df_af %>% mutate(origin = ifelse(pop %in% c('BRE', 'LE', 'OR', 'EG', 'NMRI'), 'Lab', 'Field')) 

field <- filter(df_af, origin == 'Field')
lab <- filter(df_af, origin == 'Lab')

ks <- ks.test(field$V3, lab$V3) 

#Asymptotic two-sample Kolmogorov-Smirnov test

#data:  field$V3 and lab$V3
#D = 0.45971, p-value < 2.2e-16
#alternative hypothesis: two-sided

df_af$origin <- factor(df_af$origin, levels = c('Lab', 'Field'))

p_ECDF <- ggplot(df_af, aes(V3, color = pop, linetype = origin)) + 
  stat_ecdf() +
  scale_color_manual(values = all_colors) + 
  labs(x = "Allele frequency", y = "ECDF", color = 'Population', linetype = 'Origin') +
  guides(color = guide_legend(nrow = 2, byrow = T), linetype = guide_legend(nrow = 2)) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 16), axis.line = element_line(), legend.box.just = 'center', 
        legend.position = 'bottom') +
  annotate('text', x = 0.2, y =1, label = 'Kolmogorov-Smirnov test, p < 0.001')

ggsave('Plots/ECDF_AF.png', width = 8, height = 6)


# STATS
combos <- combn(order, 2, simplify = F)

f_ks <- function(x) {
  a <- x[1]
  b <- x[2]
  
  a <- df %>% filter(pop == a)
  b <- df %>% filter(pop == b)
  
  return(ks.test(a$V3, b$V3))
}

ks_out_af <- lapply(combos, f_ks)

names(ks_out_af) <- combos

results <- data.frame(
  comparison = names(ks_out_af), 
  D = sapply(ks_out_af, "[[", "statistic"), 
  pvalue = sapply(ks_out_af, "[[", "p.value")
)

results <- arrange(results, pvalue)
write_delim(results, 'AF_stats.txt')

# Calculate AF in bins ####
breaks <- seq( 0, 1e6, 2500)
  
af_bin <- df_af %>%
  group_by(pop) %>%
  reframe(out = stats.bin(V2, V3, breaks = breaks))

BRE <- af_bin %>% filter(pop == 'BRE')
EG<- af_bin %>% filter(pop == 'EG')

ks.test(BRE$out$stats['mean',], EG$out$stats['mean',])

t <- list(BRE$out$stats['mean',], EG$out$stats['mean',])

names(t) <- c('BRE', 'EG')

t <- lapply(t, data.frame)
t <- bind_rows(t, .id = 'pop')

ggplot(t, aes(X..i.., color = pop)) +
   stat_ecdf() # ~ ld_binned$centers,


# Export plots ####

pdf('Figure3_Tajima_ECDF.pdf', width = 8, height = 10)

(p_tajima / p_ECDF)  +
  plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold'))

dev.off()


pdf('FigureSx_sfs.pdf', width = 10, height = 10)

(p_sfs_folded / p_af) +
  plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold'))

dev.off()


