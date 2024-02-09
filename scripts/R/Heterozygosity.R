# Heterozygosity  for Aim 1 - Sm Single GT project
# Author: Kathrin Jutzeler
# Date: November 16, 2023
# Last updated: February 8, 2024

library(ggridges)
library(tidyverse)
library(SNPRelate)
library(ggpubr)
library(lemon)
library(ggprism)
library(rstatix)


# Define colors
lab_colors = c("OR" = "#B22222", "BRE" = "#00A08A", "LE" = "#F98400", "EG" = "#5BBCD6", "NMRI" = '#8A2BE2')
field_colors <- c('Brazil' = "#9986A5", 'Senegal' = "#79402E", 'Niger'="#CCBA72" ,'Tanzania'='#0F0D0E')

all_colors <- c(lab_colors, field_colors)

order <- c('BRE', 'EG', 'LE', 'OR', 'NMRI', 'Brazil', 'Niger', 'Senegal', 'Tanzania')

############ Heterozygous SNPS#########
# Import the data ####

#anno_het <- read_delim("annotated_snps_nuclear_CDS_v10.het")
anno_het <- read_delim('annotated_snps_nuclear_nosex_CDS_v10.het')

df <- anno_het

# Compute the number of actual heterozygous sites (because the het file contains the observed (O) and expected (E) homozygous (HOM) sites)

df <- df %>%
  mutate(htz = N_SITES - `O(HOM)`) %>%
  mutate(adjusted = htz / 14e6 * 10e3)

# Import metadata
metadata <- read.csv('metadata.csv')

data <- df %>%
  left_join(metadata, by=c('INDV' = 'sample_id')) %>%
  #filter(!population %in% c('Uganda', 'Kenya', 'Caribbean', NA))
  #mutate(origin = ifelse(pop %in% c('OR', 'LE', 'EG', 'BRE'), 'Lab', "Field")) %>%
  filter(pop != 'Outgroup')

data %>%
  group_by(pop) %>%
  count()

# Make a summary table showing mean and sd
sumtable <-  data %>% group_by(pop, origin) %>% summarize(mean = mean(adjusted), sd = sd(adjusted))

#Stats ####
sumtable %>%
  group_by(origin) %>%
  shapiro_test(mean)

het_res <- sumtable %>% ungroup() %>% t_test(mean ~ origin) %>% mutate(pop = 'EG', origin = 'Lab')

# A tibble: 1 × 8
#.y.   group1 group2    n1    n2 statistic    df     p
#* <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl>
#  1 mean  Field  Lab        4     5      1.72  3.38 0.174

sumtable$pop <- factor(sumtable$pop, levels = order)

# Plot ####
p_het <- ggplot(sumtable, aes(pop, mean)) +
  geom_errorbar(aes(x=pop, ymin = mean-sd, ymax = mean+sd), width = 0.5) +
  facet_grid(~factor(origin, levels =c("Lab", "Field")), scales = 'free') +
  geom_col(aes(fill = pop)) +
  scale_fill_manual(values = all_colors) +
  theme_minimal() +
  labs(x= "Population", y = "Mean # of heterozygous exomic SNP/10kb") +
  geom_text(data = het_res, aes(label = paste0('T-test', ', p = ',p)), y = 12) +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
                                            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                            text = element_text(size = 18), axis.line = element_line(),
        legend.position = 'none') 

ggsave('Plots/hetdiv_exome_nuc_all_signif.png', width = 8, height = 6)



