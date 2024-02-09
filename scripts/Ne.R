# Effective population size for Aim 1 - Sm Single GT project
# Author: Kathrin Jutzeler
# Date: May 5, 2023
# Last updated: February 8, 2024
# R version 4.2.0, tidyverse version 1.3.2

library(radiator)
library(SNPRelate)
library(tidyverse)
library(RLDNe)
library(lemon)
library(ggbreak)
library(ggpubr)

setwd('./vcf/estimator')

# Read metadata
meta <- read_csv('metadata.csv')

#-----------------------------------------------------------------------#
# Import vcf files to generate a suitable input file for NeEstimator ####
#-----------------------------------------------------------------------#

setwd('./vcf/ne_colony')
vcf.list<- list.files(pattern=".vcf")

f_estimator <- function(vcf.in){
  
  name <- paste0(vcf.in, '425', '.gds')
  snpgdsVCF2GDS(vcf.in, out.fn = name, method="biallelic.only")
  
  gp_out <- genomic_converter(data = name, output = "genepop")
  #gp_tidy <- write_genepopedit(gp_out$tidy.data) # Converts 3 letter format to single letter.
  
  return(gp_out)
  
}

# For some reason can't use radiator in loop :(
#o2 <- lapply(vcf.list[1:9], f_estimator)

br_gp <- lapply(vcf.list[1], f_estimator)
BRE_gp <- f_estimator("annotated_snps_nuclear_CDS_nosex_v10_BRE.vcf")
EG_gp <- lapply(vcf.list[3], f_estimator)
LE_gp <- f_estimator("annotated_snps_nuclear_CDS_nosex_v10_LE.vcf")
ne_gp <- lapply(vcf.list[5], f_estimator)
NMRI_gp <- f_estimator("annotated_snps_nuclear_CDS_nosex_v10_NMRI.vcf")
OR_gp <- lapply(vcf.list[7], f_estimator)
sn_gp <- f_estimator('annotated_snps_nuclear_CDS_nosex_v10_sn.vcf')
tz_gp <- f_estimator('annotated_snps_nuclear_CDS_nosex_v10_tz.vcf')

gp_list <- c(br_gp, EG_gp, ne_gp,  OR_gp)
gp_list <- unlist(gp_list, recursive = F)

gp_list <- c(gp_list, BRE_gp, LE_gp, NMRI_gp, sn_gp, tz_gp)

f_tidy <- function(x){
  
  gp_out <- x
  gp_out$MARKERS <- gsub('SM_V10_[1-9]__','', gp_out$MARKERS)
  gp_out$MARKERS <- gsub('__.*', '', gp_out$MARKERS)
  
  return(gp_out)

}

gp_out <- lapply(gp_list, f_tidy)

names(gp_out) <- c('Brazil', 'EG', 'Niger', 'OR', 'BRE','LE', 'NMRI', 'Senegal', 'Tanzania')

# Write individual files for each population
lapply(seq(gp_out), function(x) write_genepop(gp_out[[x]], filename = names(gp_out)[[x]])) #genepop.header = names(x)))

# Make a list with chromosomes and positions
gp_df <- bind_rows(gp_out)
chrom <- gp_df$CHROM
locus <- gp_df$LOCUS

s <- unique(data.frame(chrom, locus))

write_delim(x = s, file = 'chr_pos.txt')

#--------------------------------------------------------#
# To run through R but probably easier to run through GUI ####
#--------------------------------------------------------#

# Function to clean up date frame and make it work with NeEstimator (can be run for entire dataframe when 
# gp_list is merged with bind_rows)


f_cleanup <- function(x) {
  gp_tidy <- x
  gp_tidy <- gp_tidy[-2:-3]
gp_tidy$SampleID <- gsub('-', '_', gp_tidy$SampleID)

gp_tidy <- meta %>%
  dplyr::select(sample_id, pop) %>%
  right_join(gp_tidy, by = c('sample_id'='SampleID'))

#gp_tidy %>%
#  filter(sample_id == 'Sm_SN_Nd103_1', colnames(gp_tidy) == 'SM_V10_1__SM_V10_1_620966')

gp_tidy$sample_id <- gsub('.hc.vcf.gz', '', gp_tidy$sample_id)

# Apply the function to each column
df <- lapply(gp_tidy[3:ncol(gp_tidy)], remove_zeros)

fin_df<- bind_cols(gp_tidy[1:2],df)
}

gp_file<-write_genepop_zlr(loci = fin_df[,3:ncol(fin_df)],pops = fin_df$pop,
                           ind.ids = fin_df$sample_id,folder = "",filename ="genepop_output.txt",missingVal = NA,
                           ncode = 2,diploid = T)

param_files<- NeV2_LDNe_create(input_file = gp_file ,param_file = "Ne_params.txt", 
                               matingsystem = 0, NE_out_file = "Ne_out.txt")

run_LDNe(LDNe_params = param_files$param_file)


#----------------------------------#
#Import file NeEstimator files ####
#----------------------------------#
setwd('annot2')
f <- list.files()

all_lines <- lapply(f, function(x) readLines(x))


f_extract <- function(x){
  start <- which(grepl("Population", x))
  x[start:(start+21)]
}

df_list <- lapply(all_lines, f_extract)

names(df_list) <- c('Brazil', 'BRE', 'EG', 'LE', 'Niger', 'NMRI', 'OR', 'Senegal', 'Tanzania')

df_list <- lapply(df_list, function(x) as.data.frame(x[c(13,16:17)]))

neest_df <- lapply(df_list, function(x) separate(x, into = c('label', '0.1', '0.05', '0.02', 'No S*', '0+'), 
                                     sep = '\\s\\s+', col =`x[c(13, 16:17)]`))

neest_df <- bind_rows(neest_df, .id = 'pop')

df_ne2 <- 
  neest_df %>%
  gather('cv', 'value', 3:7)

df_ne2$value <- gsub('Infinite', Inf, df_ne2$value)
       
df_ne2$value <- as.numeric(df_ne2$value)

df_cv <- df_ne2 %>%
            filter(cv == '0+')

df_out <- df_cv %>% spread(label, value)

write.csv(df_out, 'ne_neest.csv')

# Plot

lab_colors = c("OR" = "#B22222", "BRE" = "#00A08A", "LE" = "#F98400", "EG" = "#5BBCD6", "NMRI" = '#8A2BE2')
field_colors <- c('Brazil' = "#9986A5", 'Senegal' = "#79402E", 'Niger'="#CCBA72" ,'Tanzania'='#0F0D0E')

all_colors <- c(lab_colors, field_colors)

order <- c('BRE', 'EG', 'LE', 'OR', 'NMRI', 'Brazil', 'Niger', 'Senegal', 'Tanzania')

df_ne2$pop <- factor(df_ne2$pop, levels = order)

df_ne3 <- df_ne2 %>%
    filter(cv == "0+") %>%
  mutate(origin = ifelse(pop %in% c('BRE', 'EG', 'LE', 'OR', 'NMRI'), 'Lab', 'Field'))

ne_breaks <- c(0, 100, 200, 300, 400, 500, 3000, 3500, 4000, 4500)

# Plot final NeEstimator ####
p_ne_neest <- ggplot(subset(df_ne3, label == 'Estimated Ne^ ='), aes(x = pop, y = value, fill = pop)) +
  geom_errorbar(data = subset(df_ne3, pop %in% c('BRE', 'EG', 'LE', 'OR', 'NMRI', 'Brazil') & 
                                label == ""), 
                aes(ymin=2, ymax=value, width = 0.5)) + 
    geom_col(position = position_dodge()) +
  scale_fill_manual(values = all_colors) +
  facet_grid(~factor(origin, levels = c('Lab', 'Field')), scales = 'free_x') +
  ylab(expression(bold('Effective population size ' ~ italic((N[e]))))) +
  scale_y_continuous(limits = c(0,4500), breaks = ne_breaks, 
                     labels = c('0', '100', '200', '300', '400', '500', '3000', '3500', '4000', 'Infinite')) +
  scale_y_break(c(510,3000), scales= 'free', space = 0.2) +
  #scale_y_break(c(500,3000), scales= 'free', space = 0.2) +
  labs(x = 'Population', fill = 'Population', title = 'NeEstimator') +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, vjust = 0.7), axis.title = element_text(size = 14, face = 'bold'),
        plot.title = element_text(hjust = 0.5, size = 16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 18), axis.line = element_line(), legend.position = 'none') +
  theme(axis.text.y.right = element_blank(), axis.line.y.right = element_blank(),axis.ticks.y.right = element_blank())
  
  
ggsave("Plots/Ne_annot_010324_all_SNPs.png", width = 8, height =6)


#-----------------------------------#
# Prepare input file for COLONY  ####
#-----------------------------------#

setwd('./vcf/ne_colony')
vcf.list<- list.files(pattern=".vcf")

meta <- read_csv('metadata.csv')

f_colony <- function(vcf.maf){
  
  gds <- snpgdsVCF2GDS(vcf.maf, paste0(vcf.maf, "ne2.gds"), method="biallelic.only")

#gp_out <- genomic_converter(data = 'ne.gds', output = "genepop") # Puts out a file in a table but in 3 letter code
  colony_out <- genomic_converter(data = paste0(vcf.maf, "ne2.gds"), output = "colony")
#gp_tidy <- write_genepopedit(gp_out$tidy.data) # Converts 3 letter format to single letter. 

  colony_out$tidy.data$INDIVIDUALS <- gsub('-', '_', colony_out$tidy.data$INDIVIDUALS)

  colony_out$tidy.data <- meta %>%
    dplyr::select(sample_id, pop) %>%
    right_join(colony_out$tidy.data, by = c('sample_id'='INDIVIDUALS'))


  colony_out$tidy.data$INDIVIDUALS <- colony_out$tidy.data$sample_id 

  colony_out$tidy.data$STRATA <- colony_out$tidy.data$pop

  colony_out$tidy.data <- colony_out$tidy.data[-1:-2]

  #gp_tidy$sample_id <- gsub('.hc.vcf.gz', '', gp_tidy$sample_id)
  
  colony_out$tidy.data$MARKERS <- gsub('SM_V10_[1-7]__', '', colony_out$tidy.data$MARKERS)

  tidy_colony <- write_colony(colony_out$tidy.data, filename = paste0(colony_out$tidy.data$STRATA[1], '.txt'))
}

lapply(vcf.list, f_colony)

# Run COLONY on ranch 

#-------------------------------#
# Import COLONY output files ####
#-------------------------------#
setwd('./vcf/ne_colony')

ne.list <- list.files(pattern=".Ne")

ne_list <- lapply(ne.list, function(x) as.data.frame(readLines(x)[5:7]))

names(ne_list) <- c('Brazil','BRE', 'EG', 'LE', 'Niger', 'NMRI', 'OR', 'Senegal', 'Tanzania')

ne_list <- lapply(ne_list, function(x) x %>% separate(`readLines(x)[5:7]`, sep = '\\s+', into = c('label', 'equal', 'value')))

ne_df <-  bind_rows(ne_list, .id = 'pop')

ne_df <- ne_df[-3]

ne_df <- ne_df %>%
  spread(label, value)

ne_df[2:4] <- sapply(ne_df[2:4], as.numeric)

write.csv(ne_df, 'ne_df.csv')

ne_df$Ne  <- ifelse(ne_df$Ne == '2147483647', Inf, ne_df$Ne)

ne_df$pop <- factor(ne_df$pop, levels = order)

ne_df <- ne_df %>% mutate(origin = ifelse(pop %in% c('BRE', 'EG', 'LE', 'OR', 'NMRI'), 'Lab', 'Field'))


# Plot colony Ne ####
p_ne_sib <- ggplot(ne_df, aes(x =pop, y = Ne, fill = pop)) +
  geom_errorbar(data = subset(ne_df, origin == 'Lab'), 
                aes(ymin=Ne-2, ymax=`CI95(U)`, width = 0.5)) + 
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = all_colors) +
  facet_grid(~factor(origin, levels = c('Lab', 'Field')), scales = 'free') +
  ylab(expression(bold('Effective population size ' ~ italic((N[e]))))) +
  labs(x = 'Population',  fill = 'Population', title = 'COLONY') +
  scale_y_continuous(limits = c(0,4500), breaks = ne_breaks, 
                     labels = c('0', '100', '200', '300', '400', '500', '3000', '3500', '4000', 'Infinite')) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, vjust = 0.7), axis.title = element_text(size = 14, face = 'bold'),
        plot.title = element_text(hjust = 0.5, size = 16),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 18), axis.line = element_line(), legend.position = 'none') +
  scale_y_break(c(510,3000), scales="free", space = 0.2) +
  theme(axis.text.y.right = element_blank(), axis.line.y.right = element_blank(),axis.ticks.y.right = element_blank())

ggsave('Plots/ne_sibship_more_variants.png', width = 8, height =6)

# Combine census, colony and neestimator output ####

df_ne_combo <-df_ne3 %>% filter(label == 'Estimated Ne^ =')

table <- ne_df %>%
  full_join(hmean_df, by = c('pop' = 'Sm strain' )) %>%
  full_join(df_ne_combo, by = 'pop')

table <- table[c(-2,-3,-6, -7,-9,-10, -12)]
colnames(table) <- c('pop', 'COLONY', 'origin', 'Census', 'NeEstimator')

table <- table %>% filter(origin == 'Lab', pop != 'NMRI' )
table <- table %>% gather('tool', 'ne', c(2,5))

# Plot correlation ####
p_ne_cor <- ggplot(table, aes(Census, ne, linetype =tool)) +
  geom_point(aes(color = pop), size =2, position = position_dodge(width = 1.5)) +
  scale_color_manual(values = all_colors) +
  geom_line() +
  stat_cor() +
  theme_minimal() +
  ylab(expression(bold('Effective population size ' ~ italic((N[e]))))) +
  xlab('Census population size') +
  labs(linetype = 'Tool', color = 'Population') +
  theme(axis.text.x = element_text(size = 12, vjust = 0.7), axis.title = element_text(size = 14, face = 'bold'),
        plot.title = element_text(hjust = 0.5, size = 16),
        text = element_text(size = 12), axis.line = element_line())
  

# Export plots ####

pdf('Figure6_Ne.pdf', width = 8, height = 10)

(p_ne_neest & theme(axis.title.x = element_blank()))  / p_ne_sib +
  plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold'))

dev.off()
