library(radiator)
library(SNPRelate)
library(tidyverse)
library(RLDNe)
library(lemon)

setwd('./vcf/estimator')

vcf.maf<- 'vcf/maf05_nuclear_CDS_nosex_nodup_v102.vcf'

vcf.in <- 'vcf/annotated_snps_nuclear_CDS_admixture_v10.vcf'

vcf.sub <- 'vcf/annotated_subset.vcf'

meta <- read_csv('metadata.csv')

#f <- list.files(pattern='anno')

setwd('./vcf/ne_colony')
vcf.list<- list.files(pattern=".vcf")

f_estimator <- function(vcf.in){
  
  name <- paste0(vcf.in, '55', '.gds')
  snpgdsVCF2GDS(vcf.in, out.fn = name, method="biallelic.only")
  
  gp_out <- genomic_converter(data = name, output = "genepop")
  #gp_tidy <- write_genepopedit(gp_out$tidy.data) # Converts 3 letter format to single letter.
  
  return(gp_out)
  
}

# For some reason can't use radiator in loop :(
#o2 <- lapply(vcf.list[1:9], f_estimator)

br_gp <- lapply(vcf.list[1], f_estimator)
BRE_gp <- lapply(vcf.list[2], f_estimator)
EG_gp <- lapply(vcf.list[3], f_estimator)
LE_gp <- lapply(vcf.list[4], f_estimator)
ne_gp <- lapply(vcf.list[5], f_estimator)
NMRI_gp <- lapply(vcf.list[6], f_estimator)
OR_gp <- lapply(vcf.list[7], f_estimator)
sn_gp <- lapply(vcf.list[11], f_estimator)
tz_gp <- lapply(vcf.list[9], f_estimator)

gp_list <- c(br_gp, BRE_gp, EG_gp, LE_gp, ne_gp, NMRI_gp, OR_gp, sn_gp, tz_gp)
gp_list <- unlist(gp_list, recursive = F)

f_tidy <- function(x){
  
  gp_out <- x
  gp_out$MARKERS <- gsub('SM_V10_[1-9]__','', gp_out$MARKERS)
  gp_out$MARKERS <- gsub('__.*', '', gp_out$MARKERS)
  
  return(gp_out)

}

gp_out <- lapply(gp_list, f_tidy)

names(gp_out) <- c('Brazil', 'BRE', 'EG', 'LE', 'Niger', 'NMRI', 'OR', 'Senegal', 'Tanzania')

# Write individual files for each population
lapply(seq(gp_out), function(x) write_genepop(gp_out[[x]], filename = names(gp_out)[[x]])) #genepop.header = names(x)))

# Make a list with chromosomes and positions
gp_df <- bind_rows(gp_out)
chrom <- gp_df$CHROM
locus <- gp_df$LOCUS

s <- unique(data.frame(chrom, locus))

write_delim(x = s, file = 'chr_pos.txt')


#gds <- snpgdsVCF2GDS(vcf.in, "ne_anno2.gds", method="biallelic.only")
#gds <- snpgdsVCF2GDS(vcf.maf, "ne_maf5.gds", method="biallelic.only")
#gds <- snpgdsVCF2GDS(vcf.sub, "ne_sub.gds", method="biallelic.only")

#gp_out <- genomic_converter(data = 'ne.gds', output = "genepop") # Puts out a file in a table but in 3 letter code
#gp_out <- genomic_converter(data = 'ne_anno2.gds', output = "genepop")
#gp_tidy <- write_genepopedit(gp_out$tidy.data) # Converts 3 letter format to single letter. 


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

#Ne_estimates<-readLDNe_tab(path = param_files$Ne_out_tab)

#----------#
#IMPORT FILE
#----------#

f <- list.files()

all_lines <- lapply(f, function(x) readLines(x))

#start <- lapply(all_lines, function(x) as.list(which(grepl("Population", x))))

f_extract <- function(x){
  start <- which(grepl("Population", x))
  x[start:(start+21)]
}
#start <- as.list(which(grepl("Population", all_lines_maf)))

df_list <- lapply(all_lines, f_extract)
#df_list <- lapply(start, function(x) all_lines[x:(x+21)])

f_clean <- function(x){
  
    x <- gsub("\\s\\s+", "\t",x )
    data.frame(pop = strsplit(x[1],split = "\t")[[1]][2],
             Estimated_Ne_0.1 = strsplit(x[13],split = "\t")[[1]][2],
             Estimated_Ne_0.05 = strsplit(x[13],split = "\t")[[1]][3],
             #Estimated_Ne_0.01 = strsplit(x[13],split = "\t")[[1]][4],
             "No S*" = strsplit(x[13],split = "\t")[[1]][5],
             "Estiamted_Ne_0+" = strsplit(x[13], split = "\t")[[1]][5],
             "Expected r^2_0.1" = strsplit(x[11], split = "\t")[[1]][2],
             "Expected r^2_0.05" = strsplit(x[11], split = "\t")[[1]][3],
             "Expected r^2_0.02" = strsplit(x[11], split = "\t")[[1]][4],
             "Expected r^2_0+" = strsplit(x[11], split = "\t")[[1]][5])
}

# to test
#t <- gsub("\\s\\s+", '\t', sn_list[[1]])
#strsplit(t,split = "\t")[[13]] [5]
  
df_clean <- lapply(df_list,f_clean)

df_clean_maf <- lapply(df_list,f_clean)

df_ne <- bind_rows(df_clean)
df_ne <- bind_rows(df_clean_maf)

#df_ne <- df_ne[-8:-9,]
df_ne$pop <- c('Brazil', 'BRE', 'EG', 'LE', 'Niger', 'NMRI', 'OR', 'Senegal', 'Tanzania')
#df_ne$pop <- c('Brazil', 'BRE', 'EG', 'LE', 'Niger', 'NMRI', 'OR', 'Outgroup', 'Tanzania', 'Senegal')

df_ne2 <- df_ne %>%
  dplyr::select(1:5) %>%
    gather('ne', 'value', 2:5)

df_ld <- df_ne %>%
  dplyr::select(1,6:9) %>%
  gather('ld', 'value', 2:5)

df_ne2$value <- gsub('Infinite', Inf, df_ne2$value)
       
df_ne2$value <- as.numeric(df_ne2$value)
df_ld$value <- as.numeric(df_ld$value)
# PLOT

lab_colors = c("OR" = "#B22222", "BRE" = "#00A08A", "LE" = "#F98400", "EG" = "#5BBCD6", "NMRI" = '#8A2BE2')
field_colors <- c('Brazil' = "#9986A5", 'Senegal' = "#79402E", 'Niger'="#CCBA72" ,'Tanzania'='#0F0D0E')

all_colors <- c(lab_colors, field_colors)

order <- c('BRE', 'EG', 'LE', 'OR', 'NMRI', 'Brazil', 'Niger', 'Senegal', 'Tanzania')
df_ne2$pop <- factor(df_ne2$pop, levels = order)

ggplot(df_ne2, aes(x = ne, y = value, fill = pop)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = all_colors) +
  theme_minimal() +
  #scale_x_discrete(labels=c('0.02', '0.05', '0.1', '0+')) +
  labs(x = 'Crtitical Value', y = 'Estimated population size', fill = 'Population') +
  theme(panel.border=element_blank(), axis.line = element_line()) +
  #scale_y_log10(breaks = c(1, 10, 100, 1000, 10000), labels = c(1, 10, 100, 1000, 'Infinite')) +
  coord_capped_cart( left = 'both')

ggsave("Plots/Ne_010223.png", width = 8, height =6)
  

df_ld$pop <- factor(df_ld$pop, levels = order)

ggplot(df_ld, aes(x = ld, y = value, fill = pop)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = all_colors) +
  theme_minimal() +
  scale_x_discrete(labels=c('0.02', '0.05', '0.1', '0+')) +
  labs(x = 'Allele Frequency used', y = 'Estimated R^2', fill = 'Population') +
  theme(panel.border=element_blank(), axis.line = element_line()) +
  #scale_y_log10(breaks = c(1, 10, 100, 1000, 10000), labels = c(1, 10, 100, 1000, 'Infinite')) +
  coord_capped_cart(ylim= c(0, 0.25), left = 'both')

ggsave("Plots/LD_maf.png", width = 8, height =6)



## other file ####

lines <- readLines('Ne_outxLD_all_anno.txt')
lines <- readLines('Ne_out_redoxLD_partial.txt')

start <- which(lines == paste0(rep("-",25), collapse = ""))+1
end <-  which(lines == paste0(rep("-",37), collapse = ""))-1

df <- lines[start:end]

#samples <- df[which(grepl('[0-9]:', df))]

#samples <- gsub('\\s\\s+', '\t', samples)

df2 <- gsub('\\s\\s+', '\t', df)

start2 <- which(grepl('[0-9]:', df2))

subdf <- lapply(start2, function(x) df2[x:(x+3)])

t <- subdf[[1]]

f_cleaner <- function(x){
  data.frame(pop = strsplit(x, split = "\t")[[1]][2],
             sample_no = strsplit(x, split = "\t")[[1]][3],
             crit_value = c(strsplit(x, split = "\t")[[1]][4], strsplit(x, split = "\t")[[2]][2], 
                            strsplit(x, split = "\t")[[3]][2], strsplit(x, split = "\t")[[4]][2]),
             NE = c(strsplit(x, split = "\t")[[1]][9], strsplit(x, split = "\t")[[2]][7], 
                    strsplit(x, split = "\t")[[3]][7], strsplit(x, split = "\t")[[4]][7]),
             CI_low = c(strsplit(x, split = '\t')[[1]][10], strsplit(x, split = "\t")[[2]][8],
                        strsplit(x, split = "\t")[[3]][8], strsplit(x, split = "\t")[[4]][8]),
             CI_high = c(strsplit(x, split = '\t')[[1]][11], strsplit(x, split = "\t")[[2]][9],
                         strsplit(x, split = "\t")[[3]][9], strsplit(x, split = "\t")[[4]][9]))  
}

fin_df <- lapply(subdf, f_cleaner)
fini_df <- bind_rows(fin_df)

fini_df$NE<- as.numeric(fini_df$NE)
fini_df$CI_high <- as.numeric(fini_df$CI_high)
fini_df$CI_low <- as.numeric(fini_df$CI_low)

fini_df$pop <- factor(fini_df$pop)

fini_df$pop <- recode_factor(fini_df$pop, "5:PdV_0447_1" = 'Brazil', "4:7_S32_L003" = 'BRE', "3:0_S40_L003" = 'EG', 
              "2:4_S18_L003" = 'LE', "6:NE_Di158_1" = 'Niger',  "7:_nmri_1_S7" = 'NMRI', "1:1_S2_L003" = 'OR',
              "8:TZ_009_1_1" = 'Tanzania', "9:sm_senegal" = 'Senegal')

fini_df$pop <- factor(fini_df$pop, levels = order)

# PLOT

ggplot(fini_df, aes(x = crit_value, y = NE, fill = pop)) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), position = position_dodge()) + 
  scale_fill_manual(values = all_colors) +
  theme_minimal() +
  #scale_x_discrete(labels=c('0.02', '0.05', '0.1', '0+')) +
  labs(x = 'Allele Frequency used', y = 'Estimated population size', fill = 'Population') +
  theme(panel.border=element_blank(), axis.line = element_line()) +
  #coord_cartesian(ylim=c(-750,1500)) +
  coord_capped_cart(ylim=c(-750,1500), left = 'both')

ggsave('Plots/NE_anno2.png', width = 8, height = 6)


ggplot(fini_df, aes(x = pop, y = NE, fill = pop)) +
  geom_boxplot() +
  scale_fill_manual(values = all_colors) +
  theme_minimal() +
  #scale_x_discrete(labels=c('0.02', '0.05', '0.1', '0+')) +
  labs(x = 'Population', y = 'Estimated population size', fill = 'Population') +
  theme(panel.border=element_blank(), axis.line = element_line(), legend.position = 'none') +
  #coord_cartesian(ylim=c(-750,1500)) +
  coord_capped_cart(ylim=c(-750,1500), left = 'both')

ggsave('Plots/NE_boxplot.png', width = 8, height = 6)


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

ne <- lapply(ne.list, function(x) readLines(x)[5])

ne  <- gsub('\\s+', '\t', ne)

ne <- strsplit(ne, split = "\t")

names(ne) <- c('Brazil','BRE', 'EG', 'LE', 'Niger', 'NMRI', 'OR', 'Senegal', 'Tanzania')

o <- bind_rows(ne, .id = 'pop')

o <- o[-1:-2,]

o <- gather(o, 'pop', 'ne', 1:ncol(o))
o$ne <- as.numeric(o$ne)

o$pop <- factor(o$pop, levels = order)

ggplot(o, aes(x =pop, y = ne, fill = pop)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = all_colors) +
  theme_minimal() +
  #coord_cartesian(ylim = c(0, 1000)) +
  #scale_y_continuous(breaks = c(0, 100, 250,500,3000, 5000), labels=c('0', '100', '250', '500', '3000', 'Infinite')) +
  labs(x = 'Population', y = 'Estimated population size', fill = 'Population') +
  scale_y_log10(breaks = c(0, 100, 1e5, 1e9), labels = c('0', '100', '1e+05', 'Infinite')) +
  theme(panel.border=element_blank(), axis.line = element_line(), legend.position = 'none', text =
          element_text(size = 12)) +
  coord_capped_cart(left = 'top')

ggsave('Plots/ne_sibship_more_variants.png', width = 8, height =6)
