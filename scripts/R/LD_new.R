# Linkage disequilibrium for Aim 1 - Sm Single GT project
# Author: Kathrin Jutzeler
# Date: May 5, 2023
# Last updated: November 30, 2023
# R version 4.2.0, tidyverse version 1.3.2

library(tidyverse)
library(fields)
library(pracma)
library(patchwork)
library(rstatix)

lab_colors = c("OR" = "#B22222", "BRE" = "#00A08A", "LE" = "#F98400", "EG" = "#5BBCD6", "NMRI" = '#8A2BE2')
field_colors <- c('Brazil' = "#9986A5", 'Senegal' = "#79402E", 'Niger'="#CCBA72" ,'Tanzania'='#0F0D0E')

all_colors <- c(lab_colors, field_colors)

order <- c('BRE', 'EG', 'LE', 'OR', 'NMRI', 'Brazil', 'Niger', 'Senegal', 'Tanzania')


#---------------#
# Functions ####
#---------------#

# Define bins to calculate the mean
#breaks <- seq(1, 1e5, length.out = 1000)

breakslog <- pracma::logseq(1, 1000000, n=1000)

# This function calculates the distance between variants and the R2 mean in each bin.

f_ld <- function(df) {
  
  # Get the distance between SNPs
  df <- df %>%
    group_by(CHR_A) %>%
    mutate(dist = BP_B-BP_A)
  
  df <- df %>%
     filter(dist <= 1000000)
  
  
  #bin r2 values and calculate stats
  ld_binned <- stats.bin(df$dist, df$R2, breaks = breakslog)
  
  # Create a new data frame 
  output <- data.frame(centers = ld_binned$centers, 
                       mean = ld_binned$stats["mean",])
  
  
  # Smooth the data and add a new column with predicted values
  
  loessMod  <- loess(output$mean ~ log(output$centers))
  output$smooth <- predict(loessMod, newdata = output)
  
  return(output)
  
}

# This function calculates the position at which LD decays to 0.5.

f_decay <- function(df){
  
  decay <- as.data.frame(approx(df$centers, df$smooth, n =120000))
  
  decay <- decay[which.min(abs(0.5-decay$y)),]
  return(decay)
  
}

#--------------------#
# Import the data ####
#--------------------#
# Part for whole genome lab and exon field data
setwd('D:/Documents/Anderson_lab/Diversity_paper/Analysis/post_review/ld/main')

#files <- list.files()
maf_files <- list.files()#grep(list.files(), pattern='ann', invert=TRUE, value=TRUE)

#files <- grep(list.files(path="."), pattern='maf', invert=TRUE, value=TRUE)

#LD_list <- lapply(files, read.table, header =T)
LD_list_maf <- lapply(maf_files[1:9], read.table, header =T)

#LD_list <- lapply(files, read.table, header =T)

# For LD adjusted (WG lab, CDS field)
names(LD_list_maf) <- c('Brazil', 'Niger', 'Senegal', 'Tanzania', 'BRE',  'EG', 'LE', 'NMRI',
                        'OR')

#-------------------------#
# Manipulate the data ####
#------------------------#

# Apply the function to each data set
LD_list_maf <- lapply(LD_list_maf, f_ld)
#LD_list_an <- lapply(LD_list, f_ld)

# Combine data sets
LD_df_maf <- bind_rows(LD_list_maf, .id = 'pop')
#LD_df_an <- bind_rows(LD_list_an, .id = 'pop')

# Order population 
LD_df_maf$pop <- factor(LD_df_maf$pop, levels = order)
#LD_df_an$pop <- factor(LD_df_an$pop, levels = order)

LD_df_maf <- LD_df_maf %>% 
  mutate(origin = ifelse(pop %in% c('BRE', 'EG', 'LE', 'OR', 'NMRI'), 'Laboratory', 'Field'))

# Get LD decay for each population 
decay_list_maf <- lapply(LD_list_maf, f_decay)
#decay_list_an <- lapply(LD_list_an, f_decay)

decay_df <- bind_rows(decay_list_maf, .id = 'pop')
#decay_df <- bind_rows(decay_list_an, .id = 'pop')

decay_df$pop <- factor(decay_df$pop, levels = order)

decay_df <- decay_df %>% mutate(origin = ifelse(pop %in% c('BRE', 'EG', 'LE', 'OR', 'NMRI'),'Laboratory', 'Field'))
decay_df$origin <- factor(decay_df$origin, levels = c('Laboratory', 'Field'))

#-----------#
# Plots  ####
#-----------#

LD_df_maf$origin  <- factor(LD_df_maf$origin, levels = c('Laboratory', 'Field'))


write.csv(LD_df_maf, "main.csv")

LD_df_maf <- read.csv('main.csv')
## Plot LD across 1e6 base pairs ####
p_ld_100 <- 
  ggplot(LD_df_maf, aes(centers, mean, color= pop ,linetype = origin)) +
  #geom_point(data = unld_df, aes(dist, R2),alpha =0.02, color = 'black') +
  geom_smooth(method = 'loess', se =F) +
  geom_point(alpha =0.1) +
  scale_y_continuous(name   = expression(bold(paste('Mean ', italic("r")^2))), 
                            expand = c(0,0), 
                            limits = c(0,1)) +
  scale_color_manual(values = all_colors) +
  xlab('Distance between SNPs (base pairs)') +
  scale_x_log10(limits = c(1e1, 1e6), breaks= c(1E1, 1e2, 1e3, 1e4, 1e5, 1e6), expand =c(0,0.1)) +
  labs(color = 'Population', linetype = 'Origin') +
  guides(color = guide_legend(nrow = 2, byrow = T), linetype = guide_legend(nrow = 2, 
                                                                            override.aes = list(color = 'black'))) +
  theme_minimal() +
  geom_segment(data = decay_df, aes(x = x, xend = x, y = 0, yend = 0.5), linetype = 'dashed') +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 12), axis.line = element_line(), legend.box.just = 'center', 
        legend.position = 'bottom' ) +
  annotation_logticks(sides = 'b') 

ggsave('Plots/LD_maf_log.png',width = 8, height = 6)
  

## Plot LD decay at R2 = 0.5 ####
decay_df %>%
  group_by(origin) %>%
  shapiro_test(x)

decay_res <- decay_df %>% wilcox_test(x ~ origin) %>% mutate(pop = 'EG', origin = 'Laboratory')

p_decay_100 <-
  ggplot(decay_df, aes(pop,  x, fill = pop)) +
  geom_col() +
  facet_grid(~factor(origin, levels = c('Laboratory', 'Field')), scales ='free') +
  scale_fill_manual(values = all_colors) +
  labs(y = expression(bold(paste("Position (bp) at ", italic("r")^2," = 0.5"))), x = "Population") +
  theme_minimal() +
  #scale_y_log10() + 
  #ylim(0, 16000) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 18), axis.line = element_line() ) +
  geom_text(data =decay_res, aes(label = paste0('Wilcoxon test', ', p = ', format(round(p,3),nsmall = 3))), y = 150000)


ggsave('Plots/LD_maf_decay.png',width = 8, height = 6)

#---------------#
# Exon data ####
#---------------#
#--------------------#
## Import the data ####
#--------------------#
#setwd('./LD_exon_100')

setwd('D:/Documents/Anderson_lab/Diversity_paper/Analysis/post_review/ld')
exon_files <- list.files(pattern = 'maf')#grep(list.files(), pattern='ann', invert=TRUE, value=TRUE)

#files <- grep(list.files(path="."), pattern='maf', invert=TRUE, value=TRUE)

#LD_list <- lapply(files, read.table, header =T)
LD_list_exon <- lapply(exon_files[1:9], read.table, header =T)

#LD_list <- lapply(files, read.table, header =T)

# For LD adjusted (WG lab, CDS field)
names(LD_list_exon) <- c('Brazil', 'BRE', 'EG', 'LE', 'Niger', 'NMRI', 'OR', 'Senegal', 'Tanzania')


#-------------------------#
## Manipulate the data ####
#------------------------#

# Apply the function to each data set
LD_list_exon <- lapply(LD_list_exon, f_ld)
#LD_list_an <- lapply(LD_list, f_ld)

# Combine data sets
LD_df_exon <- bind_rows(LD_list_exon, .id = 'pop')
#LD_df_an <- bind_rows(LD_list_an, .id = 'pop')

# Order population 
LD_df_exon$pop <- factor(LD_df_exon$pop, levels = order)
#LD_df_an$pop <- factor(LD_df_an$pop, levels = order)

LD_df_exon <- LD_df_exon %>% 
  mutate(origin = ifelse(pop %in% c('BRE', 'EG', 'LE', 'OR', 'NMRI'), 'Laboratory', 'Field'))

LD_df_exon$origin  <- factor(LD_df_exon$origin, levels = c('Laboratory', 'Field'))
# Get LD decay for each population 
decay_list_exon <- lapply(LD_list_exon, f_decay)
#decay_list_an <- lapply(LD_list_an, f_decay)

decay_df_exon <- bind_rows(decay_list_exon, .id = 'pop')
#decay_df <- bind_rows(decay_list_an, .id = 'pop')
decay_df_exon$pop <- factor(decay_df_exon$pop, levels = order)

decay_df_exon <- decay_df_exon %>% mutate(origin = ifelse(pop %in% c('BRE', 'EG', 'LE', 'OR', 'NMRI'),'Laboratory', 'Field'))
decay_df_exon$origin <- factor(decay_df_exon$origin, levels = c('Laboratory', 'Field'))

write.csv(LD_df_exon,'exon_df.csv')
write.csv(decay_df_exon,'exon_decay_df.csv')

decay_df_exon <- read.csv("exon_decay_df.csv")
# Just for BRE CDS
BRE <- read.table(exon_files[2], header = T)

BRE <- BRE %>%
  group_by(CHR_A) %>%
  mutate(dist = BP_B-BP_A)

BRE <- mutate(BRE, pop = 'BRE', origin = 'Laboratory')

LD_df_exon <- read.csv("exon_df.csv")
LD_df_exon$origin  <- factor(LD_df_exon$origin, levels = c('Laboratory', 'Field'))


# Plot exon LD ####

p_LD_exon_100 <-
ggplot(subset(LD_df_exon, pop != 'BRE'), aes(centers, mean, color= pop ,linetype = origin)) +
#ggplot(BRE, aes(dist, R2, color= pop ,linetype = origin)) +
  geom_point(alpha =0.1) +
 # geom_point(data = BRE, aes(dist, R2),alpha =0.02) +
  geom_smooth(method = 'loess', se =F) +
  #geom_smooth(data = BRE, aes(dist,R2), method = 'loess', se =F) +
  
  scale_y_continuous(name   = expression(bold(paste('Mean ', italic("r")^2))), 
                     expand = c(0,0), 
                     limits = c(0,1)) +
  scale_color_manual(values = all_colors, breaks = order) +
  xlab('Distance between SNPs (base pairs)') +
  scale_x_log10(limits = c(1e1, 1e6), breaks= c( 1E1, 1e2, 1e3, 1e4, 1e5, 1e6), expand =c(0,0.1)) +
  labs(color = 'Population', linetype = 'Origin') +
  guides(color = guide_legend(nrow = 2, byrow = T, order =1), linetype = guide_legend(nrow = 2, 
                                                                            override.aes = list(color = 'black'))) +
  theme_minimal() +
  #geom_segment(data = decay_df_exon, aes(x = x, xend = x, y = 0, yend = 0.5), linetype = 'dashed') +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 12), axis.line = element_line(), legend.box.just = 'center', 
        legend.position = 'bottom' ) +
  annotation_logticks(sides = 'b') 

p_LD_exon <- p_LD_exon_100 +
  geom_point(data = BRE, aes(dist, R2),alpha =0.02) +
  geom_smooth(data = BRE, aes(dist,R2), method = 'loess', se =F) 

decay_df_exon %>%
  group_by(origin) %>%
  shapiro_test(x)

decay_exon_res <- decay_df_exon %>% wilcox_test(x ~ origin) %>% mutate(pop = 'EG', origin = 'Laboratory')


p_decay_exon_100 <-
  ggplot(decay_df_exon, aes(pop,  x, fill = pop)) +
  geom_col() +
  facet_grid(~factor(origin, levels = c('Laboratory', 'Field')), scales ='free') +
  scale_fill_manual(values = all_colors) +
  labs(y = expression(bold(paste("Position (bp) at ", italic("r")^2," = 0.5"))), x = "Population") +
  theme_minimal() +
  #scale_y_log10() + 
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 18), axis.line = element_line() ) +
  geom_text(data =decay_exon_res, aes(label = paste0('Wilcoxon test', ', p = ', 
                                                     format(round(p,3),nsmall = 3))), y = 250000)

# Post review for unlinked SNPs ####
library(dplyr)
library(fuzzyjoin)

#setwd('D:/Documents/Anderson-lab/Diversity_paper/Analysis/post_review/ld/unlinked')

setwd('D:/Documents/temp')

unld <- list.files()

#EG_df <- read.table(unld[1], header =T)
#EG_pos <- read.table(unld[3], header =F)

#colnames(EG_pos) = c('chr', 'start', 'end') 

# Fuzzy join on BP_A
f_filter <- function(input_df, pos_file){
  ld_df <- read.table(input_df, header =T)
  pos_df <- read.table(pos_file, header =F)
  
  colnames(pos_df) = c('chr', 'start', 'end') 
  
EG_A <- fuzzy_left_join(
  ld_df,
  pos_df,
  by = c("BP_A" = "start", "BP_A" = "end"),
  match_fun = list(`>=`, `<=`)
) %>%
  mutate(CHR_A = ifelse(!is.na(chr), chr, "other")) %>%
  select(-chr, -start, -end)

# Fuzzy join on BP_B
EG_AB <- fuzzy_left_join(
  EG_A,
  pos_df,
  by = c("BP_B" = "start", "BP_B" = "end"),
  match_fun = list(`>=`, `<=`)
) %>%
  mutate(CHR_B = ifelse(!is.na(chr), chr, "other")) %>%
  select(-chr, -start, -end)

# Done
final_df <- EG_AB %>%
  filter(CHR_A != CHR_B)
return(final_df)

}

EG <- f_filter(unld[12], unld[5])
br <- f_filter(unld[1], unld[3])
BRE <- f_filter(unld[2], unld[4])
LE <- f_filter(unld[13], unld[6])
ng <- f_filter(unld[14], unld[7])
NMRI <- f_filter(unld[15], unld[8])
OR <- f_filter(unld[16], unld[9])
sn <- f_filter(unld[17], unld[10])
tz <- f_filter(unld[19], unld[11])

unld_list <- list(br, BRE, EG, LE, ng, NMRI, OR, sn, tz)

f_unld <- function(df) {
  df <- df %>%
  mutate(dist = BP_B-BP_A) %>%
  dplyr::select(dist, R2)
}

t <- lapply(unld_list, f_unld)

names(t) <- c('Brazil', 'BRE', 'EG', 'LE', 'Niger', 'NMRI',  'OR', 'Senegal', 'Tanzania')

unld_df <- bind_rows(t, .id = 'pop')

unld_df <- unld_df %>% mutate(origin = ifelse(pop %in% c('BRE', 'EG', 'LE', 'OR', 'NMRI'),'Laboratory', 'Field'))


unld_df %>%
  group_by(pop) %>%
  count()

unld_df <- read.csv('unlinked_snps.csv')

unld_df$pop <- factor(unld_df$pop, levels = order)
unld_df$origin <- factor(unld_df$origin, levels = c('Laboratory', 'Field'))

p_unld <- ggplot(unld_df, aes(pop, R2, fill = pop)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(~origin, scales ='free') +
  scale_fill_manual(values = all_colors, breaks = order) +
  labs(fill = 'Population', x = "Population") + theme_minimal() +
  #geom_segment(data = decay_df_exon, aes(x = x, xend = x, y = 0, yend = 0.5), linetype = 'dashed') +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 12), axis.line = element_line(), legend.box.just = 'center', 
        legend.position = 'bottom' )

#unld_df <- read.csv('unlinked_snps.csv')
# 
# unld_df <- unld_df %>% filter(dist > 1E6)
# 
# # Plot unlinked LD ####
# p_unld <- 
# ggplot(unld_df, aes(dist, R2, color= pop ,linetype = origin)) +
#   geom_smooth(method = 'loess', se =F) +
#   geom_point(alpha =0.1) +
#   scale_y_continuous(name   = expression(bold(paste('Mean ', italic("r")^2))), 
#                      expand = c(0,0), 
#                      limits = c(0,1)) +
#   scale_color_manual(values = all_colors, breaks = order) +
#   xlab('Distance between SNPs (base pairs)') +
#   scale_x_log10(limits = c(1e5, 1e7), breaks= c( 1E1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7), expand =c(0,0.1)) +
#   labs(color = 'Population', linetype = 'Origin') +
#   guides(color = guide_legend(nrow = 2, byrow = T, order =1), linetype = guide_legend(nrow = 2, 
#                                                                             override.aes = list(color = 'black'))) +
#   theme_minimal() +
#   #geom_segment(data = decay_df_exon, aes(x = x, xend = x, y = 0, yend = 0.5), linetype = 'dashed') +
#   theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         text = element_text(size = 12), axis.line = element_line(), legend.box.just = 'center', 
#         legend.position = 'bottom' ) +
#   annotation_logticks(sides = 'b') 
# 
# p_ld_100 + p_unld  

# Export plots ####
tiff("Figure4_LD_100.tiff", width = 8, height = 10, units = 'in', res = 600)

(p_ld_100 / p_decay_100)  +
  plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold', size = 18))

dev.off()

tiff("FigureS3_LD_exon_100.tiff", width = 8, height = 10, units = 'in', res = 600)

(p_LD_exon / p_decay_exon_100)  +
  plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold', size = 18))

dev.off()


tiff("UNL_LD.tiff", width = 8, height = 5, units = 'in', res = 600)
p_unld
dev.off()
# Get ranges to determine fold change ##

decay_range <- decay_df %>%
  #filter(pop != 'BRE') %>%
  group_by(origin) %>%
  reframe(range = range(x, na.rm = T))

lab_mat <- matrix(c(decay_range$range[1:2], decay_range$range[2:1]), 2,2)
field_mat <- matrix(c(decay_range$range[3:4], decay_range$range[3:4]), 2,2)

lab_mat / field_mat
