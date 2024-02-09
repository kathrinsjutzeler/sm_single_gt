# Linkage disequilibrium for Aim 1 - Sm Single GT project
# Author: Kathrin Jutzeler
# Date: May 5, 2023
# Last updated: February 8 2023
# R version 4.2.0, tidyverse version 1.3.2

library(tidyverse)
library(fields)
library(pracma)

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
  
  decay <- as.data.frame(approx(df$centers, df$smooth, n =100000))
  
  decay <- decay[which.min(abs(0.5-decay$y)),]
  return(decay)
  
}

#--------------------#
# Import the data ####
#--------------------#
setwd('./ld_v10')

#files <- list.files()
maf_files <- list.files( pattern = 'maf')

LD_list_maf <- lapply(maf_files, read.table, header =T)

#LD_list <- lapply(files, read.table, header =T)

names(LD_list_maf) <- c('Brazil', 'BRE', 'EG', 'LE', 'Niger', 'NMRI', 'OR', 'Senegal', 'Tanzania') 

#-------------------------#
# Manipulate the data ####
#------------------------#

# Apply the function to each data set
LD_list_maf <- lapply(LD_list_maf, f_ld)

# Combine data sets
LD_df_maf <- bind_rows(LD_list_maf, .id = 'pop')

# Order population 
LD_df_maf$pop <- factor(LD_df_maf$pop, levels = order)

# Assign origin
LD_df_maf <- LD_df_maf %>% 
  mutate(origin = ifelse(pop %in% c('BRE', 'EG', 'LE', 'OR', 'NMRI'), 'Lab', 'Field'))

# Get LD decay for each population 
decay_list_maf <- lapply(LD_list_maf, f_decay)

decay_df <- bind_rows(decay_list_maf, .id = 'pop')

decay_df$pop <- factor(decay_df$pop, levels = order)

decay_df <- decay_df %>% mutate(origin = ifelse(pop %in% c('BRE', 'EG', 'LE', 'OR', 'NMRI'),'Lab', 'Field'))
decay_df$origin <- factor(decay_df$origin, levels = c('Lab', 'Field'))

#-----------#
# Plots  ####
#-----------#

LD_df_maf$origin  <- factor(LD_df_maf$origin, levels = c('Lab', 'Field'))

## Plot LD across 1e6 base pairs ####
p_ld <- 
  ggplot(LD_df_maf, aes(centers, mean, color= pop, linetype = origin)) +
  geom_smooth(method = 'loess', se =F) +
  geom_point(alpha =0.1) +
  scale_y_continuous(name   = expression(bold(paste('Mean ', italic("r")^2))), 
                            expand = c(0,0), 
                            limits = c(0,1)) +
  scale_color_manual(values = all_colors) +
  xlab('Distance between SNPs (base pairs)') +
  scale_x_log10(limits = c(1e1, 1e6), breaks= c( 1E1, 1e2, 1e3, 1e4, 1e5, 1e6)) +
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

decay_res <- decay_df %>% t_test(x ~ origin) %>% mutate(pop = 'EG', origin = 'Lab')


p_decay <-
  ggplot(decay_df, aes(pop,  x, fill = pop)) +
  geom_col() +
  facet_grid(~factor(origin, levels = c('Lab', 'Field')), scales ='free') +
  scale_fill_manual(values = all_colors) +
  labs(y = expression(bold(paste("Position (bp) at ", italic("r")^2," = 0.5"))), x = "Population") +
  theme_minimal() +
  #scale_y_log10() + 
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 18), axis.line = element_line() ) +
  geom_text(data =decay_res, aes(label = paste0('T-test', ', p = ', format(round(p,3)),1)), y = 2.5e+05)


ggsave('Plots/LD_maf_decay.png',width = 8, height = 6)

# Export plots ####
pdf("Figure5_LD.pdf", width = 8, height = 10)

(p_ld / p_decay)  +
  plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold'))

dev.off()
