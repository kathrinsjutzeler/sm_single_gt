# Tajima's D
# Author: Kathrin Jutzeler
# Date: February 16, 2023
# Last updated: February 8, 2024
# R version 4.2.0, tidyverse version 1.3.2

library(ggridges)
library(tidyverse)
library(stats)
library(lemon)

# Define colors
lab_colors = c("OR" = "#B22222", "BRE" = "#00A08A", "LE" = "#F98400", "EG" = "#5BBCD6", "NMRI" = '#8A2BE2')
field_colors <- c('Brazil' = "#9986A5", 'Senegal' = "#79402E", 'Niger'="#CCBA72" ,'Tanzania'='#0F0D0E')

all_colors <- c(lab_colors, field_colors)

order <- c('BRE', 'EG', 'LE', 'OR', 'NMRI', 'Brazil', 'Niger', 'Senegal', 'Tanzania')

setwd('./tajima_v10')


# Import the data ---------------------------------------------------------------

files <- list.files(pattern = 'ann') # This is for annotated_snp file

tajima_list <- lapply(files, read.delim, header =T)

names(tajima_list) <- c('Brazil', 'BRE', 'EG', 'LE', 'Niger', 'NMRI', 'OR', 'Senegal', 'Tanzania') 

tajima_df <- bind_rows(tajima_list, .id = 'pop')

tajima_df$pop <- factor(tajima_df$pop, levels = order)

tajima_df <- tajima_df %>% mutate(origin = ifelse(pop %in% c('BRE', 'EG', 'LE', 'OR', 'NMRI'),'Lab', 'Field'))

tajima_df$origin <- factor(tajima_df$origin, levels = c('Lab', 'Field'))

tajima_df %>%
  group_by(pop) %>%
  na.omit() %>%
  count()

wilcox.test(TajimaD ~ origin, data = tajima_df)

# Plot boxplot ----------------------------------------------------------------------
ggplot(tajima_df, aes(pop,y = TajimaD, fill=pop )) +
  geom_boxplot() + 
  facet_grid(~origin, scales ='free') +
  #scale_y_log10() +
  scale_fill_manual(values = all_colors)+
  theme_minimal() +
  labs(y="Tajima's D", x='Population') +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 18), axis.line = element_line(),
        legend.position = 'none') +
  
  coord_capped_cart(bottom = 'both')

ggsave('Plots/tajima_boxplot_v10.png', width = 8, height = 6)


# Plot barplot -------------------------------------------------------------------
avg_tajima <- tajima_df %>%
  group_by(pop, origin) %>%
  summarise(mean = mean(TajimaD, na.rm =T), sd = sd(TajimaD, na.rm = T), se = plotrix::std.error(TajimaD, na.rm = T))

avg_tajima %>% 
  group_by(origin) %>%
  shapiro_test(mean)

tajima_res <- avg_tajima %>% ungroup() %>% t_test(mean ~ origin) %>% mutate(pop = 'EG', origin = 'Lab')

p_tajima <- ggplot(avg_tajima, aes(pop,y = mean, fill=pop )) +
  geom_errorbar(data = subset(avg_tajima, origin == "Lab"), aes(ymin=mean-0.01, ymax=mean+sd), width = 0.5) +
  geom_errorbar(data = subset(avg_tajima, origin == "Field"), aes(ymin=mean+0.01, ymax=mean-sd), width = 0.5) +
  geom_col() + 
  scale_fill_manual(values = all_colors) +
  facet_grid(~factor(origin, levels = c('Lab', 'Field')), scales ='free') +
  theme_minimal() +
  labs(y="Tajima's D", x='Population') +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 18), axis.line = element_line(), legend.position = 'none') +
  geom_text(data = tajima_res, aes(label = paste0('T-test', ', p = ',round(p,3)), y = 2.2))

ggsave('Plots/tajima_barplot_v10_se.png', width = 8, height = 6)

setwd('..')

# Generate figure --------------------------------------------------------------

pdf("Figure3.pdf", width = 10, height = 6)

(p_tajima) # +
  #plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold'))

dev.off()
