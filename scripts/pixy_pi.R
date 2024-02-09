# Pi statistics for Sm single gentoype project
# Author: Kathrin Bailey
# Date: October 13, 2022
# Last updated: February 8, 2024
# R version 4.2.0, tidyverse version 1.3.2

library(ggridges)
library(tidyverse)
library(stats)
library(rstatix)


# Define colors
lab_colors = c("OR" = "#B22222", "BRE" = "#00A08A", "LE" = "#F98400", "EG" = "#5BBCD6", "NMRI" = '#8A2BE2')
field_colors <- c('Brazil' = "#9986A5", 'Senegal' = "#79402E", 'Niger'="#CCBA72" ,'Tanzania'='#0F0D0E')

all_colors <- c(lab_colors, field_colors)

order <- c('BRE', 'EG', 'LE', 'OR', 'NMRI', 'Brazil', 'Niger', 'Senegal', 'Tanzania')

# Import the data ####
pixy_pi <- read.delim("single_gt_pi.txt")

pixy_pi <- pixy_pi %>% mutate(origin = ifelse(pop %in% c('BRE', 'EG', 'LE', 'OR', 'NMRI'),'Lab', 'Field'))

# Calculate mean
pixy_avg_pi <- pixy_pi %>%
  group_by(pop, origin) %>% 
  summarise(mean = mean(avg_pi, na.rm =T),  sd = sd(avg_pi, na.rm =T))

pixy_pi[c(1:5000),] %>% 
  group_by(origin) %>%
  shapiro_test(avg_pi)

pi_res <- pixy_pi %>%  ungroup() %>% wilcox_test(avg_pi ~ origin) %>% mutate(pop = 'EG', origin = 'Lab')

pixy_avg_pi$pop <- factor(pixy_avg_pi$pop, levels = order)
pixy_avg_pi$origin <- factor(pixy_avg_pi$origin, levels = c('Lab', 'Field'))

### Plot barplot ####
p_pi_barplot <- ggplot(pixy_avg_pi, aes(pop, mean )) +
  geom_errorbar(aes(x=pop, ymin = mean-0.0001, ymax = mean+sd), width = 0.5) +
  facet_grid(~factor(origin, levels = c('Lab', 'Field')), scales ='free') +
  geom_col(aes(fill = pop)) + 
  scale_fill_manual(values = all_colors) +
  theme_minimal() +
  labs(y=expression(bold('Nucleotide diversity'~(pi))), x='Population') +
  geom_text(data = pi_res, aes(label = paste0('Wilcoxon', ', p < 0.001')), y = 0.004, x =3) +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 18), axis.line = element_line(), legend.position = 'none')

# Plot ECDF ####
pi_field <- filter(temp_pi, origin == 'Field')
pi_lab <- filter(temp_pi, origin == 'Lab')

ks <- 
  ks.test(pi_field$avg_pi, pi_lab$avg_pi) 

ggplot(temp_pi, aes(avg_pi, color = pop, linetype = origin)) + 
  stat_ecdf() +
  scale_color_manual(values = all_colors) + 
  labs(x = "Avg Pi / 25kb window", y = "ECDF", color = 'Population', linetype = 'Origin') +
  guides(color = guide_legend(nrow = 2, byrow = T), linetype = guide_legend(nrow = 2)) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 16), axis.line = element_line(), legend.box.just = 'center', 
        legend.position = 'bottom') +
  annotate('text', x = 0.05, y =0.75, label = 'Kolmogorov-Smirnov test, p < 0.001')


# Export plots ####
pdf('Figure3_nucleotide_diversity.pdf', width = 8, height = 10)

(p_het & theme(axis.title.x = element_blank())) / p_pi_barplot  +
  plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold'))

dev.off()

