# Pi statistics for Sm single gentoype project
# Author: Kathrin Bailey
# Date: October 13, 2022
# Last updated: February 8, 2024
# R version 4.2.0, tidyverse version 1.3.2

library(ggridges)
library(tidyverse)
library(stats)
library(rstatix)
library(patchwork)


# Define colors
lab_colors = c("OR" = "#B22222", "BRE" = "#00A08A", "LE" = "#F98400", "EG" = "#5BBCD6", "NMRI" = '#8A2BE2')
field_colors <- c('Brazil' = "#9986A5", 'Senegal' = "#79402E", 'Niger'="#CCBA72" ,'Tanzania'='#0F0D0E')

all_colors <- c(lab_colors, field_colors)

order <- c('BRE', 'EG', 'LE', 'OR', 'NMRI', 'Brazil', 'Niger', 'Senegal', 'Tanzania')

# Import the data ####
pixy_pi <- read.delim("pi/single_gt_pi.txt")

pixy_pi <- pixy_pi %>% mutate(origin = ifelse(pop %in% c('BRE', 'EG', 'LE', 'OR', 'NMRI'),'Laboratory', 'Field'))

# Calculate mean
pixy_avg_pi <- pixy_pi %>%
  group_by(pop, origin) %>% 
  summarise(mean = mean(avg_pi, na.rm =T),  sd = sd(avg_pi, na.rm =T))

pixy_pi[c(1:5000),] %>% 
  group_by(origin) %>%
  shapiro_test(avg_pi)

#pi_res <- pixy_pi %>%  ungroup() %>% wilcox_test(avg_pi ~ origin) %>% mutate(pop = 'EG', origin = 'Lab')

pixy_avg_pi$pop <- factor(pixy_avg_pi$pop, levels = order)
pixy_avg_pi$origin <- factor(pixy_avg_pi$origin, levels = c('Laboratory', 'Field'))

# Plot pi barplot ####
p_pi_barplot <- ggplot(pixy_avg_pi, aes(pop, mean )) +
  geom_errorbar(aes(x=pop, ymin = mean-0.0001, ymax = mean+sd), width = 0.5) +
  facet_grid(~factor(origin, levels = c('Laboratory', 'Field')), scales ='free') +
  geom_col(aes(fill = pop)) + 
  scale_fill_manual(values = all_colors) +
  theme_minimal() +
  labs(y=expression(bold('Nucleotide diversity'~(pi))), x='Population') +
  geom_text(data = pi_res, aes(label = paste0('Kolmogorov-Smirnov', ', p < 0.001')), y = 0.004, x =3) +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 18), axis.line = element_line(), legend.position = 'none')


pixy_pi$pop <- factor(pixy_pi$pop, levels = order)
pixy_pi$origin <- factor(pixy_pi$origin, levels = c('Laboratory', 'Field'))

# Plot pi boxplot ####
p_pi_boxplot <- ggplot(pixy_pi, aes(pop, avg_pi, fill = pop)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(~origin, scales = 'free') +
  #facet_grid(~factor(origin, levels = c('Laboratory', 'Field')), scales ='free') +
  scale_fill_manual(values = all_colors) +
  ylim(0, 0.003) +
  #stat_summary(fun.y=mean, geom="line", size=14, color="red")
  theme_minimal() +
  labs(y=expression(bold('Average nucleotide diversity'~(pi))), x='Population') +
  geom_text(data = subset(pixy_pi, origin == 'Laboratory'), aes(label = paste0('Kolmogorov-Smirnov', ', p < 0.001')), y = 0.003, x =3) +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 16, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 16), axis.line = element_line(), legend.position = 'none')

# Plot pi across genome ####

lab_pi <- read.delim("pi/single_gt_lab_pi.txt")

lab_pi <- lab_pi %>% mutate(origin = ifelse(pop %in% c('BRE', 'EG', 'LE', 'OR', 'NMRI'),'Lab', 'Field'))
lab_pi$pop <- factor(lab_pi$pop, levels = order)

lab_pi <- lab_pi %>% mutate(smooth = runmed(avg_pi, 21))

t <- lab_pi %>%
  distinct(smooth, .keep_all = T)

chrom.labs <- c('1', '2', '3', '4', '5', '6', '7')
names(chrom.labs) <- levels(factor(lab_pi$chromosome))

p_pi_across <- ggplot(lab_pi, aes(window_pos_1, avg_pi, color = pop)) +
  geom_point(size = 1, alpha = 0.5) +
  facet_grid(pop~chromosome, switch = "x", scales = 'free_x', space = 'free_x',
             labeller = labeller(chromosome = chrom.labs)) +
  scale_color_manual(values = all_colors) +
  geom_smooth(method = 'loess', se = F, color = 'black') +
  #scale_y_log10() +
  labs(y = expression(bold('Average nucleotide diversity' ~(pi))), x = 'Chromosome') +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.line.x = element_line(), 
        axis.title = element_text(size = 16, face = 'bold'), text = element_text(size =12),
        axis.line.y = element_line(), panel.grid = element_blank(), legend.position = 'none') 
    

# Plot SNP density across genome ####
density_files <- list.files('./SNPdensity')

snpden <- lapply(paste0("./SNPdensity/",density_files), read.delim)

names(snpden) <- c('BRE', 'EG', 'LE', 'NMRI', 'OR')

snpden_df <- bind_rows(snpden, .id = 'pop')
snpden_df <- snpden_df %>% filter(!CHROM %in% c('SM_V10_Z', 'SM_V10_WSR', 'SM_V10_MITO'))

snpden_df <- snpden_df %>% group_by(pop, CHROM) %>% mutate(smooth = runmed(VARIANTS.KB, 11))

p_snp_density <- ggplot(snpden_df, aes(BIN_START,smooth, color = pop)) +
  geom_point(size = 1, alpha = 0.5) +
  facet_grid(pop~CHROM, scales = 'free_x', switch = "x") +
  scale_color_manual(values = all_colors) +
  geom_smooth(method = 'loess', se = F, color = 'black') +
  theme_minimal() +
  theme(axis.text.x = element_blank(), legend.position = 'none', panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(y = "SNP density", x = 'Genomic position')


# Plot ECDF ####
pi_field <- filter(pixy_pi, origin == 'Field', pop != 'Tanzania')
pi_lab <- filter(pixy_pi, origin == 'Laboratory',  pop != 'BRE')

#ks <- 
  ks.test(pi_field$avg_pi, pi_lab$avg_pi) 

pixy_pi$pop <- factor(pixy_pi$pop, levels = order)

p_pi_ECDF <- ggplot(pixy_pi, aes(avg_pi, color = pop, linetype = origin)) + 
  stat_ecdf() +
  scale_color_manual(values = all_colors) + 
  #facet_grid(pop~ chromosome, scale = 'free_x')
  labs(x = "Avg Pi / 25kb window", y = "ECDF", color = 'Population', linetype = 'Origin') +
  guides(color = guide_legend(nrow = 2, byrow = T), linetype = guide_legend(nrow = 2)) +
  theme_minimal() +
  guides(color = guide_legend(nrow = 2, byrow = T), linetype = guide_legend(nrow = 2, 
                                                                            override.aes = list(color = 'black'))) +
  
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 12), axis.line = element_line(), legend.box.just = 'center', 
        legend.position = 'bottom') +
  annotate('text', x = 0.05, y =0.75, label = 'Kolmogorov-Smirnov test, p < 0.001')


# Export plots ####
#pdf('Figure3_nucleotide_diversity.pdf', width = 8, height = 10)

#(p_het & theme(axis.title.x = element_blank())) / p_pi_barplot  +
#  plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold'))

#dev.off()

jpeg('Figure3_nucleotide_diversity.jpg', width = 8, height = 10, units = 'in', res = 300)

(p_pi_across / p_pi_boxplot)   +
    plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold', size =18))

dev.off()

pdf('FigureS_ECDF_nucleotide_diversity.pdf', width = 8, height = 6)

(p_pi_ECDF)

dev.off()

pdf('FigureS_SNP_density.pdf', width = 8, height = 6)

(p_snp_density)

dev.off()

jpeg('FigureS_Pi_across.jpg', width = 8, height = 6, units = 'in', res = 300)

(p_pi_boxplot)

dev.off()

ggsave('pi_across_smooth.png', width = 8, height = 6)
