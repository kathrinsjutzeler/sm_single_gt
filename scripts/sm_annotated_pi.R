# Pi statistics for Aim 1 - Sm Single GT project
# Author: Kathrin Bailey
# Date: February 16, 2023
# Last updated
# R version 4.2.0, tidyverse version 1.3.2

library(ggridges)
library(tidyverse)
library(stats)

# Import the data
pi_OR <- read.delim("annotated-OR.pi", header = FALSE)
pi_OR$pop <- "OR"
pi_EG <- read.delim("annotated-EG.pi", header = FALSE)
pi_EG$pop <- "EG"
pi_BRE <- read.delim("annotated-BRE.pi", header = FALSE)
pi_BRE$pop <- "BRE"
pi_LE <- read.delim("annotated-LE.pi", header = FALSE)
pi_LE$pop <- "LE"

pi_df_list <- list(pi_OR, pi_EG, pi_BRE, pi_LE)

f_setup_pi <- function(data) {
  
  data <- data %>%
    separate(V9, c("ID", "coverage", "sequence_id", "ORF", "extraCN", "CN"), sep = ';') %>%
    select(1:5,7:8, ID, pop) 
  data$ID <- str_remove(data$ID, "ID=")
  
  setNames(data.frame(data), c("CHROM", "START", "END", "N_VARIANTS", "PI", "GENE_START", "GENE_END","ID", "pop"))
}

pi_data_fixed <- lapply(pi_df_list, f_setup_pi)

pi_data <- bind_rows(pi_data_fixed)

# Put chromosomes in actual order
pi_data$CHROM <- factor(pi_data$CHROM, 
                                 levels= c("HE601624.3", "HE601625.3", "HE601626.3", "HE601627.3", "HE601628.3", 
                                           "HE601629.3", "HE601630.3", "OU426849.1", "OU426848.1", "OU426850.1",
                                           "HE601631.3"))

# Rename chromosomes
chrom.labs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "PAR1", "Z", "PAR2", "W")
names(chrom.labs) <- c("HE601624.3", "HE601625.3", "HE601626.3", "HE601627.3", "HE601628.3", 
                       "HE601629.3", "HE601630.3", "OU426849.1", "OU426848.1", "OU426850.1",
                       "HE601631.3")

pi_data$pop <- factor(pi_data$pop, 
                      levels= c("EG", "OR", "LE", "BRE"))

# Smoothing
pi_data <- pi_data %>% 
  group_by(pop) %>%
  mutate(smooth_pi = runmed(PI, 21))

# Get median for horizonal line
#pi_data_median <- pi_data %>%
#  group_by(pop) %>%
#  summarise(median = median(PI, na.rm = TRUE))

# A tibble: 4 × 2
#pop      median
#<chr>     <dbl>
#1 BRE   0.0000523
#2 EG    0.0000472
#3 LE    0.0000392
#4 OR    0.0000465

# get 99th percentile
pi_threshold <- pi_data %>%
  group_by(pop) %>%
  summarize(top = quantile(smooth_pi, probs = 0.99), bottom = quantile(smooth_pi, probs = 0.01))

# A tibble: 4 × 3
#pop        top     bottom
#<fct>    <dbl>      <dbl>
#1 EG    0.000250 0.00000737
#2 OR    0.000186 0.00000905
#3 LE    0.000172 0.00000653
#4 BRE   0.000191 0.0000123

#Add threshold to the data table
pi_data <- pi_data %>% left_join(pi_threshold, by = 'pop')

# Try to space out x axis by chromosome
#p_pretty_pi <-
  ggplot(pi_data, aes(START,y = smooth_pi, color=CHROM)) +
  geom_point(size = 0.5, alpha = .75) +
  facet_grid(pop~CHROM, switch = "x", scales = "free_x", labeller = labeller(CHROM = chrom.labs)) +
  theme_minimal() +
  geom_hline(aes(yintercept=top), linetype="dashed", size=0.75) +
  #geom_hline(aes(yintercept=bottom), linetype="dashed", size=0.75) +
  scale_colour_cyclical(values = c("#3030D0", "#D12959")) +
  labs(x= "Genomic Position", y = "Nucleotide Diversity (Pi)") +
  scale_x_continuous(breaks=seq(0, 5e+8, 1E+8)) +
  scale_y_log10() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank())
          #strip.text.x = element_blank())
  
# Boxplot 
ggplot(pi_data, aes(pop,y = PI, fill=pop )) +
  geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual(values = population_colors) +
  theme_minimal() +
  theme(legend.position="none") +
  stat_compare_means()

# Subset all values above 99th percentile
pi_subset <- pi_data %>%  
  filter(smooth_pi >= top)
  
ggplot(pi_subset, aes(GENE_START,y = smooth_pi, color=CHROM )) +
  geom_point(size = 0.5,alpha = .75) +
  facet_grid(pop~CHROM, switch = "x", scales = "free_x", labeller = labeller(CHROM = chrom.labs)) +
  theme_minimal() +
  geom_hline(aes(yintercept=top), linetype="dashed", size =0.5) +
  #geom_text(aes(label = ID),check_overlap = TRUE, size =3) +
  scale_colour_cyclical(values = c("#3030D0", "#D12959")) +
  labs(x= "Genomic Position", y = "Mean Fst") +
  scale_y_log10() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank())
#scale_x_continuous(breaks=seq(0, 4e+8, 1E+8))
  
# Export subset table
pi_export <- pi_subset %>%
  group_by(CHROM, ID) %>%
  filter(!duplicated(ID))

write_excel_csv2(pi_export, "Pi_SNPs_above_threshold.csv")

#SmOR-SULT
pi_Smp_089320 <- 
  pi_data %>%
  filter(ID == "Smp_089320") %>%
  group_by(pop, CHROM, GENE_START) %>%
  summarise(mean= mean(smooth_pi))

ggplot(pi_Smp_089320, aes(GENE_START,y = smooth_pi, color=CHROM )) +
  geom_point(size = 1,alpha = .75, position = position_dodge(1)) +
  facet_grid(pop~CHROM, switch = "x", scales = "free_x", labeller = labeller(CHROM = chrom.labs)) +
  theme_minimal() +
  geom_hline(aes(yintercept=top), linetype="dashed", size =0.5) +
  #geom_text(aes(label = ID),check_overlap = TRUE) +
  scale_colour_cyclical(values = c("#3030D0", "#D12959")) +
  labs(x= "Genomic Position", y = "Mean Fst", title = "SmOR-SULT") +
  scale_y_log10() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank())

## POPULATION COMPARISONS WITH RATIO
df_BRE <- pi_data_fixed[[3]]
df_LE <- pi_data_fixed[[4]]
df_OR <- pi_data_fixed[[1]]
df_EG <- pi_data_fixed[[2]]

BREvLE <- df_BRE %>%
  inner_join(df_LE, by= c("CHROM", 'START', 'END')) %>%
  mutate(pi_ratio = log(PI.x/PI.y), comparison = "BRE_v_LE") %>%
  select(CHROM, START, pi_ratio, comparison)

EGvBRE <- df_EG %>%
  inner_join(df_BRE, by= c("CHROM", 'START', 'END')) %>%
  mutate(pi_ratio = log(PI.x/PI.y), comparison = "EG_v_BRE") %>%
  select(CHROM, START, pi_ratio, comparison)

EGvLE <- df_EG %>%
  inner_join(df_LE, by= c("CHROM", 'START', 'END')) %>%
  mutate(pi_ratio = log(PI.x/PI.y), comparison = "EG_v_LE") %>%
  select(CHROM, START, pi_ratio, comparison)

EGvOR <- df_EG %>%
  inner_join(df_OR, by= c("CHROM", 'START', 'END')) %>%
  mutate(pi_ratio = log(PI.x/PI.y), comparison = "EG_v_OR") %>%
  select(CHROM, START, pi_ratio, comparison)

ORvBRE <- df_OR %>%
  inner_join(df_BRE, by= c("CHROM", 'START', 'END')) %>%
  mutate(pi_ratio = log(PI.x/PI.y), comparison = "OR_v_BRE") %>%
  select(CHROM, START, pi_ratio, comparison)

ORvLE <- df_OR %>%
  inner_join(df_LE, by= c("CHROM", 'START', 'END')) %>%
  mutate(pi_ratio = log(PI.x/PI.y), comparison = "OR_v_LE") %>%
  select(CHROM, START, pi_ratio, comparison)

pi_ratio <- bind_rows(BREvLE, EGvBRE, EGvLE, EGvOR, ORvBRE, ORvLE)

pi_ratio$CHROM <- factor(pi_ratio$CHROM, 
                        levels= c("HE601624.3", "HE601625.3", "HE601626.3", "HE601627.3", "HE601628.3", 
                                  "HE601629.3", "HE601630.3", "OU426849.1", "OU426848.1", "OU426850.1",
                                  "HE601631.3"))

pi_ratio <- pi_ratio %>%
  group_by(comparison) %>%
  mutate(smooth = runmed(pi_ratio, 31))

# get 99th percentile
pi_ratio_99 <- pi_ratio %>%
  group_by(comparison) %>%
  summarize(top = quantile(smooth, probs = 0.99), bottom = quantile(smooth, probs = 0.01))

# A tibble: 6 × 3
#comparison   top bottom
#<chr>      <dbl>  <dbl>
#1 BRE_v_LE    1.83 -0.794
#2 EG_v_BRE    1.28 -1.88 
#3 EG_v_LE     1.99 -1.77 
#4 EG_v_OR     1.53 -1.79 
#5 OR_v_BRE    1.32 -1.38 
#6 OR_v_LE     1.79 -1.39

#Add threshold to the data table
pi_ratio <- pi_ratio %>% left_join(pi_ratio_99, by = 'comparison')
  
pi_ratio_subset <- pi_ratio %>%  
  filter(smooth >= top | smooth <= bottom)

p +ggplot(pi_ratio_subset, aes(START,y = smooth, color=CHROM)) +
  geom_point(size = 0.5, alpha = .75) +
  facet_grid(comparison~CHROM, switch = "x", scales = "free_x", labeller = labeller(CHROM = chrom.labs)) +
  theme_minimal() +
  geom_hline(aes(yintercept=top), linetype="dashed", size=0.5) +
  geom_hline(aes(yintercept=bottom), linetype="dashed", size=0.5) +
  scale_colour_cyclical(values = c("#3030D0", "#D12959")) +
  labs(x= "Genomic Position", y = "log(Ratio of Nucleotide Diversity (Pi))") +
  scale_x_continuous(breaks=seq(0, 5e+8, 1E+8)) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank())


### Look at antigens
# Isolate genes
IPSE_genes <- c("Smp_112110", "Smp_334350", "Smp_334690", "Smp_335420", "Smp_335430", "Smp_335440", "Smp_335450")
kappa_genes <- c("Smp_335460", "Smp_335470", "Smp_335480", "Smp_335490", "Smp_335510", "Smp_344300")
omega_genes <- c("Smp_158430", "Smp_333870", "Smp_333930", "Smp_334170", "Smp_334240", "Smp_345790")

IPSE_pi <- pi_data %>%
  filter(ID %in% IPSE_genes)

ggplot(IPSE_pi, aes(GENE_START,y = smooth_pi, color=pop, shape = ID)) +
  geom_point(position = position_jitter(width = 4E2, height = 0)) +
  theme_minimal() +
  scale_color_manual(values = population_colors) +
  labs(x= "Genomic Position", y = "Nucleotide diversity (Pi)", title = "IPSE") +
  theme(panel.grid = element_blank(), axis.text.x = element_blank()) 

kappa_pi <- pi_data %>%
  filter(ID %in% kappa_genes)

ggplot(kappa_pi, aes(GENE_START,y = smooth_pi, color=pop, shape = ID)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = population_colors) +
  labs(x= "Genomic Position", y = "Nucleotide diversity (Pi)", title = "kappa-5") +
  theme(panel.grid = element_blank(), axis.text.x = element_blank()) 

omega_pi <- pi_data %>%
  filter(ID %in% omega_genes)

ggplot(omega_pi, aes(GENE_START,y = smooth_pi, color=pop, shape = ID)) +
  geom_point(position = position_jitter(w = 4E3, h =0)) +
  theme_minimal() +
  scale_color_manual(values = population_colors) +
  labs(x= "Genomic Position", y = "Nucleotide diversity (Pi)", title = "Omega-1") +
  theme(panel.grid = element_blank(), axis.text.x = element_blank()) 


####### lab v field
# Define colors 
all_colors <- c(population_colors, 'Brazil' = "olivedrab3", 'Senegal' = "red", 'Niger'="brown4" ,'Tanzania'='gray37')

# Import data (note these are CDS only!)
pi_OR <- read.delim("C:/Users/kbailey/Downloads/temp/OR.list_hardfilter-pi.windowed.pi")
pi_OR$pop <- "OR"
pi_EG <- read.delim("C:/Users/kbailey/Downloads/temp/EG.list_hardfilter-pi.windowed.pi")
pi_EG$pop <- "EG"
pi_BRE <- read.delim("C:/Users/kbailey/Downloads/temp/BRE.list_hardfilter-pi.windowed.pi")
pi_BRE$pop <- "BRE"
pi_LE <- read.delim("C:/Users/kbailey/Downloads/temp/LE.list_hardfilter-pi.windowed.pi")
pi_LE$pop <- "LE"
pi_br <- read.delim("C:/Users/kbailey/Downloads/temp/br.list_hardfilter-pi.windowed.pi")
pi_br$pop <- "Brazil"
pi_tz <- read.delim("C:/Users/kbailey/Downloads/temp/tz.list_hardfilter-pi.windowed.pi")
pi_tz$pop <- "Tanzania"
pi_sn <- read.delim("C:/Users/kbailey/Downloads/temp/sn.list_hardfilter-pi.windowed.pi")
pi_sn$pop <- "Senegal"
pi_ne <- read.delim("C:/Users/kbailey/Downloads/temp/ne.list_hardfilter-pi.windowed.pi")
pi_ne$pop <- "Niger"

pi_df_list <- list(pi_OR, pi_EG, pi_BRE, pi_LE, pi_br, pi_tz, pi_sn, pi_ne)

# Combine data frames
pi_data <- bind_rows(pi_df_list)
pi_data$pop <- factor(pi_data$pop, levels = c("BRE", "EG", "LE", "OR", "Brazil", "Senegal", "Niger", "Tanzania"))

# Plot boxplot 
ggplot(pi_data, aes(pop,y = PI, fill=pop )) +
  geom_boxplot() + 
  #scale_y_log10() +
  scale_fill_manual(values = all_colors)+
  theme_minimal() +
  labs(y='Nucleotide diversity (Pi)', x='Population') +
  theme(legend.position="none", text = element_text(size =18))

ggsave('Plots/labvfield_pi.png', width = 8, height = 6)

# Plot barplot 
avg <- pi_data %>%
  group_by(pop) %>%
  summarise(mean = mean(PI)) %>%
  mutate(origin = c(rep('Lab', 4), rep('Field',4)))

ggplot(avg, aes(pop,y = mean, fill=pop )) +
  geom_col() + 
  scale_fill_manual(values = all_colors) +
  facet_grid(~origin, scales ='free') +
  theme_minimal() +
  labs(y='Nucleotide diversity (Pi)', x='Population') +
  theme(legend.position="none", text = element_text(size =18),
        panel.grid = element_blank())

ggsave('Plots/labvfield_pi_bar.png', width = 8, height = 6)

# Use Neal's data
pi_neal <- read.delim('pi_neal.txt', sep = ' ', header = FALSE)

pi <- pi_neal %>%
  filter(!V1 %in% c('kenya', 'guadeloupe', 'oman', 'na', 'puerto_rico', 'uganda') )

ggplot(pi, aes(V1,V2, fill=V1 )) +
  geom_col() + 
  scale_fill_manual(values = all_colors)+
  facet_grid(~V3, scales ='free') +
  theme_minimal() +
  labs(y='Mean nucleotide diversity (Pi)', x='Population') +
    theme(legend.position="none", text = element_text(size =18), panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45))

ggsave('Plots/neal_pi.png', width = 8, height = 6)

