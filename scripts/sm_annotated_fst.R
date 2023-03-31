# Fst statistics for Aim 1 - Sm Single GT project
# Author: Kathrin Bailey
# Date: February 16, 2023
# Last updated
# R version 4.2.0, tidyverse version 1.3.2, ggplot2 3.3.6      ✔ purrr   0.3.4 
#✔ tibble  3.1.8      ✔ dplyr   1.0.10
#✔ tidyr   1.2.0      ✔ stringr 1.4.1 
#✔ readr   2.1.2      ✔ forcats 0.5.2 

library(ggridges)
library(tidyverse)
library(stats)

# Import the data
fst_EGvBRE <-read.delim("annotated-EGvBRE.fst", header = FALSE)
fst_EGvBRE$comparison <- "EG_v_BRE"
fst_EGvOR <-read.delim("annotated-EGvOR.fst", header = FALSE)
fst_EGvOR$comparison <- "EG_v_OR"
fst_EGvLE <-read.delim("annotated-EGvLE.fst", header = FALSE)
fst_EGvLE$comparison <- "EG_v_LE"
fst_BREvLE <-read.delim("annotated-BREvLE.fst", header = FALSE) 
fst_BREvLE$comparison <- "BRE_v_LE"
fst_ORvLE <-read.delim("annotated-ORvLE.fst", header = FALSE) 
fst_ORvLE$comparison <- "OR_v_LE"
fst_ORvBRE <-read.delim("annotated-ORvBRE.fst", header = FALSE)
fst_ORvBRE$comparison <- "OR_v_BRE"

fst_EG <- read.delim("EGannotated.fst", header = FALSE)
fst_EG$comparison <- "EG"
fst_BRE <- read.delim("BREannotated.fst", header = FALSE)
fst_BRE$comparison <- "BRE"
fst_OR <- read.delim("ORannotated.fst", header = FALSE)
fst_OR$comparison <- "OR"
fst_LE <- read.delim("LEannotated.fst", header = FALSE)
fst_LE$comparison <- "LE"

df_list <- list(fst_EGvBRE, fst_EGvLE, fst_EGvOR, fst_BREvLE, fst_ORvBRE,fst_ORvLE, fst_EG, fst_OR, fst_LE, fst_BRE)

f_setup_df <- function(data) {
  
  data <- data %>%
    separate(V10, c("ID", "coverage", "sequence_id", "ORF", "extraCN", "CN"), sep = ';') %>%
    select(1:6,8:9, ID, comparison) 
  data$ID <- str_remove(data$ID, "ID=")
  
  setNames(data.frame(data), c("CHROM", "START", "END", "N_VARIANTS", "WEIGHTED_FST","MEAN_FST",'GENE_START', 'GENE_END',
                               "ID", "comparison"))
}

all_data <- lapply(df_list, f_setup_df)

data <- bind_rows(all_data)

# Add position to make plotting easier
#data <- data %>%
#  group_by(comparison) %>%
#  mutate(position = 1:n())

# Put chromosomes in actual order
data$CHROM <- factor(data$CHROM, 
                             levels= c("HE601624.3", "HE601625.3", "HE601626.3", "HE601627.3", "HE601628.3", 
                                       "HE601629.3", "HE601630.3", "OU426849.1", "OU426848.1", "OU426850.1",
                                       "HE601631.3"))

# Rename chromosomes
chrom.labs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "PAR1", "Z", "PAR2", "W")
names(chrom.labs) <- c("HE601624.3", "HE601625.3", "HE601626.3", "HE601627.3", "HE601628.3", 
                       "HE601629.3", "HE601630.3", "OU426849.1", "OU426848.1", "OU426850.1",
                       "HE601631.3")


# This part is to smooth the data
data <- data %>%
  group_by(comparison) %>%
  mutate(smooth_fst = runmed(MEAN_FST, 11))

#get median for intercept
#data %>%
#  group_by(comparison) %>%
#  summarize(median = median(MEAN_FST))#median(sort(data$MEAN_FST))

## A tibble: 6 × 2
#comparison median
#<chr>       <dbl>
#1 BRE_v_LE   0.0558
#2 EG_v_BRE   0.202 
#3 EG_v_LE    0.219 
#4 EG_v_OR    0.158 
#5 OR_v_BRE   0.183 
#6 OR_v_LE    0.195

# get 99th percentile
threshold <- data %>%
  group_by(comparison) %>%
  summarize(quantile = quantile(smooth_fst, probs = .99))

# A tibble: 6 × 2
#comparison quantile
#<chr>         <dbl>
#1 BRE_v_LE      0.270
#2 EG_v_BRE      0.495
#3 EG_v_LE       0.552
#4 EG_v_OR       0.419
#5 OR_v_BRE      0.460
#6 OR_v_LE       0.522

#Add threshold to the data table
data <- data %>% left_join(threshold, by = 'comparison')

# list of single populations (within comparison)
within_data <- c("EG", "OR", "BRE", "LE")

# Plot by chromosome
p_smooth_99 <-  
  ggplot(subset(data, comparison!=c("EG", "OR", "BRE", "LE")), aes(GENE_START,y = smooth_fst, color=CHROM )) +
  geom_point(size = 0.5,alpha = .75) +
  facet_grid(comparison~CHROM, switch = "x", scales = "free_x", labeller = labeller(CHROM = chrom.labs)) +
  theme_minimal() +
  geom_hline(aes(yintercept=quantile), linetype="dashed", size=0.75) +
  scale_colour_cyclical(values = c("#3030D0", "#D12959")) +
  labs(x= "Genomic Position", y = "Mean Fst") +
  theme(axis.text.x = element_blank(), panel.grid = element_blank())

ggsave("pretty_fst.png", width = 8, height = 6)

# Boxplot 
ggplot(subset(data, !(comparison %in% within_data)), aes(comparison,y = smooth_fst, color=comparison )) +
  geom_boxplot() + 
  scale_color_manual(values = comp_colors) +
  theme_minimal() +
  theme(legend.position="none") +
  stat_compare_means() +
  labs(x= "Population", y = "Mean Fst") 

# Subset all values above 99th percentile
subset <- data %>%  
  filter(smooth_fst >= quantile)

p <- ggplot(subset, aes(START,y = smooth_fst, color=CHROM )) +
  geom_point(size = 0.5,alpha = .75) +
  facet_grid(comparison~CHROM, switch = "x", scales = "free_x", labeller = labeller(CHROM = chrom.labs)) +
  theme_minimal() +
  geom_hline(aes(yintercept =quantile), linetype="dashed", size =0.5) +
  #geom_text(aes(label = ID),check_overlap = TRUE) +
  scale_colour_cyclical(values = c("#3030D0", "#D12959")) +
  labs(x= "Genomic Position", y = "Mean Fst") +
  scale_x_continuous(breaks=seq(0, 4e+8, 1E+8))

# SmOR-SULT
Smp_089320 <- data %>%
    filter(ID == "Smp_089320")

ggplot(Smp_089320, aes(GENE_START,y = smooth_fst, color=CHROM )) +
  geom_point(size = 1,alpha = .75) +
  facet_grid(comparison~CHROM, switch = "x", scales = "free_x", labeller = labeller(CHROM = chrom.labs)) +
  theme_minimal() +
  #geom_hline(aes(yintercept =quantile), linetype="dashed", size =0.5) +
  #geom_text(aes(label = ID),check_overlap = TRUE) +
  scale_colour_cyclical(values = c("#3030D0", "#D12959")) +
  labs(x= "Genomic Position", y = "Mean Fst", title = "SmOR-SULT") +
  theme(panel.grid = element_blank(), axis.text.x = element_blank())

# Export subset table
export <- subset %>%
  group_by(CHROM, ID) %>%
  filter(!duplicated(ID))

write_excel_csv2(export, "SNPs_above_threshold.csv")

### Look at antigens
# Set colors
comp_colors <- c("#D8B70A", "#3B9AB2", "#A2A475", "#81A88D", "#972D15", "#F4B5BD")

# Isolate genes
IPSE_genes <- c("Smp_112110", "Smp_334350", "Smp_334690", "Smp_335420", "Smp_335430", "Smp_335440", "Smp_335450")
kappa_genes <- c("Smp_335460", "Smp_335470", "Smp_335480", "Smp_335490", "Smp_335510", "Smp_344300")
omega_genes <- c("Smp_158430", "Smp_333870", "Smp_333930", "Smp_334170", "Smp_334240", "Smp_345790")

IPSE <- data %>%
  filter(ID %in% IPSE_genes)

ggplot(IPSE, aes(GENE_START,y = smooth_fst, color=comparison, shape = ID)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = comp_colors) +
  labs(x= "Genomic Position", y = "Mean Fst", title = "IPSE") +
  theme(panel.grid = element_blank(), axis.text.x = element_blank()) 

kappa <- data %>%
  filter(ID %in% kappa_genes)

ggplot(kappa, aes(GENE_START,y = smooth_fst, color=comparison, shape = ID)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = comp_colors) +
  labs(x= "Genomic Position", y = "Mean Fst", title = "kappa-5") +
  theme(panel.grid = element_blank(), axis.text.x = element_blank()) 

omega <- data %>%
  filter(ID %in% omega_genes)

ggplot(omega, aes(GENE_START,y = smooth_fst, color=comparison, shape = ID)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = comp_colors) +
  labs(x= "Genomic Position", y = "Mean Fst", title = "Omega-1") +
  theme(panel.grid = element_blank(), axis.text.x = element_blank()) 

### Look at within populations

  