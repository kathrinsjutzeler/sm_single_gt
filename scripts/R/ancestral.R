### Derived Variants 7/11/2025 ####
library(tidyverse)
library(rstatix)

setwd('D:/Documents/Anderson_lab/Diversity_paper/Analysis/post_review/ancest')

lab_freq <- read.delim("lab_freq.txt", header = FALSE)
colnames(lab_freq) <- c('chrom', 'pos', 'freq_lab')
field_freq <- read.delim("field_freq.txt", header = FALSE)
colnames(field_freq) <- c('chrom', 'pos', 'freq_field')
outgroup_freq <- read.delim("outgroup_freq.txt", header = FALSE)
colnames(outgroup_freq) <- c('chrom', 'pos', 'freq_outgroup')

freq <- list.files(pattern = 'txt')

freqs <- lapply(freq[c(1:3, 6:12)], read.delim, header = FALSE)

names(freqs) <- c('Brazil', 'BRE', 'EG', 'LE', 'Niger', 'NMRI', 'OR', 'Outgroup', 'Senegal', 'Tanzania')

df <- bind_rows(freqs, .id = "pop")

colnames(df) <- c('pop', 'chrom', 'pos', 'AF')

df2 <- df %>%
  spread(pop, AF)

BRE <- df2 %>%
  dplyr::select(chrom, pos, BRE, Outgroup) %>% filter(Outgroup == 0, BRE >= 1)

LE <- df2 %>%
  dplyr::select(chrom, pos, LE, Outgroup) %>% filter(Outgroup == 0, LE >= 1)

OR <- df2 %>%
  dplyr::select(chrom, pos, OR, Outgroup) %>% filter(Outgroup == 0, OR >= 1)

NMRI <- df2 %>%
  dplyr::select(chrom, pos, NMRI, Outgroup) %>% filter(Outgroup == 0, NMRI >= 1)

EG <- df2 %>%
  dplyr::select(chrom, pos, EG, Outgroup) %>% filter(Outgroup == 0, EG >= 1)

br <- df2 %>%
  dplyr::select(chrom, pos, Brazil, Outgroup) %>% filter(Outgroup == 0, Brazil >= 1)

sn <- df2 %>%
  dplyr::select(chrom, pos, Senegal, Outgroup) %>% filter(Outgroup == 0, Senegal >= 1)

tz <- df2 %>%
  dplyr::select(chrom, pos, Tanzania, Outgroup) %>% filter(Outgroup == 0, Tanzania >= 1)

ng <- df2 %>%
  dplyr::select(chrom, pos, Niger, Outgroup) %>% filter(Outgroup == 0, Niger >= 1)


df3 <- bind_rows(OR, NMRI, BRE, EG, LE, br, sn, tz, ng)

df4 <- df3 %>% dplyr::select(-Outgroup) %>% gather("pop", "af", 3:11)

chrom.labs <- c('1', '2', '3', '4', '5', '6', '7')
names(chrom.labs) <- levels(factor(df3$chrom))

# Define colors
lab_colors = c("OR" = "#B22222", "BRE" = "#00A08A", "LE" = "#F98400", "EG" = "#5BBCD6", "NMRI" = '#8A2BE2')
field_colors <- c('Brazil' = "#9986A5", 'Senegal' = "#79402E", 'Niger'="#CCBA72" ,'Tanzania'='#0F0D0E')

all_colors <- c(lab_colors, field_colors)

order <- c('BRE', 'EG', 'LE', 'OR', 'NMRI', 'Brazil', 'Niger', 'Senegal', 'Tanzania')

df4$pop <- factor(df4$pop, levels = order)

df4 <- df4 %>% mutate(origin = ifelse(pop %in% c('BRE', 'EG', 'LE', 'OR', 'NMRI'),'Laboratory', 'Field'))

df4$origin <- factor(df4$origin, levels = c('Laboratory', 'Field'))

counts2 <- df4 %>%
  group_by(pop, origin) %>%
  na.omit() %>%
  count()

derived_res <- counts2 %>% ungroup() %>% wilcox_test(n ~ origin) %>% mutate(pop = 'EG', origin = 'Laboratory')

fixed_SNV <- ggplot(counts2, aes(pop, n,  fill=pop )) +
  geom_col() +
  scale_fill_manual(values = all_colors) +
  facet_grid(~factor(origin, levels = c('Laboratory', 'Field')), scales ='free') +
  theme_minimal() +
  labs(y="Number of fixed SNVs", x='Population') +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 16, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 16), axis.line = element_line(), legend.position = 'none') +
  geom_text(data = derived_res, aes(label = paste0('Wilcoxon test', ', p = ',round(p,3)), y = 2000))


wilcox.test(n ~ origin, data = counts2)


write.csv(der, 'derived2.csv')
