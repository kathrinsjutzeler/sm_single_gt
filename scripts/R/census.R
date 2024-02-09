# Cercarial shedding 
# Author: Kathrin Jutzeler
# Date: May 15, 2023
# Last updated
# R version 4.2.0, tidyverse version 1.3.2, ggplot2 3.3.6      ✔ purrr   0.3.4 
#✔ tibble  3.1.8      ✔ dplyr   1.0.10
#✔ tidyr   1.2.0      ✔ stringr 1.4.1 
#✔ readr   2.1.2      ✔ forcats 0.5.2 

library(tidyverse)
library(readxl)
library(rstatix)
library(ggpubr)


# Import data ####
worms <- readxl::read_xlsx('Archive/Snail-rodent infection log.xlsx', sheet = 'Infections')

worms$`Due date...5` <- as.numeric(worms$`Due date...5`)
worms$`Due date...5` <- as.Date(worms$`Due date...5`, origin = "1900/01/01")

worms$`Date of infection...2` <- as.numeric(worms$`Date of infection...2`)
worms$`Date of infection...2` <- as.Date(worms$`Date of infection...2`, origin = "1900/01/01")

worms_df <- worms %>%
  filter(`Rodent type` == "Hamster") %>%
  mutate(year = ifelse(`Date of infection...2` < '2016-01-01' & `Date of infection...2` > '2014-12-31', '2015', 
                       ifelse(`Date of infection...2` < '2017-01-01' & `Date of infection...2` > '2015-12-31', '2016',
                       ifelse(`Date of infection...2` < '2018-01-01' & `Date of infection...2` > '2016-12-31', '2017',
                       ifelse(`Date of infection...2` < '2019-01-01' & `Date of infection...2` > '2017-12-31', '2018',
                       ifelse(`Date of infection...2` < '2020-01-01' & `Date of infection...2` > '2018-12-31', '2019',
                       ifelse(`Date of infection...2` < '2021-01-01' & `Date of infection...2` > '2019-12-31', '2020',
                       ifelse(`Date of infection...2` < '2022-01-01' & `Date of infection...2` > '2020-12-31', '2021',
                       ifelse(`Date of infection...2` < '2023-01-01' & `Date of infection...2` > '2021-12-31', '2022',
                       ifelse(`Date of infection...2` > '2022-12-31',"2023", "other"))))))))),
)

cen <- worms_df %>%
  filter(`Sm strain` %in% c('SmLE', 'SmEG', 'SmOR', 'SmBRE'))

cen2 <- cen %>%
  dplyr::select(1,4:7, 8, year)

cen2$year <- as.numeric(cen2$year)
cen2$`Number of miracidia` <- as.numeric(cen2$`Number of miracidia`)
cen2$`Number of surviving snails` <- as.numeric(cen2$`Number of surviving snails`)

cen2 <-  filter(cen2, year >= 2016, `Number of miracidia` !=1)

cen2 <- cen2 %>% 
  mutate(gt_high = `Number of miracidia` * `Number of snails infected` , 
         gt_low = 1 * `Number of snails infected`,
         IR = `Number of snails infected` / `Number of surviving snails`,
         lambda = -log(1-IR))


# Plot poisson distribution ####

## Data processing #####
for (i in c(0:5)) {
  col_name <- as.character(i)
  cen2[[col_name]] <- dpois(i, cen2$lambda)
}

for (i in c(6, 8:11)) {
  col_name <- as.character(i)
  cen2[[col_name]] <- ifelse(cen2$`Number of miracidia` == i, dpois(i, cen2$lambda), NA)
}

pois2 <- gather(cen2, 'miracidia', 'probability', 12:22)

pois2$miracidia <- as.numeric(pois2$miracidia)

pois3 <- pois2 %>%
  mutate(snails = probability*`Number of surviving snails`) %>%
  mutate(gt = snails * miracidia)

pois4 <- pois3 %>%
  group_by(`Due date...5`, `Sm strain`) %>%
  summarize(sum = sum(gt, na.rm = T)) %>%
  filter(!sum == 0)

## Plot across time #####
p_census <- ggplot(pois4, aes(`Due date...5`, sum, color = `Sm strain`)) +
  geom_line() +
  geom_point() + 
  #geom_smooth(se =F) +
  theme_minimal() +
  scale_color_manual(values = c('SmOR' = "#B22222", 'SmBRE' = "#00A08A", 'SmLE'= "#F98400", 'SmEG'= "#5BBCD6"),
                     labels = c('BRE', 'EG', 'LE', 'OR')) +
  labs(x = 'Shedding date', y = 'No. of possible genotypes') +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 18), axis.line = element_line(), legend.position = 'bottom') +
  labs(color = 'Population')

ggsave('Plots/Poisson_census.png', width = 8, height =  6)

# Get harmonic mean
hmean_df <- pois4 %>%
  group_by(`Sm strain`) %>%
  reframe(h.mean = DescTools::Hmean(sum,conf.level = 0.95))

replacements <- c("SmBRE" = "BRE", "SmEG" = "EG", 'SmOR' = 'OR', 'SmLE' = 'LE')

hmean_df$`Sm strain` <- str_replace_all(hmean_df$`Sm strain`, replacements)
hmean_df$label <- rep(c('mean', 'CI (L)', 'CI (U)'),4)

hmean_df <- hmean_df %>% spread(label, h.mean)

write.csv(hmean_df, 'census_df.csv')

## Plot harmonic mean #####
p_census_hm <- ggplot(hmean_df, aes(`Sm strain`,mean, fill=`Sm strain` )) +
  geom_errorbar(aes(ymin=10, ymax=`CI (U)`, width = 0.5)) + 
  geom_col() + 
  scale_fill_manual(values = lab_colors[1:4]) +
  theme_minimal() +
  labs(y="Census population size", x='Population') +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 18), axis.line = element_line(), legend.position = 'none')


ggsave('Plots/Poisson_census_barplot.png', width = 8, height =  6)

# Export plots ####

pdf('FigureS_census.pdf', width = 8, height = 8)

(p_census/ (p_census_hm + p_ne_cor)) +
  plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold', size = 18))

dev.off()





