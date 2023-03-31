# Admixture plots for Aim 1 - Sm Single GT project
# Author: Kathrin Bailey
# Date: February 16, 2023
# Last updated
# R version 4.2.0, tidyverse version 1.3.2library

library(ggsci)
library(tidyverse)
library(stats)


# Admixture plot
tbl=read.table("final_lab_nuclear_DNA.4.Q")
tbl=read.table("lab_nuclear_nosex.4.Q")
meta <- read.table("sample_order.txt")
meta <- read.table("sample_order_nosex.txt")

# Add samples to groups
OR <- c("5", "1", "6", "2", "24", "7", "25", "14", "26", "17")
BRE <- c("27", "32", "51", "38", "68", "48", "72", "49", "80", "55")
LE <- c("84", "92", "86", "104", "90", "107", "97", "109", "112", "111")
EG <- c("138", "130", "154", "135", "161", "142","167", "146", "169", "163")        


samples <- as.data.frame(str_split(meta$V1, "_", simplify = TRUE))

meta <- samples %>%
  dplyr::select(V1) %>%
  mutate(pop = ifelse(samples$V1 %in% OR, "OR",
                                          ifelse(samples$V1 %in% BRE, "BRE",
                                                 ifelse(samples$V1 %in% LE, "LE", "EG"))))

admix <- cbind(tbl, meta)
colnames(admix) <- c("v1", "v2","v3",'v4', "sample", "pop")

# Create a dataframe listing sample, population, and all probabilities 
data <- gather(admix, "col", "prob", 1:4)


ggplot(subset(data, sample != "5"),aes(sample,prob, fill=col)) +
  geom_col(color = "gray", size = 0.1)+
  facet_grid(~fct_inorder(pop), switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(title = "K=4", y = "Ancestry", x = NULL) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 0.7)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid = element_blank(),
        text = element_text(size = 18)) +
  scale_fill_npg(guide = "none")

ggsave("Plots/admix.png", width = 8, height = 6)

#### lab vs field
tbl=read.table("labvfield_nuclear_CDS.7.Q")
meta <- read.table("samples.txt")

new_admix <- cbind(tbl, meta)

data <- gather(new_admix, "col", "prob", 1:7)

data <- data %>% 
  filter(!V2 %in% c('Uganda', 'Kenya', 'Caribbean'))
  
colnames(data) <- c('sample', "pop", 'col', 'prob')

data$pop <- factor(data$pop, levels = c("BRE", "EG", "LE", "OR", "Brazil", "Niger", "Senegal",
                                        "Tanzania"))

ggplot(subset(data, sample != '5_S29_L003'),aes(sample,prob, fill=col)) +
  #geom_col(color = "gray", size = 0.1) +
  geom_col(size = 0.1) +
  facet_grid(~pop, switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(title = "K=7", y = "Ancestry", x = NULL) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 0.7)) +
  theme(panel.spacing.x = unit(0.1, "lines"),
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid = element_blank(),
        text = element_text(size = 18)) +
  scale_fill_manual(values = c('red',"#BF812D", 'olivedrab3', '#D6604D', 'brown4', '#BCBDDC', 'gray37'), guide = 'none') # Q7
    #scale_fill_npg(guide = "none")

ggsave("Plots/fieldadmix_q7.png", width = 8, height = 6)

