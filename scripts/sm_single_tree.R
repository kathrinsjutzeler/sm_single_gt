# Tree for Aim 1 - Sm Single GT project
# Author: Kathrin Bailey
# Date: February 16, 2023
# Last updated
# R version 4.2.0, tidyverse version 1.3.2

library(ape)
library(ggtree)
library(tidyverse)

population_colors = c("OR" = "#BF812D", "BRE" = "#FED976", "LE" = "#BCBDDC", "EG" = "#D6604D")
# Define colors 
field_colors <- c('Brazil' = "olivedrab3", 'Senegal' = "red", 'Niger'="brown4" ,'Tanzania'='gray37',
                'Caribbean' = "lightskyblue1", 'Egypt' = "steelblue3", 'Kenya'="palevioletred",
                'Uganda' ="aquamarine2", 'Outgroup' = 'pink', '0' = 'green')

#### TREE

tree <- ape::read.tree(paste0("outgroup_nuclear_DNA.min4.phy.treefile"))
tree <- ape::read.tree(paste0("filtered_nuclear_DNA.min4.phy.treefile"))


# Optionally set an outgroup.
outgroup <- c("Sro_female_1.1_CCATCCTC", 
              "Sro_female_1.2_CCGACAAC",
              "Sro_female_2.1_CCTAATCC",
              "Sro_female_2.2_CCTCTATC",
              "Sro_male_1.1_ATCATTCC",
              "Sro_male_1.2_ATTGGCTC",
              "Sro_male_2.1_CAAGGAGC",
              "Sro_male_2.2_CACCTTAC")

tree <- root(tree,outgroup = outgroup, resolve.root = T)

# Add samples to groups
OR <- c("5", "1", "6", "2", "24", "7", "25", "14", "26", "17")
BRE <- c("27", "32", "51", "38", "68", "48", "72", "49", "80", "55")
LE <- c("84", "92", "86", "104", "90", "107", "97", "109", "112", "111")
EG <- c("138", "130", "154", "135", "161", "142","167", "146", "169", "163")        

tree_samples <- tree$tip.label

tree_samples <- ifelse(grepl("^[1-9]", tree_samples), str_split(tree_samples, "_", simplify = TRUE), tree_samples)

tree$tip.label <- tree_samples

#Set group info
groupInfo <- split(tree$tip.label, ifelse(tree$tip.label %in% OR, "OR",
                                          ifelse(tree$tip.label %in% BRE, "BRE",
                                                 ifelse(tree$tip.label %in% LE, "LE",
                                                        ifelse(tree$tip.label %in% EG, "EG", "Outgroup")))))
                                          
tree <- groupOTU(tree, groupInfo)

p_nuclear_tree <- 
  ggtree(tree, aes(color=group), layout = 'circular') + 
  geom_tiplab(hjust = 0.1, size=3, align = TRUE) +
  geom_treescale() +
  labs(colour= "Populations", title = "Single GT infections - Nuclear DNA tree") +
  scale_color_manual(values = population_colors)

######## Lab vs field tree
# Tree with Uganda, Kenya, Guadeloupe
tree <- ape::read.tree(paste0("nuclear_out_field_CDS.min4.phy.treefile"))
# Tree w/o Uganda, Kenya, Guadeloupe
tree <- ape::read.tree(paste0("fortree.treefile"))

metadata <- read.csv('samples.csv')

f_label <- function(label) {
  x <- grep(label, metadata$sample)
  print(metadata[x,2])
}

# Optionally set an outgroup.
outgroup <- c("Sro_female_1.1_CCATCCTC", 
              "Sro_female_1.2_CCGACAAC",
              "Sro_female_2.1_CCTAATCC",
              "Sro_female_2.2_CCTCTATC",
              "Sro_male_1.1_ATCATTCC",
              "Sro_male_1.2_ATTGGCTC",
              "Sro_male_2.1_CAAGGAGC",
              "Sro_male_2.2_CACCTTAC")

groups <- lapply(tree$tip.label, f_label)

groups <- as.character(groups)

y <- which(groups %in% "c(\"Brazil\", \"Brazil\")")

groups <- replace(groups, y, 'Brazil')

groupInfo <- split(tree$tip.label, groups)

tree <- root(tree,outgroup = outgroup, resolve.root = T)

tree <- groupOTU(tree, groupInfo)

# To get node numbers
#ggtree(tree, aes(color=group)) + 
  #geom_tiplab(hjust = 0.1, size=3, align = TRUE) +
#  geom_text2(aes(subset=!isTip, label=node), hjust=1)

# Plot tree
p <- 
  ggtree(tree, aes(color=group)) + #, layout = 'circular') + 
  #geom_tiplab(hjust = 0.1, size=3, align = TRUE) +
  #geom_text2(aes(subset =!isTip, label=label), hjust=-.2) +
  geom_treescale() +
  #labs(colour= "none") +
  scale_color_manual(values = field_colors, limits = force) +
  theme(legend.position="left", text = element_text(size =18)) 
  #geom_hilight(node = c(219, 245, 336), alpha = .6, color = 'black', fill ="gray")
  #geom_hilight(node = c(211, 208), alpha = .6, color = 'black', fill ="gray")

# Collapse nodes
p2 <- collapse(p, 211, 'max', fill = 'steelblue3' ) %>%
  collapse(328, 'max', fill = 'brown4' )  %>%
  collapse(337, 'max', fill = 'red' ) %>%
  collapse(362, 'max', fill = 'lightskyblue1' ) %>%
  collapse(280, 'max', fill = 'olivedrab3' ) %>%
  collapse(207, 'max', fill= 'olivedrab3')
  
# Add bootstrap values 
p2 + geom_text2(aes(subset=!isTip & group != "Tanzania" & group != "Outgroup", label=label), hjust=0, vjust =1.25) 

# Save
ggsave('Plots/labvfield_tree_simp.png', width = 13, height = 9)

### BRE v LE parentals
tree <- ape::read.tree(paste0("parentals.treefile"))

metadata <- read.csv('samples.csv')

f_label <- function(label) {
  x <- grep(label, metadata$sample)
  print(metadata[x,5])
}

# Optionally set an outgroup.
outgroup <- c("Sro_female_1.1_CCATCCTC", 
              "Sro_female_1.2_CCGACAAC",
              "Sro_female_2.1_CCTAATCC",
              "Sro_female_2.2_CCTCTATC",
              "Sro_male_1.1_ATCATTCC",
              "Sro_male_1.2_ATTGGCTC",
              "Sro_male_2.1_CAAGGAGC",
              "Sro_male_2.2_CACCTTAC")

groups <- lapply(tree$tip.label, f_label)

groups <- as.character(groups)

groupInfo <- split(tree$tip.label, groups)

tree <- root(tree,outgroup = outgroup, resolve.root = T)

tree <- groupOTU(tree, groupInfo)

p <- ggtree(tree, aes(color=group)) + 
  #geom_tiplab(hjust = 0.1, size=3, align = TRUE) +
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.3)
  geom_treescale() +
  labs(colour= "Population") +
  scale_color_manual(values = c("BRE" = "#FED976", "LE" = "#BCBDDC", 'Outgroup' ='pink')) +
  theme(legend.position="bottom", text = element_text(size =18))

collapse(p, node = 20, 'max', fill = "#FED976") %>%
  collapse(22, 'min', fill='#BCBDDC')   

ggsave('Plots/parental_tree.png', width = 8, height = 6)

### BRE v LE all samples 
tree <- ape::read.tree(paste0("BRE_LE_outgroup.treefile"))

metadata <- read.csv('samples.csv')

f_label <- function(label) {
  x <- grep(label, metadata$sample)
  print(metadata[x,5])
}

# Optionally set an outgroup.
outgroup <- c("Sro_female_1.1_CCATCCTC", 
              "Sro_female_1.2_CCGACAAC",
              "Sro_female_2.1_CCTAATCC",
              "Sro_female_2.2_CCTCTATC",
              "Sro_male_1.1_ATCATTCC",
              "Sro_male_1.2_ATTGGCTC",
              "Sro_male_2.1_CAAGGAGC",
              "Sro_male_2.2_CACCTTAC")

groups <- lapply(tree$tip.label, f_label)

groups <- as.character(groups)

groupInfo <- split(tree$tip.label, groups)

tree <- root(tree,outgroup = outgroup, resolve.root = T)

tree <- groupOTU(tree, groupInfo)

p <- ggtree(tree, aes(color=group)) + 
  #geom_tiplab(hjust = 0.1, size=3, align = TRUE) +
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.3)
  geom_treescale() +
  labs(colour= "Population") +
  scale_color_manual(values = c("BRE" = "#FED976", "LE" = "#BCBDDC", 'Outgroup' ='pink')) +
  theme(legend.position="bottom", text = element_text(size =18))

p2 <- collapse(p, node = 44, 'max', fill = "#BCBDDC") 


ggsave('Plots/BRE-LE_tree.png', width = 8, height = 6)
