library(SNPRelate)
library(gdsfmt)
library(tidyverse)

# Define colors
population_colors = c("OR" = "#BF812D", "BRE" = "#FED976", "LE" = "#BCBDDC", "EG" = "#D6604D")
field_colors <- c('Brazil' = "olivedrab3", 'Senegal' = "red", 'Niger'="brown4" ,'Tanzania'='gray37',
                  'Caribbean' = "lightskyblue1", 'Egypt' = "steelblue3")

# Import the data
vcf.in <- "C:/Users/kbailey/Downloads/temp/final_lab_nuclear_DNA.vcf"

gds <- snpgdsVCF2GDS(vcf.in, "final_nuclear_DNA.gds", method="biallelic.only")

genofile <- snpgdsOpen(gds)
#snpgdsClose(genofile)

pca <- snpgdsPCA(genofile, num.thread=2, autosome.only = F)

samples <- as.data.frame(pca$sample.id)

colnames(samples) <- "name"
pca_samples <- separate(samples, name, c("sample", "ID", "lane"))

OR <- c("5", "1", "6", "2", "24", "7", "25", "14", "26", "17")
BRE <- c("27", "32", "51", "38", "68", "48", "72", "49", "80", "55")
LE <- c("84", "92", "86", "104", "90", "107", "97", "109", "112", "111")
EG <- c("138", "130", "154", "135", "161", "142"," 167", "146", "169", "163")        

metadata <- pca_samples %>%
  mutate(group = ifelse(sample %in% OR, "OR",
                        ifelse(sample %in% BRE, "BRE", 
                               ifelse(sample %in% LE, "LE", "EG"))))
                               
data <-
  data.frame(sample.id = pca$sample.id,
             EV1 = pca$eigenvect[,1],
             EV2 = pca$eigenvect[,2],
             EV3 = pca$eigenvect[,3],
             EV4 = pca$eigenvect[,4],
             EV5 = pca$eigenvect[,5],
             EV6 = pca$eigenvect[,6],
             group = metadata$group,
             stringsAsFactors = FALSE)

#p_pca_nuclear_DNA <-
  ggplot(subset(data, sample.id != "5_S29_L003"), aes(EV1, EV2, col = group)) +
  geom_point(size = 3) +
  #geom_text(size=4) +
  theme_bw() +
  labs(x = paste0("PC1 variance: ",round(pca$varprop[1]*100,digits=2),"%"),
     y = paste0("PC2 variance: ",round(pca$varprop[2]*100,digits=2),"%"),
      color = "Population") +
  scale_color_manual(values = population_colors) +
  ggtitle(expression(atop(paste(italic("S. mansoni"), " Single Genotypes")))) +
  theme(text = element_text(size = 18))


ggsave("Plots/PCA.png", width = 8, height = 6)


### Plot with field samples
vcf.in <- "C:/Users/kbailey/Downloads/temp/labvfield_nuclear_annotated.vcf"
vcf.in <- "C:/Users/kbailey/Downloads/temp/labvfield_norm.vcf"
metadata <- read.csv('samples_pca.csv')

gds <- snpgdsVCF2GDS(vcf.in, "labvfield.gds", method="biallelic.only")

genofile <- snpgdsOpen(gds)
#snpgdsClose(genofile)

pca <- snpgdsPCA(genofile, num.thread=2, autosome.only = F)

samples <- as.data.frame(pca$sample.id)

colnames(samples) <- "name"

metadata <- metadata %>%
  full_join(samples, by = c('sample' = 'name'))

data <-
  data.frame(sample.id = pca$sample.id,
             EV1 = pca$eigenvect[,1],
             EV2 = pca$eigenvect[,2],
             EV3 = pca$eigenvect[,3],
             EV4 = pca$eigenvect[,4],
             EV5 = pca$eigenvect[,5],
             EV6 = pca$eigenvect[,6],
             group = metadata$population,
             origin = metadata$origin,
             stringsAsFactors = FALSE)

data <- data %>%
  filter(!sample.id %in% c('5_S29_L003', "sm_guadeloupe_3", "sm_kenya", "sm_pr", "sm_uganda","sm_guadeloupe_1" ))

#p_pca_nuclear_DNA <-
ggplot(data, aes(EV1, EV2, col = group)) +
  geom_point(size = 3) +
  #geom_text(size=4) +
  theme_bw() +
  labs(x = paste0("PC1 variance: ",round(pca$varprop[1]*100,digits=2),"%"),
       y = paste0("PC2 variance: ",round(pca$varprop[2]*100,digits=2),"%"),
       color = "Population", shape = "Sample origin") +
  scale_color_manual(values = all_colors) +
  ggtitle(expression(atop(paste(italic("S. mansoni"), " Single Genotypes")))) +
  theme(text = element_text(size =18)) +
  stat_ellipse(aes(group = origin), type ='norm')

ggsave("Plots/pca_field.png", width = 8, height = 6)

set.seed(1231)

result2 <- sample(palette,8)

pie(rep(1, 8), col = result)
pie(rep(1, 9), col = result2)
