# Bottleneck simulation for Aim 1 - Sm Single GT project
# Author: Kathrin Jutzeler
# Date: February 15, 2024
# Last updated: February 16, 2024
# R version 4.2.0, tidyverse version 1.3.2

library(radiator)
library(vcfR)
library(tidyverse)


# Generate input file ####
vcf.in <- read.vcfR('maf05_CDS_autosome_br_v10.vcf')

#loci <- data.frame(vcfR::getID(vcf.BRE)[1:10])

gt <- data.frame(t(vcfR::extract.gt(vcf.in)))

dim <- as.data.frame(dim(gt))
colnames(dim) <- NULL
dim[2,1] <- 10000

write.table(dim, 'botsim_input_br.txt', quote = F, row.names = F)

gt[is.na(gt)] <- "? ?"

gt <- gt %>%
  mutate(across(everything(), ~str_replace(., "[/|]"," "))) 

rownames(gt) <- gsub('.hc.vcf.gz', '', rownames(gt))

rm <- c(13,16, 34)

gt10000 <- gt[-rm,sample(1:ncol(gt),10000)]

write.table(gt10000, 'botsim_input_br.txt', quote = F, append = T)

# Run Bottlesim.exe
#---------------------------------------------------------------------------#
## Percentage Ho ####
f_Ho <- function(x){
  start <- which(grepl("Start of Ho data", x))
  end <- which(grepl("End of Ho data", x))
  df <- x[start:end]
  df <- df[-1:-2]
  df <- head(df, -2)
  return(as.data.frame(df))
}

# Import output file(s) ####
#bottlesim <- list.files(pattern="bs_nmri")
bottlesim <- list.files(pattern="nooverlap")
bottlesim2 <- list.files(pattern="br_overlap")

bs_list <- lapply(bottlesim, readLines)
bs_list2 <- lapply(bottlesim2, readLines)

#t <- lapply(bs_list, readLines)

bs_df <- lapply(bs_list, f_Ho)
bs_df2 <- lapply(bs_list2, f_Ho)

names(bs_df) <- c('10', '100', '200', '25', '400', '5', '50')
names(bs_df2) <- c('10', '100', '200', '25', '400', '5', '50')

column_names <- c('blank', 'Year', as.character(seq(1:length(gt10000))), 'Average', 'SD', 'SE', 'Percentage')

#BRE_df <- as.data.frame(BRE_df)
t <- lapply(bs_df, function(x) separate(x, col = df, sep = '\\s+', into = column_names))

t <- lapply(t, function(x) x %>% mutate(across(everything(), as.numeric)))

bs_out <- bind_rows(t, .id = 'size')
bs_out$size <- as.factor(bs_out$size)

t2 <- lapply(bs_df2, function(x) separate(x, col = df, sep = '\\s+', into = column_names))

t2 <- lapply(t2, function(x) x %>% mutate(across(everything(), as.numeric)))

bs_out2 <- bind_rows(t2, .id = 'size')
bs_out2$size <- as.factor(bs_out2$size)

bs_out$feature <- 'no_overlap'
bs_out2$feature <- 'overlap'

p <- bind_rows(bs_out, bs_out2)

p$size <- as.numeric(p$size)

# Plot Ho ####
#p_Ho <-
ggplot(bs_out, aes(Year, Percentage, group = size, color = size)) +
  geom_line() +
  labs(x = 'Generations', y = '%' ~ italic(Ho)~ 'retained', color = "Population size") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, vjust = 0.7), axis.title = element_text(size = 14, face = 'bold'),
        plot.title = element_text(hjust = 0.5, size = 16),
        text = element_text(size = 12), axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  geom_hline(yintercept =49, alpha = 0.2, linetype = 'dashed')
  

ggplot(p, aes(Year, Percentage, group = interaction(size,feature), linetype =  feature)) +
  geom_line() +
  labs(x = 'Generations', y = '%' ~ italic(Ho)~ 'retained', linetype = "Maintenance") +
  scale_linetype_manual(values = c('solid', 'dashed'), labels = c('no_overlap' = 'No overlap', 'overlap' = 'Overlap')) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, vjust = 0.7), axis.title = element_text(size = 14, face = 'bold'),
        plot.title = element_text(hjust = 0.5, size = 16),
        text = element_text(size = 12), axis.line = element_line(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.justification = 'top')  +
  geom_hline(yintercept =49, alpha = 0.3, linetype = 'dashed') #+ 
  #ggrepel::geom_text_repel(aes(label = paste0('N = ', size)), subset(p, Year == 80 & size < 50)) +
  #ggrepel::geom_text_repel(aes(label = paste0('N = ', size)), subset(p, Year == 400 & size > 49),
             #              position = position_dodge(0.5)) 

ggsave('bottlesim.jpg', width = 8, heigh =6, background = 'white')


#ggsave('Bottlesim.png', unit = 'in', width = 8, height = 6)

# Get a specific number
df5 <- filter(bs_out2, size == "5")
df100 <- filter(bs_out2, size == "100")

range(df100[df100$Year == 400, "Percentage"])
range(df5[df5$Percentage >= 1, "Year"])

range(df5[df5$Year == 56, "Percentage"])

which(df5$Percentage >= 1, 'Year')

