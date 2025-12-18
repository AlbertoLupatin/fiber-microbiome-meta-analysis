library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggembl)

metadata <- readRDS("../data/preprocessing/Metadata.RDS")

shotgun_studies <- read.csv("../data/query/Query_Final.csv") %>% 
  filter(Notes == "Shotgun") %>% 
  mutate(Study = gsub("-", "", First.Author)) %>% 
  pull(Study)

sts_studies <- read.csv("../data/query/Query_Final.csv") %>% 
  filter(Notes == "16s") %>% 
  mutate(Study = gsub("-", "", First.Author)) %>% 
  pull(Study)

readcounts <- NULL

for (study_name in shotgun_studies) {
  
  print(study_name)
  meta <- metadata[[study_name]]
  
  if (study_name == "RehnerJ_2023") {
      samples_toselect <- gsub(".*2023_", "", meta$Run)
    } else if (study_name == "WastykHC_2021") {
      samples_toselect <- gsub(".*2021_", "", meta$Run)
    } else {
      samples_toselect <- gsub(".*_", "", meta$Run)
    }
  
  
  study_name_path <- gsub("(.*\\d{4}).*$", "\\1", study_name)
  
  tax <- readRDS(paste0("../data/query/", study_name_path, "/", study_name_path, "_Taxonomic.rds"))
  
  if (study_name == "KoremT_2017") {colnames(tax) <- gsub(".singles", "", colnames(tax))}
  
  tax <- tax[, samples_toselect]
  
  readcounts_tmp <- data.frame(colSums(tax))
  readcounts_tmp$Study <- study_name
  readcounts_tmp$SeqMethod <- "Shotgun"
  
  if (is.null(readcounts) == TRUE) {
    readcounts <- readcounts_tmp
  } else {
    readcounts <- bind_rows(readcounts, readcounts_tmp)
    }

}

for (study_name in sts_studies) {
  
  if (study_name == "NishimotoY_2023") next
  if (study_name == "InoueR_2025") next
  print(study_name)
  meta <- metadata[[study_name]]
  
  samples_toselect <- gsub(".*_", "", meta$Run)
  
  study_name_path <- gsub("(.*\\d{4}).*$", "\\1", study_name)
  
  tax <- readRDS(paste0("../data/query/", study_name_path, "/", study_name_path, "_Taxonomic.rds"))
  
  if (study_name == "NishimotoY_2023") {colnames(tax) <- gsub(".singles", "", colnames(tax))}
  
  tax <- tax[, samples_toselect]
  
  readcounts_tmp <- data.frame(colSums(tax))
  readcounts_tmp$Study <- study_name
  readcounts_tmp$SeqMethod <- "16S"

  readcounts <- bind_rows(readcounts, readcounts_tmp)
  
}

ggplot(readcounts, aes(x = SeqMethod, y = colSums.tax.)) +
  geom_boxplot() +
  stat_compare_means(method = "t.test", label = "p.format", label.x = 1.4) +
  ylab("Read Count") +
  xlab("Sequencing Method") +
  theme_publication() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) 
ggsave("../figures/Read Counts.pdf")
