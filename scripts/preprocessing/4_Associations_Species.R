library(tidyverse)
library(Interventions)
library(SIAMCAT)

shotgun_studies <- read.csv("../../data/query/Query_Final.csv") %>% 
  filter(Notes == "Shotgun") %>% 
  pull(First.Author) %>% 
  gsub("-", "", .)

Metadata_list <- readRDS("../../data/preprocessing/Metadata.RDS")[shotgun_studies]
Phyloseq_list <- readRDS("../../data/preprocessing/Phyloseq.RDS")[shotgun_studies]


# SEPARATED
Associations_separated <- list()

for (study_name in names(Phyloseq_list)) {
  cat(study_name, "\n")
  
  phylo <- Phyloseq_list[[study_name]]
  
  label <- create.label(meta = sample_data(phylo), label = "Timepoint", case = "After", verbose = 0)
  
  sc <- siamcat(phyloseq = phylo, label = label, verbose = 0)
  sc <- filter.features(sc, filter.method = "prevalence", cutoff = 0.1, verbose = 0)
  sc <- normalize.features(sc, norm.method = "log.std", verbose = 0)
  
  sc <- check.associations(sc, formula = "feat~label+(1|Individual)", verbose = 0)
  Associations_separated[[study_name]] <- associations(sc)
}

cat("Associations Separated - Operation Completed\n")
saveRDS(Associations_separated, "../../data/preprocessing/Associations_Species_Separated.rds")


# COMBINED
for(study_name in names(Metadata_list)) {
  
  Metadata_list[[study_name]] <- Metadata_list[[study_name]] %>%
    mutate(Individual = paste0(study_name, as.character(Individual)),
           Study = study_name) 
}

metadata_merged <- purrr::reduce(.x = Metadata_list, merge, all = T) %>% 
  column_to_rownames("Run") %>%
  filter(Study %in% shotgun_studies)

phyloseq_merged <- do.call(merge_phyloseq, Phyloseq_list)


label_merged <- create.label(meta = metadata_merged, label = "Timepoint", case = "After")

sc_merged <- siamcat(feat = phyloseq_merged@otu_table, meta = metadata_merged, label = label_merged, verbose = 0)

sc_merged <- filter.features(sc_merged, filter.method = "prevalence", cutoff = 0.1, verbose = 0)
sc_merged <- normalize.features(sc_merged, norm.method = "log.std")

sc_merged <- check.associations(sc_merged, formula = "feat~label+(1|Study)+(1|Individual)", verbose = 0)

cat("Associations Combined - Operation Completed\nSaving and Exiting\n")
saveRDS(associations(sc_merged), "../../data/preprocessing/Associations_Species_Combined.rds")


