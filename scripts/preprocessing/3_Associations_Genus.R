options(warn=-1)

library(tidyverse)
library(Interventions)
library(SIAMCAT)

Metadata_list <- readRDS("../../data/preprocessing/Metadata.RDS")
Phyloseq_all <- readRDS("../../data/preprocessing/Phyloseq.RDS")

shotgun_studies <- read.csv("../../data/query/Query_Final.csv") %>% 
  filter(Notes == "Shotgun") %>% 
  pull(First.Author) %>% 
  gsub("-", "", .)

Phyloseq_list_genera <- list()
Associations_list <- list()

for (i in seq(length(Metadata_list))) {

  study_name <- names(Metadata_list)[i]
  
  if (study_name %in% shotgun_studies) next
  cat(paste0(i, " - ", study_name, "\n"))

  meta <- data.frame(Metadata_list[[i]])
  rownames(meta) <- meta$Run
  tax_ps_rel <- Phyloseq_all[[i]]
  
  
  OTU_genera <- otu_table(tax_ps_rel) %>% 
    as.data.frame() %>% 
    filter(grepl("g__", rownames(.))) 
  rownames(OTU_genera) <- gsub(".*g__", "", rownames(OTU_genera))
  OTU_genera <- as.matrix(OTU_genera)
  OTU_genera <- otu_table(OTU_genera, taxa_are_rows = T)
  
  tax_ps_rel_genera <- phyloseq(OTU_genera, sample_data(meta))
  Phyloseq_list_genera[[study_name]] <- tax_ps_rel_genera
  
  # SIAMcat object + Associations
  label <- create.label(meta = sample_data(tax_ps_rel_genera), label = "Timepoint", case = "After", verbose = 0)
  
  sc <- siamcat(phyloseq = tax_ps_rel_genera, label = label, verbose = 0)
  sc <- filter.features(sc, filter.method = "prevalence", cutoff = 0.1, verbose = 0)
  sc <- normalize.features(sc, norm.method = "log.std", verbose = 0)
  
  sc <- check.associations(sc, formula = "feat~label+(1|Individual)", verbose = 0)
  Associations_list[[study_name]] <- associations(sc)
  
  if (i == length(Metadata_list)) {
    cat("Associations Separated - Operation Completed\n\n")
  }
  
}

saveRDS(Phyloseq_list_genera, "../../data/preprocessing/Phyloseq_Genus_Separated.rds")
saveRDS(Associations_list, "../../data/preprocessing/Associations_Genus_Separated.rds")

# Merge metadata dataframes and phyloseq objects
Metadata_list_toMerge <- list()
for (study_name in names(Metadata_list)) {
  Metadata_list_toMerge[[study_name]] <- Metadata_list[[study_name]] %>% 
    mutate(Study = study_name,
           Individual = str_c(Study, Individual, sep = "_"))
}

Phyloseq_merged <- do.call(merge_phyloseq, Phyloseq_list_genera)

Metadata_merged <- purrr::reduce(.x = Metadata_list_toMerge,
                          merge, all = T) %>% 
  column_to_rownames("Run") %>% 
  filter(SeqMethod == "16S") %>% 
  select(-X)

# Compute associations
label_merged <- create.label(meta = Metadata_merged, label = "Timepoint", case = "After")

SC_merged <- siamcat(feat = Phyloseq_merged@otu_table, meta = Metadata_merged, label = label_merged, verbose = 0)

SC_merged <- filter.features(SC_merged, filter.method = "prevalence", cutoff = 0.1, verbose = 0)
SC_merged <- normalize.features(SC_merged, norm.method = "log.std")

SC_merged <- check.associations(SC_merged, formula = "feat~label+(1|Study)+(1|Individual)", verbose = 0)
Associations_merged <- associations(SC_merged)

cat("Associations Merged - Operation Completed\nSaving and Exiting\n")

saveRDS(Associations_merged, "../../data/preprocessing/Associations_Genus_Combined.rds")

