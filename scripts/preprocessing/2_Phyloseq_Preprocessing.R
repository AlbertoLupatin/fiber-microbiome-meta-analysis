library(common)
library(tidyverse)
library(phyloseq)
library(Interventions)

Metadata_list <- readRDS("../../data/preprocessing/Metadata.RDS")

studies_16S <- read.csv("../../data/query/Query_Final.csv") %>% 
  filter(Notes == "16s") %>% 
  pull(First.Author) %>% 
  gsub("-", "", .)

Phyloseq_all <- list()

for (i in seq(length(Metadata_list))) {
  
  study_name <- names(Metadata_list)[i]
  cat(paste0(i, " - ", study_name, "\n"))
  
  meta <- Metadata_list[[i]]
  rownames(meta) <- meta$Run
  
  # Load abundance matrix
  study_name_rds <- gsub("^(([^_]+)_([^_]+))_.*", "\\1", study_name)
  
  path <- paste0("../../data/query/", study_name_rds, "/", study_name_rds, "_Taxonomic.rds") 
  study <- readRDS(path) %>%
    as.data.frame() 
  

  if (study_name %in% studies_16S) {
    study <- study[grep("g__", rownames(study)),]
  } else {
    study <- study[grep("ref_mOTU", rownames(study)),]
  }

  study <- study %>% 
    rownames_to_column("taxon") %>%
    pivot_longer(-taxon) %>%
    dplyr::rename(Run = name, Counts = value) %>% 
    mutate(Run = str_c(study_name, Run, sep = "_")) %>% 
    na.omit(.)
  
  combined <- study %>%
    group_by(taxon, Run) %>%
    mutate(Counts = sum(Counts)) %>%
    ungroup() %>%
    distinct() %>%
    pivot_wider(names_from = Run, values_from = Counts) %>%
    column_to_rownames("taxon") %>%
    filter(rownames(.) != "unassigned") %>%
    as.matrix()
  
  if (grepl("Korem", study_name)) {colnames(combined) <- gsub(".singles", "", colnames(combined))}
  combined <- combined[, colnames(combined) %in% rownames(meta)]
  
  # Create taxonomy table
  taxonomy <- as.data.frame(rownames(combined)) %>%
    rename("Full_taxonomy" =  1) %>%
    filter(Full_taxonomy != "unassigned") %>%
    separate_wider_delim(Full_taxonomy, "|", names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Specie", "ref_mOTU"), too_few = "align_start") %>%
    mutate(full_name = rownames(combined)) %>%
    column_to_rownames(var = "full_name") %>%
    as.matrix()
  
  
  # Create phyloseq object
  class(combined) <- "numeric"
  OTU <- otu_table(combined, taxa_are_rows = T)
  TAX <- tax_table(taxonomy)

  meta <- as.data.frame(meta)
  rownames(meta) <- meta$Run # For some reason I have to assign the row names again

  tax_ps <- phyloseq(OTU, TAX, sample_data(meta))
  tax_ps@otu_table <- tax_ps@otu_table[rownames(tax_ps@otu_table) != "unassigned", ]

  # Convert to absolute counts
  tax_ps_rel <- microbiome::transform(tax_ps, "compositional")

  Phyloseq_all[[study_name]] <- tax_ps_rel

  if (i == length(Metadata_list)) {
    cat("Operation Completed\nSaving and Exiting\n")
  }
  
}

all_genera_16S <- lapply(Phyloseq_all[studies_16S], function(ph) {
  ph@tax_table@.Data[, 6]
})
common_genera_16S <- Reduce(intersect, all_genera_16S)

Phyloseq_16S_filt <- lapply(Phyloseq_all[studies_16S], function(ph) {subset_taxa(ph, Genus %in% common_genera_16S)})
Phyloseq_all[studies_16S] <- Phyloseq_16S_filt

saveRDS(common_genera_16S, "../../data/preprocessing/Common Genera.rds")
saveRDS(Phyloseq_all, "../../data/preprocessing/Phyloseq.rds")
