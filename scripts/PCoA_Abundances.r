library(tidyverse)
library(Interventions)
library(vegan)
library(ggplot2)
library(ggembl)

Metadata_list <- readRDS("../data/preprocessing/Metadata.RDS")
Metadata_list_toMerge <- list()
for (study_name in names(Metadata_list)) {
  Metadata_list_toMerge[[study_name]] <- Metadata_list[[study_name]] %>% 
    mutate(Study = study_name,
           Individual = str_c(Study, Individual, sep = "_"))
}

Metadata_merged <- purrr::reduce(.x = Metadata_list_toMerge,
                                 merge, all = T) %>% 
  column_to_rownames("Run") %>% 
  select(-X)

shotgun_studies <- read.csv("../data/query/Query_Final.csv") %>% 
  filter(Notes == "Shotgun") %>% 
  pull(First.Author) %>% 
  gsub("-", "", .)

common_genera <- readRDS("../data/preprocessing/Common Genera.rds")

abundances <- NULL

for (name in names(Metadata_list)) {
  
  if (grepl("Nishimoto", name)) next
  cat(name, "\n\n")
  
  meta_tmp <- Metadata_list[[name]]

  name_path <- gsub("(\\d{4})\\d*.*$", "\\1", name)
  path <- paste0("../data/query/", name_path, "/", name_path, "_Taxonomic.rds")
  
  taxonomy_tmp <- as.data.frame(readRDS(path)) 
  colnames(taxonomy_tmp) <- paste0(name, "_", colnames(taxonomy_tmp))
  colnames(taxonomy_tmp) <- gsub(".singles", "", colnames(taxonomy_tmp))
  if (name == "InoueR_2025") {meta_tmp <- meta_tmp %>% filter(!(Individual %in% meta_tmp$Individual[which(meta_tmp$Run %in% setdiff(meta_tmp$Run, colnames(taxonomy_tmp)))]))}
  taxonomy_tmp <- taxonomy_tmp[, meta_tmp$Run]

  if (name %in% shotgun_studies) {
    
    taxonomy_tmp <- taxonomy_tmp %>%
      rownames_to_column("Taxon") %>%
      pivot_longer(-Taxon, names_to = "Run", values_to = "Abundance") %>%
      mutate(Specie = clearNames(Taxon),
             Genus = gsub("\\s.*", "", Specie),
             Genus = paste0("g__", Genus)) %>%
      filter(Genus %in% common_genera) %>%
      group_by(Run, Genus) %>%
      summarize(Abundance = sum(Abundance)) %>%
      ungroup() %>%
      pivot_wider(names_from = Run, values_from = Abundance) %>%
      column_to_rownames("Genus")
    taxonomy_tmp <- taxonomy_tmp[rowSums(taxonomy_tmp[]) > 0,]
    
  } else {
    
    # taxonomy_tmp <- taxonomy_tmp[grep("g__", rownames(taxonomy_tmp)),] 
    rownames(taxonomy_tmp) <- gsub(".*g__", "g__", rownames(taxonomy_tmp))
    taxonomy_tmp <- taxonomy_tmp[which(rownames(taxonomy_tmp) %in% common_genera),]
    
  }
  
  # Absolute -> Relative counts (equivalent to microbiome::transform - compositional)
  taxonomy_rel_tmp <- apply(taxonomy_tmp, 2, function(x) {x / max(sum(x, 1e-32))})
  taxonomy_rel_tmp <- as.data.frame(t(taxonomy_rel_tmp))
  
  if (is.null(abundances)) {
    abundances <- taxonomy_rel_tmp
  } else {
    abundances <- bind_rows(abundances, taxonomy_rel_tmp)
  }
}


# RETRIEVE COMMON GENERA
abundances <- abundances[, !colSums(is.na(abundances)), drop = F]


# COMPUTE DISTANCES
bc <- vegdist(abundances, method = "bray")
pcoa <- cmdscale(bc, eig = T, add = T)

position <- pcoa$points
colnames(position) <- c("PCOA1", "PCOA2")
position <- as_tibble(position, rownames = "samples") %>%
  merge(Metadata_merged, by.x = "samples", by.y = 0) %>%
  mutate(SequencingType = ifelse(Study %in% shotgun_studies, "Shotgun", "16S"))


explained_variance <- round(100 * pcoa$eig / sum(pcoa$eig), digits = 2)
# PLOT PCoA
ggplot(position, aes(x = PCOA1, y = PCOA2, color = SequencingType)) +
  geom_point() +
  scale_color_manual(values = c("Shotgun" = "black", "16S" = "grey")) +
  xlab(paste0("PCoA 1 (", explained_variance[1], "%)")) +
  ylab(paste0("PCoA 2 (", explained_variance[2], "%)")) +
  theme_publication() +
  theme(legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        aspect.ratio = 1)

ggsave("../figures/PCoA Seq Method.pdf", width = 7, height = 7)
