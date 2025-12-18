library(tidyverse)
library(Interventions)

shotgun_studies <- read.csv("../../data/query/Query_final.csv") %>% 
  filter(Notes == "Shotgun" & !(grepl("Omar", First.Author))) %>% 
  mutate(First.Author = gsub("-", "", First.Author)) %>% 
  pull(First.Author)
author_names <- shotgun_studies %>% gsub("^(.*?\\d{4}).*$", "\\1", .) %>% unique()

Meta_all <- readRDS("../../data/preprocessing/Metadata.RDS")[shotgun_studies]


for(study_name in names(Meta_all)) {
  Meta_all[[study_name]] <- Meta_all[[study_name]] %>%
    mutate(Individual = paste0(study_name, "_", as.character(Individual)),
           Study = study_name)
}
all_meta <- merge_metadata(Meta_all)
all_meta$Run <- gsub(".*_20[0-9]{2}_", "", rownames(all_meta))
all_meta$Run <- gsub("Pea_|PeaFructan_|Extruded_|Complete_|Orange_", "", all_meta$Run)


Feat_all <- list()

for (study in author_names) {
  
  cat(study, "\n")
  
  tmp_feat <- read.delim(file.path("../../data/query", study, paste0(study, "_CAZy.txt"))) %>% 
    filter(!(row_number() %in% c(1, 2, 3))) %>% 
    column_to_rownames("feature")
  
  if (study == "DelannoyBrunoO_2022") {
    tmp_feat_pea <- tmp_feat[, all_meta %>% 
                               filter(Study == "DelannoyBrunoO_2022_Pea") %>%
                               pull(Run)]
    colnames(tmp_feat_pea) <- paste0("DelannoyBrunoO_2022_Pea_", colnames(tmp_feat_pea))
    Feat_all[[paste0(study, "_Pea")]] <- tmp_feat_pea
    
    tmp_feat_orange <- tmp_feat[, all_meta %>% 
                                  filter(Study == "DelannoyBrunoO_2022_Orange") %>% 
                                  pull(Run)]
    colnames(tmp_feat_orange) <- paste0("DelannoyBrunoO_2022_Orange_", colnames(tmp_feat_orange))
    Feat_all[[paste0(study, "_Orange")]] <- tmp_feat_orange
    
  } else if (study == "DelannoyBrunoO_2021") {
    
    extra_sample <- read.delim(file.path("../../Data/Query", study, "ERR5194045_CAZy.txt")) %>% 
      filter(!(row_number() %in% c(1, 2, 3))) %>%
      filter(feature %in% rownames(tmp_feat)) %>% 
      column_to_rownames("feature") %>% 
      dplyr::rename(ERR5194045 = combined_rpkm) %>% 
      select(ERR5194045)
    
    extra_sample[setdiff(rownames(tmp_feat), rownames(extra_sample)), ] <- NA
    
    tmp_feat <- bind_cols(tmp_feat, extra_sample)
    
    tmp_feat_extruded <- tmp_feat[, all_meta %>% 
                                    filter(Study == "DelannoyBrunoO_2021_Extruded") %>%
                                    filter(!grepl("DelannoyBrunoO_2021_Extruded_14|09", Individual)) %>% 
                                    pull(Run)]
    colnames(tmp_feat_extruded) <- paste0("DelannoyBrunoO_2021_Extruded_", colnames(tmp_feat_extruded))
    Feat_all[[paste0(study, "_Extruded")]] <- tmp_feat_extruded
    
    tmp_feat_complete <- tmp_feat[, all_meta %>% 
                                    filter(Study == "DelannoyBrunoO_2021_Complete") %>%
                                    pull(Run)]
    colnames(tmp_feat_complete) <- paste0("DelannoyBrunoO_2021_Complete_", colnames(tmp_feat_complete))
    Feat_all[[paste0(study, "_Complete")]] <- tmp_feat_complete
    
    tmp_feat_peafructan <- tmp_feat[, all_meta %>% 
                                      filter(Study == "DelannoyBrunoO_2021_PeaFructan") %>%
                                      pull(Run)]
    colnames(tmp_feat_peafructan) <- paste0("DelannoyBrunoO_2021_PeaFructan_", colnames(tmp_feat_peafructan))
    Feat_all[[paste0(study, "_PeaFructan")]] <- tmp_feat_peafructan
  }
  
  else {
    tmp_feat <- tmp_feat[, all_meta %>% 
                           filter(Study == study) %>% 
                           pull(Run)]
    colnames(tmp_feat) <- paste0(study, "_", colnames(tmp_feat))
    Feat_all[[study]] <- tmp_feat
  }
}

cat("\nOperation Completed\nSaving and Exiting\n\n")
saveRDS(Feat_all, "../../data/preprocessing/CAZymes_List.RDS")