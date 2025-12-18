library(tidyverse)
setwd("Rodriguez Curated Studies")

collated_profiling <- as.data.frame(readRDS("Collated_res_IDTaxa.rds"))
collated_metadata <- read.delim("Collated_Metadata.txt")


# DAHL
dahl_runs <- collated_metadata %>% 
  filter(grepl("dahl", sample_title)) %>% 
  select(run_accession, sample_title)

dahl_meta <- read.csv("DahlWJ_2015.txt", sep = "\t") %>% 
  filter(treatment == "fiber") %>% 
  mutate(Intervention = "RS4",
         Individual = as.numeric(gsub("dahl_", "", subject_id)),
         Timepoint = as.factor(ifelse(timepoint == "before", "Before", "After"))) %>% 
  select(sampleid, Individual, Timepoint, Intervention) %>% 
  
  merge(dahl_runs, by.x = "sampleid", by.y = "sample_title", all.x = T) %>% 
  rename(Run = run_accession) %>% 
  select(-sampleid)
rownames(dahl_meta) <- dahl_meta$Run
write.csv(dahl_meta, "../DahlWJ_2015/DahlWJ_2015_metadata.csv")

dahl_cluster_meta <- read.csv("DahlWJ_2015.txt", sep = "\t") %>% 
  merge(dahl_runs, by.x = "sampleid", by.y = "sample_title") %>% 
  rename(Run = run_accession) 
rownames(dahl_cluster_meta) <- dahl_cluster_meta$Run
write.csv(dahl_cluster_meta, "../DahlWJ_2015/DahlWJ_2015_SRA.csv")

dahl_profiling <- collated_profiling[, dahl_meta$Run]
saveRDS(dahl_profiling, "../DahlWJ_2015/DahlWJ_2015_Taxonomic.rds")


# MORALES
morales_runs <- collated_metadata %>% 
  filter(str_detect(sample_title, "morales")) %>% 
  select(run_accession, sample_title)

morales_meta <- read.csv("MoralesP_2016.txt", sep = "\t") %>% 
  filter(fiber_type == "oligofructose") %>% 
  mutate(Intervention = "Inulin",
         Individual = as.numeric(gsub("morales_", "", subject_id)),
         Timepoint = as.factor(ifelse(timepoint == "before", "Before", "After"))) %>% 
  select(sampleid, Individual, Timepoint, Intervention) %>% 

  merge(morales_runs, by.x = "sampleid", by.y = "sample_title", all.x = T) %>% 
  rename(Run = run_accession) %>% 
  select(-sampleid)
rownames(morales_meta) <- morales_meta$Run
write.csv(morales_meta, "../MoralesP_2016/MoralesP_2016_metadata.csv")

morales_cluster_meta <- read.csv("MoralesP_2016.txt", sep = "\t") %>% 
  merge(morales_runs, by.x = "sampleid", by.y = "sample_title") %>% 
  rename(Run = run_accession) 
rownames(morales_cluster_meta) <- morales_cluster_meta$Run
write.csv(morales_cluster_meta, "../MoralesP_2016/MoralesP_2016_SRA.csv")

morales_profiling <- collated_profiling[, morales_meta$Run]
saveRDS(morales_profiling, "../MoralesP_2016/MoralesP_2016_Taxonomic.rds")