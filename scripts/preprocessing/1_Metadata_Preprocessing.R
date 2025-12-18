options(warn=-1)

library(openxlsx)
library(tidyverse)

setwd("../../data/query/")

Meta_all <- list()


# Barber
cat("BarberC_2021\n")
Barber_meta <- read.csv("BarberC_2021/BarberC_2021_metadata.tsv", sep = "\t") %>%
  separate(sample_title, sep = " ", into = c("tmp", "Individual", "Timepoint")) %>%
  mutate(Timepoint = ifelse(Timepoint == "WD", "Before", "After"),
         Run = str_c("BarberC_2021", run_accession, sep = "_"),
         Individual = as.character(Individual)) %>%
  select(Run, Individual, Timepoint) %>% 
  mutate(Intervention = "General Fibers") %>%
  group_by(Individual) %>%
  filter(n() == 2)
rownames(Barber_meta) <- Barber_meta$Run

Meta_all[["BarberC_2021"]] <- Barber_meta


# Baxter
cat("BaxterNT_2019\n")
Baxter_meta <- read.delim("BaxterNT_2019/BaxterNT_2019_metadata.txt") %>% 
  group_by(subject_id, timepoint) %>% 
  slice_max(timepoint_numeric, n = 1) %>% 
  slice(1) %>% 
  ungroup() %>% 
  mutate(Timepoint = ifelse(timepoint == "before", "Before", "After"),
         Intervention = case_when(fiber_type %in% c("potato", "himaize") ~ "RS2",
                                 fiber_type == "inulin" ~ "Fructans",
                                 fiber_type == "starch_control" ~ "Control"),
         Individual = as.character(gsub("baxter_U", "", subject_id))) %>% 
  select(sampleid, Individual, Timepoint, Intervention)

# RS2
Baxter_meta_RS2 <- Baxter_meta %>% 
  filter(Intervention == "RS2") %>% 
  mutate(Run = gsub("baxter_", "BaxterNT_2019_RS2_", sampleid)) %>% 
  select(-sampleid)

rownames(Baxter_meta_RS2) <- Baxter_meta_RS2$Run
Meta_all[["BaxterNT_2019_RS2"]] <- Baxter_meta_RS2
  
# Fructans
Baxter_meta_inulin <- Baxter_meta %>% 
  filter(Intervention == "Fructans") %>% 
  mutate(Run = gsub("baxter_", "BaxterNT_2019_Fructans_", sampleid)) %>% 
  select(-sampleid)

rownames(Baxter_meta_inulin) <- Baxter_meta_inulin$Run
Meta_all[["BaxterNT_2019_Fructans"]] <- Baxter_meta_inulin


# BenPaez
cat("BenitezPaezA_2019\n")
BenitezPaezA_2019_meta <- read.csv("BenitezPaezA_2019/BenitezPaezA_2019_metadata.tsv", sep = "\t") %>%
  separate(sample_alias, sep = "_", into = c("tmp", "Individual", "Timepoint")) %>%
  mutate(Timepoint = ifelse(Timepoint == "01", "Before", "After"),
         Intervention = "Arabinoxylan",
         Run = str_c("BenitezPaezA_2019", run_accession, sep = "_"),
         Individual = as.character(Individual)) %>%
  group_by(Individual, Timepoint) %>%
  slice_head(n = 2) %>%
  ungroup() %>% 
  select(Run, Individual, Timepoint, Intervention)
rownames(BenitezPaezA_2019_meta) <- BenitezPaezA_2019_meta$Run

Meta_all[["BenitezPaezA_2019"]] <- BenitezPaezA_2019_meta


# DahlWK_2015
cat("DahlWJ_2015\n")
DahlWJ_2015_meta <- read.csv("DahlWJ_2015/DahlWJ_2015_metadata.csv") %>% 
  mutate(Run = str_c("DahlWJ_2015", Run, sep = "_"),
         Individual = as.character(Individual))
rownames(DahlWJ_2015_meta) <- DahlWJ_2015_meta$Run

Meta_all[["DahlWJ_2015"]] <- DahlWJ_2015_meta


# Deehan_2020
cat("DeehanEC_2020\n")
Deehan_2020_meta <- read.delim("DeehanEC_2020/DeehanEC_2020_metadata.txt") %>% 
  filter(timepoint_id %in% c("before_1", "after_5") & treatment != "control") %>% 
  mutate(Run = gsub("deehan_", "DeehanEC_2020_", sampleid),
         Individual = str_extract(sample.name.original, "(?<=^[MF])\\d+(?=W)"),
         Timepoint = ifelse(timepoint == "before", "Before", "After"),
         Intervention = "RS4") %>% 
  select(Run, Individual, Timepoint, Intervention)
rownames(Deehan_2020_meta) <- Deehan_2020_meta$Run

Meta_all[["DeehanEC_2020"]] <- Deehan_2020_meta


# Deehan_2024
cat("DeehanEC_2024\n")
Deehan_2024_meta <- read.csv("DeehanEC_2024/DeehanEC_2024_metadata.csv") %>% 
  filter(treatment == "Arabinoxylan" & timepoint != "Week1") %>% 
  mutate(Individual = str_extract(source_material_id, "(?<=^[0-9][MF])\\d+"),
         Timepoint = ifelse(timepoint == "Baseline", "Before", "After"),
         Intervention = "Arabinoxylan",
         Run = str_c("DeehanEC_2024", Run, sep = "_")) %>% 
  select(Run, Individual, Timepoint, Intervention)
rownames(Deehan_2024_meta) <- Deehan_2024_meta$Run

Meta_all[["DeehanEC_2024"]] <- Deehan_2024_meta


# DelBruno 2021
cat("DelannoyBrunoO_2021\n")
DelannoyBrunoO_2021_meta <- read.csv("DelannoyBrunoO_2021/DelannoyBrunoO_2021_metadata.csv") %>%
  filter(host_scientific_name == "Homo sapiens" & LibrarySelection == "RANDOM")


# Extruded Pea Fiber
DelannoyBrunoO_2021_Extruded_meta <- DelannoyBrunoO_2021_meta %>% 
  filter(startsWith(Submitter_Id, "P")) %>% 
  separate(Submitter_Id, sep = "_", into = c("tmp", "Individual", "Timepoint")) %>% 
  mutate(Intervention = "Arabinoxylan + Pectins",
         Individual = as.character(gsub("S." ,"", Individual)),
         Timepoint = as.integer(gsub("Day", "", Timepoint))) %>% 
  filter(Diet %in% c("HiSF-LoFV (pre-intervention)", "Pea fiber (3 snacks a day)")) %>% 
  group_by(Individual, Diet) %>% 
  slice_max(Timepoint, n = 1) %>% 
  ungroup() %>% 
  mutate(Timepoint = ifelse(Timepoint < 16, "Before", "After"),
         Run = str_c("DelannoyBrunoO_2021_Extruded", Run, sep = "_")) %>% 
  select(Run, Individual, Timepoint, Intervention)

rownames(DelannoyBrunoO_2021_Extruded_meta) <- DelannoyBrunoO_2021_Extruded_meta$Run
Meta_all[["DelannoyBrunoO_2021_Extruded"]] <- DelannoyBrunoO_2021_Extruded_meta


# Pea fiber + Fructans
DelannoyBrunoO_2021_PeaFructan_meta <- DelannoyBrunoO_2021_meta %>% 
  filter(startsWith(Submitter_Id, "H")) %>% 
  separate(Submitter_Id, sep = "_", into = c("tmp1", "tmp2", "Individual", "Timepoint")) %>% 
  mutate(Intervention = "Arabinoxylan + Pectins + Fructans",
         Individual = as.character(gsub("S." ,"", Individual)),
         Timepoint = as.numeric(gsub("Day", "", Timepoint))) %>% 
  filter(Timepoint > 1 & Timepoint <= 25) %>% 
  group_by(Individual) %>% 
  arrange(Timepoint) %>% 
  filter(row_number() == 1 | row_number() == n()) %>% 
  mutate(Timepoint = ifelse(Timepoint < 13, "Before", "After"),
         Run = str_c("DelannoyBrunoO_2021_PeaFructan", Run, sep = "_")) %>% 
  select(Run, Individual, Timepoint, Intervention)

rownames(DelannoyBrunoO_2021_PeaFructan_meta) <- DelannoyBrunoO_2021_PeaFructan_meta$Run
Meta_all[["DelannoyBrunoO_2021_PeaFructan"]] <- DelannoyBrunoO_2021_PeaFructan_meta


# Complete Snack
DelannoyBrunoO_2021_Complete_meta <-  DelannoyBrunoO_2021_meta %>% 
  filter(startsWith(Submitter_Id, "H")) %>% 
  separate(Submitter_Id, sep = "_", into = c("tmp1", "tmp2", "Individual", "Timepoint")) %>% 
  mutate(Intervention = "Arabinoxylan + Pectins + Fructans + β-glucan",
         Individual = as.character(gsub("S." ,"", Individual)),
         Timepoint = as.numeric(gsub("Day", "", Timepoint))) %>% 
  filter(Timepoint > 1 & Timepoint > 27) %>% 
  group_by(Individual) %>% 
  arrange(Timepoint) %>% 
  filter(row_number() == 1 | row_number() == n()) %>% 
  mutate(Timepoint = ifelse(Timepoint < 40, "Before", "After"),
         Run = str_c("DelannoyBrunoO_2021_Complete", Run, sep = "_")) %>% 
  select(Run, Individual, Timepoint, Intervention)

rownames(DelannoyBrunoO_2021_Complete_meta) <- DelannoyBrunoO_2021_Complete_meta$Run
Meta_all[["DelannoyBrunoO_2021_Complete"]] <- DelannoyBrunoO_2021_Complete_meta


# DelBruno 2022
cat("DelannoyBrunoO_2022\n")
library(wakefield)
set.seed(100)
DelBruno2022_meta <- read.csv("DelannoyBrunoO_2022/DelannoyBrunoO_2022_metadata.csv") %>%
  filter(Assay.Type == "OTHER") %>%
  separate(Submitter_Id, sep = "_", into = c("Individual", "Timepoint", "Intervention")) %>%
  group_by(Individual, Intervention) %>%
  mutate(Individual = as.character(sample(1:10000, 1, replace = TRUE)))

## Pea
DelannoyBrunoO_2022_Pea_meta <- DelBruno2022_meta %>%
  filter(Intervention == "Study1PeaFiber" & (Timepoint == "Week1" | Timepoint == "Week9")) %>%
  mutate(Timepoint = ifelse(Timepoint == "Week1", "Before", "After"),
         Intervention = "Arabinoxylan + Pectins",
         Run = str_c("DelannoyBrunoO_2022_Pea", Run, sep = "_")) %>%
  select(Run, Individual, Timepoint, Intervention)
rownames(DelannoyBrunoO_2022_Pea_meta) <- DelannoyBrunoO_2022_Pea_meta$Run

Meta_all[["DelannoyBrunoO_2022_Pea"]] <- DelannoyBrunoO_2022_Pea_meta

## Orange
DelannoyBrunoO_2022_Orange_meta <- DelBruno2022_meta %>%
  filter(Intervention == "Study2OrangeFiber" & (Timepoint == "Week1" | Timepoint == "Week8")) %>%
  mutate(Timepoint = ifelse(Timepoint == "Week1", "Before", "After"),
         Intervention = "Arabinoxylan + Pectins",
         Run = str_c("DelannoyBrunoO_2022_Orange", Run, sep = "_")) %>%
  select(Run, Individual, Timepoint, Intervention)
rownames(DelannoyBrunoO_2022_Orange_meta) <- DelannoyBrunoO_2022_Orange_meta$Run

Meta_all[["DelannoyBrunoO_2022_Orange"]] <- DelannoyBrunoO_2022_Orange_meta


# DeMartino
DeMartino_meta <- read.csv("DeMartinoP_2022/DeMartinoP_2022_intervention_metadata.csv") %>% 
  merge(read.csv("DeMartinoP_2022/DeMartinoP_2022_sample_metadata.csv"), by.x = "Group", by.y = "Sample.Name") %>% 
  filter((Diet.x == "BL" & Week == "wk1") | (Diet.x == "Pot" & Week == "wk4")) %>% 
  mutate(Individual = as.character(gsub("POT ", "", StudyID.x)),
         Timepoint = ifelse(Week == "wk1", "Before", "After"),
         Intervention = "RS2",
         Run = str_c("DeMartinoP_2022", Run, sep = "_")) %>% 
  select(Run, Individual, Timepoint, Intervention)
rownames(DeMartino_meta) <- DeMartino_meta$Run

Meta_all[["DeMartinoP_2022"]] <- DeMartino_meta


# Devarakonda
cat("DevarakondaSLS_2024\n")
DevarakondaSLS_2024_meta <- read.csv("DevarakondaSLS_2024/DevarakondaSLS_2024_metadata.csv") %>% 
  separate(Sample.Name, sep = "TP", into = c("Individual", "Timepoint")) %>% 
  separate_wider_delim(Timepoint, ".", names = c("Timepoint", "Version"), too_few = "align_start") %>% 
  mutate(Version = replace_na(Version, "1")) %>% 
  group_by(Individual, Timepoint) %>% 
  slice_max(Version, n = 1) %>% 
  ungroup() %>% 
  merge(read.csv("DevarakondaSLS_2024/DevarakondaSLS_2024_Intervention.csv"), by.x = "Individual", by.y = "participant_id") %>% 
  filter(Timepoint %in% c("02", "04", "10")) %>% 
  rename(Intervention = Int1) %>% 
  group_by(Individual) %>% filter(n() == 3) %>% 
  ungroup() %>% 
  select(Run, Individual, Timepoint, Intervention) %>% 
  mutate(Intervention = case_when(Timepoint == "10" & Intervention == "RS2" ~ "RS4",
                                  Timepoint == "10" & Intervention == "RS4" ~ "RS2",
                                  .default = Intervention),
         Individual = as.character(Individual))

## RS2
DevarakondaSLS_2024_RS2_meta <- DevarakondaSLS_2024_meta %>% 
  filter(Intervention == "RS2" | Timepoint == "02") %>% 
  mutate(Intervention = "RS2",
         Timepoint = ifelse(Timepoint == "02", "Before", "After"),
         Run = str_c("DevarakondaSLS_2024_RS2", Run, sep = "_"))

rownames(DevarakondaSLS_2024_RS2_meta) <- DevarakondaSLS_2024_RS2_meta$Run
Meta_all[["DevarakondaSLS_2024_RS2"]] <- DevarakondaSLS_2024_RS2_meta

## RS4
DevarakondaSLS_2024_RS4_meta <- DevarakondaSLS_2024_meta %>% 
  filter(Intervention == "RS4" | Timepoint == "02") %>% 
  mutate(Intervention = "RS4",
         Timepoint = ifelse(Timepoint == "02", "Before", "After"),
         Run = str_c("DevarakondaSLS_2024_RS4", Run, sep = "_"))

rownames(DevarakondaSLS_2024_RS4_meta) <- DevarakondaSLS_2024_RS4_meta$Run
Meta_all[["DevarakondaSLS_2024_RS4"]] <- DevarakondaSLS_2024_RS4_meta


# Healey 
cat("HealeyG_2018\n")
Healey_meta <- read_delim("HealeyG_2018/HealeyG_2018_metadata.txt") %>% 
  filter(treatment == "fiber") %>% 
  mutate(Individual = as.character(gsub("healey_", "", subject_id)),
         Timepoint = ifelse(timepoint == "before", "Before", "After"),
         Intervention = "Fructans",
         Run = str_c("HealeyG_2018", sample_id_2, sep = "_")) %>% 
  select(Run, Individual, Timepoint, Intervention)
rownames(Healey_meta) <- Healey_meta$Run

Meta_all[["HealeyG_2018"]] <- Healey_meta


# Inoue
Inoue_meta <- read.csv("InoueR_2025/InoueR_2025_Metadata.csv") %>% 
  filter(grepl("Treatment", replicate)) %>%
  filter(!grepl("2wk", replicate)) %>% 
  separate(replicate, sep = "-", into = c("Timepoint", "tmp")) %>% 
  separate(tmp, sep = "_", into = c("tmp", "Individual")) %>% 
  mutate(Timepoint = ifelse(Timepoint == "4wk", "After", "Before"),
         Intervention = "General Fibers",
         Run = str_c("InoueR_2025", Run, sep = "_"),
         Individual = as.character(Individual)) %>% 
  select(Run, Individual, Timepoint, Intervention)
rownames(Inoue_meta) <- Inoue_meta$Run

Meta_all[["InoueR_2025"]] <- Inoue_meta


# Kordowski
cat("KordowskiA_2022\n")
Kordowski_meta <- read.csv("KordowskiA_2022/KordowskiA_2022_metadata.csv") %>%
  separate(Submitter_Id, sep = "(?<=^\\d)(?=\\d+)", into = c("Timepoint", "Individual")) %>%
  mutate(Intervention = "Fructans",
         Timepoint = ifelse(Timepoint %in% c("1", "2"), "Before", "After"),
         Run = str_c("KordowskiA_2022", Run, sep = "_"),
         Individual = as.character(Individual)) %>%
  group_by(Individual) %>% filter(n() == 2) %>%
  ungroup() %>%
  select(Run, Individual, Timepoint, Intervention)
rownames(Kordowski_meta) <- Kordowski_meta$Run

Meta_all[["KordowskiA_2022"]] <- Kordowski_meta


# Korem
cat("KoremT_2017\n")
KoremT_2017_meta <- read.csv("KoremT_2017/KoremT_2017_metadata.tsv", sep = "\t") %>%
  filter(sample_alias != "Bread16S") %>% 
  separate(sample_alias, sep = "_", into = c("Individual", "tmp", "Intervention", "Timepoint")) %>%
  mutate(Individual = as.character(gsub("Subj", "", Individual)),
         Timepoint = ifelse(Timepoint == "Start", "Before", "After"),
         Run = str_c("KoremT_2017", run_accession, sep = "_")) %>%
  select(Run, Individual, Timepoint, Intervention) %>% 
  filter(Intervention == "Srdgh") %>%
  mutate(Intervention = "Arabinoxylan") %>% 
  group_by(Individual) %>%
  filter(n() == 2)
rownames(KoremT_2017_meta) <- KoremT_2017_meta$Run

Meta_all[["KoremT_2017"]] <- KoremT_2017_meta


# Morales
cat("MoralesP_2016\n")
MoralesP_2016_meta <- read.csv("MoralesP_2016/MoralesP_2016_metadata.csv") %>% 
  mutate(Run = str_c("MoralesP_2016", Run, sep = "_"),
         Intervention = "Fructans",
         Individual = as.character(Individual))
rownames(MoralesP_2016_meta) <- MoralesP_2016_meta$Run

Meta_all[["MoralesP_2016"]] <- MoralesP_2016_meta


# Nishimoto 
cat("NishimotoY_2023\n")
Nishimoto_meta <- read.csv("NishimotoY_2023/NishimotoY_2023_Metadata.csv") %>% 
  filter(grepl("T", Sample_name)) %>% 
  separate(Sample_name, sep = "_", into = c("tmp", "Individual", "Timepoint")) %>% 
  mutate(Run = str_c("NishimotoY_2023", Run, sep = "_"),
         Timepoint = ifelse(Timepoint == "0w", "Before", "After"),
         Intervention = "Mushroom",
         Individual = as.character(Individual)) %>% 
  group_by(Individual) %>% filter(n() == 2) %>% 
  select(Run, Individual, Timepoint, Intervention)
rownames(Nishimoto_meta) <- Nishimoto_meta$Run

Meta_all[["NishimotoY_2023"]] <- Nishimoto_meta


# NiY
cat("NiY_2023\n")
NiY_2023_meta <- read.csv("NiY_2023/NiY_2023_metadata.csv") %>%
  separate(Filename, sep = "V", remove = F, into = c("tmp", "Timepoint")) %>%
  separate(tmp, sep = "(?<=S)", into = c("Intervention", "Individual")) %>%
  group_by(Individual) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  select(Run, Individual, Intervention, Timepoint) %>% 
  filter(Intervention == "RS") %>%
  mutate(Timepoint = ifelse(Timepoint == "1", "Before", "After"),
         Run = str_c("NiY_2023", Run, sep = "_"),
         Intervention = "RS2",
         Individual = as.character(Individual))
rownames(NiY_2023_meta) <- NiY_2023_meta$Run

Meta_all[["NiY_2023"]] <- NiY_2023_meta


# Rehner
cat("RehnerJ_2023\n")
RehnerJ_2023_meta <- read.xlsx("RehnerJ_2023/RehnerJ_2023_metadata.xlsx")

RehnerJ_2023_meta_expanded <- expand.grid(Individual = RehnerJ_2023_meta$Sample_ID, Timepoint = c("00", "12")) %>%
  filter(!(Individual %in% c("AB01", "AB32", "AB33", "AB26"))) %>% # Removing individuals that don't have all 4 timepoints  
  mutate(Run = paste0("RehnerJ_2023", "_", Individual, "_", Timepoint)) %>%
  merge(RehnerJ_2023_meta, by.x = "Individual", by.y = "Sample_ID", all.x = T) %>%
  filter(Group != "omnivorous") %>%
  rename(Intervention = Group) %>%
  mutate(Timepoint = ifelse(Timepoint == "00", "Before", "After"),
         Intervention = "General Fibers",
         Individual = as.character(Individual)) %>%
  select(Run, Individual, Intervention, Timepoint)
rownames(RehnerJ_2023_meta_expanded) <- RehnerJ_2023_meta_expanded$Run

Meta_all[["RehnerJ_2023"]] <- RehnerJ_2023_meta_expanded


# Roager
cat("RoagerHM_2019\n")
RoagerHM_2019_meta <- read.csv("RoagerHM_2019/RoagerHM_2019_metadata.csv") %>%
  filter(Assay.Type == "WGS" & comment == "Completer") %>%
  mutate(Timepoint = as.numeric(gsub(".*_", "", Sample.Name)),
         Intervention = case_when((Diet == "a" & Timepoint <= 2.1 ) ~ "Whole_Grain_1",
                                            (Diet == "a" & Timepoint > 2.1 ) ~ "Whole_Grain_2",
                                            (Diet == "b" & Timepoint <= 2.1 ) ~ "Refined_Grain_1",
                                            (Diet == "b" & Timepoint > 2.1 ) ~ "Refined_Grain_2"),
         Intervention = gsub("_", "", Intervention),
         Intervention = gsub("\\d+", "", Intervention),
         Run = str_c("RoagerHM_2019", Run, sep = "_")) %>%
  rename(Individual = individual) %>%
  select(Run, Individual, Timepoint, Intervention) %>% 
  filter(Intervention == "WholeGrain" & Individual != "5385") %>%
  mutate(Timepoint = case_when(Timepoint == "1" ~ "Before",
                               Timepoint == "3" ~ "Before",
                               Timepoint == "2" ~ "After",
                               Timepoint == "4" ~ "After"),
         Intervention = "Arabinoxylan",
         Individual = as.character(Individual))
rownames(RoagerHM_2019_meta) <- RoagerHM_2019_meta$Run

Meta_all[["RoagerHM_2019"]] <- RoagerHM_2019_meta


# Venkataraman
cat("VenkataramanA_2016\n")
Venkataraman_meta <- read_delim("VenkataramanA_2016/VenkataramanA_2016_metadata.txt") %>% 
  filter(treatment == "fiber") %>% 
  group_by(subject_id) %>% 
  filter(row_number() == 1 | row_number() == n()) %>% 
  ungroup() %>% 
  mutate(Individual = as.character(gsub("venkataraman_U0", "", subject_id)),
         Timepoint = ifelse(timepoint == "before", "Before", "After"),
         Intervention = "RS2",
         Run = str_c("VenkataramanA_2016", sample_id_2, sep = "_")) %>% 
  select(Run, Individual, Timepoint, Intervention)
rownames(Venkataraman_meta) <- Venkataraman_meta$Run

Meta_all[["VenkataramanA_2016"]] <- Venkataraman_meta


# Vital
cat("VitalM_2018\n")
Vital_samples <- read.xlsx("VitalM_2018/VitalM_2018_metadata_intervention.xlsx") 
  
Vital_meta <- read.csv("VitalM_2018/VitalM_2018_metadata.csv") %>% 
  separate(Submitter_Id, sep = "_", into = c("sample", "tmp1", "tmp2", "tmp3", "Timepoint")) %>% 
  merge(Vital_samples, by.x = "sample", by.y = "Sample.name") %>% 
  separate(PAIENT.ID, sep = "-", into = c("Intervention", "Individual", "tmp4")) %>% 
  filter(DIET_KEY == "highRS" & Timepoint %in% c("001", "004")) %>% 
  mutate(Timepoint = ifelse(Timepoint == "001", "Before", "After"),
         Intervention = "RS2",
         Run = str_c("VitalM_2018", Run, sep = "_"),
         Individual = as.character(Individual)) %>% 
  select(Run, Individual, Timepoint, Intervention)
rownames(Vital_meta) <- Vital_meta$Run

Meta_all[["VitalM_2018"]] <- Vital_meta


# Walsh
cat("WalshLH_2023\n")
WalshLH_2023_meta <- read.csv("WalshLH_2023/WalshLH_2023_metadata.tsv", sep = "\t") %>%
  filter(library_strategy == "OTHER") %>% 
  separate(sample_alias, sep = "_", into = c("tmp", "Intervention", "Timepoint", "Individual")) %>%
  mutate(Timepoint = ifelse(Timepoint == "pre", "Before", "After"),
         Run = str_c("WalshLH_2023", run_accession, sep = "_")) %>%
  filter(Intervention == "inulin") %>%
  mutate(Intervention = "Fructans",
         Individual = as.character(Individual)) %>% 
  group_by(Individual) %>%
  filter(n() == 2) %>%
  ungroup() %>% 
  select(Run, Individual, Intervention, Timepoint)
rownames(WalshLH_2023_meta) <- WalshLH_2023_meta$Run

Meta_all[["WalshLH_2023"]] <- WalshLH_2023_meta


# Wastyk
cat("WastykHC_2021\n")
WastykHC_2021_meta <- as(read.csv("WastykHC_2021/WastykHC_2021_metadata.csv"), "data.frame") %>%
  filter(Assay.Type == "WGS") %>% # Filter only shotgun samples
  separate(Sample.Name, sep = "_", into = c("tmp1", "Individual", "Timepoint", "tmp2")) %>%
  rename(Intervention = host_diet) %>%
  mutate(Run = str_c("WastykHC_2021", Library.Name, sep = "_")) %>%
  select(Run, Individual, Intervention, Timepoint) %>%
  filter(!startsWith(Run, "Pilot") & (Timepoint == 1 | Timepoint == 7) & Intervention == "Fiber" & Individual != 8001) %>%
  group_by(Individual, Timepoint) %>%
  slice(1) %>%
  ungroup() %>% 
  mutate(Timepoint = ifelse(Timepoint == 1, "Before", "After"),
         Intervention = "General Fibers",
         Individual = as.character(Individual)) %>%
  group_by(Individual) %>%
  filter(n_distinct(Timepoint) == 2) %>%
  ungroup()
rownames(WastykHC_2021_meta) <- WastykHC_2021_meta$Run

Meta_all[["WastykHC_2021"]] <- WastykHC_2021_meta



shotgun_studies <- read.csv("Query_Final.csv") %>% 
  filter(Notes == "Shotgun") %>% 
  pull(First.Author) %>% 
  gsub("-", "", .)

for (study_name in names(Meta_all)) {
  
  Meta_all[[study_name]] <- Meta_all[[study_name]] %>% 
    mutate(SeqMethod = ifelse(study_name %in% shotgun_studies, "Shotgun", "16S"))
  
}



saveRDS(Meta_all, "../preprocessing/data/Metadata.RDS")
cat("\nOperation Completed\nSaving and Exiting\n")
