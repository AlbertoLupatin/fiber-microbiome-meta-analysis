library(tidyverse)
library(ggplot2)
library(ggembl)

fiber_data <- read.csv("../data/query/Query_Final.csv") %>% 
  filter(Filter == "Y" & !grepl("Omary", First.Author))

tmp <- readRDS("../data/preprocessing/Metadata.RDS")

N_individuals <- lapply(tmp, function(df) length(unique(df$Individual))) %>% 
  data.frame(.) %>%
  t(.) %>% 
  as.data.frame(.) %>% 
  rownames_to_column("Study") %>% 
  rename("N_Individuals" = 2)

N_samples <- lapply(tmp, function(df) length(df$Run)) %>% 
  data.frame(.) %>%
  t(.) %>% 
  as.data.frame(.) %>% 
  rownames_to_column("Study") %>% 
  rename("N_Samples" = 2)

fiber_data$N <- N_individuals$N_Individuals
fiber_data$N_Samples <- N_samples$N_Samples


N_samples_df <- fiber_data %>%
  select(First.Author, Fiber.Type, N_Samples, Notes) %>% 
  separate_longer_delim(Fiber.Type, delim = " + ") %>% 
  mutate(Study_Compact = case_when(grepl("Baxter", First.Author) ~ "B,NT_2019",
                                   grepl("Benitez", First.Author) ~ "BP,A_2019",
                                   grepl("Dahl", First.Author) ~ "D,WJ_2015",
                                   grepl("DeehanEC_2020", First.Author) ~ "D,EC_2020",
                                   grepl("DeehanEC_2024", First.Author) ~ "D,EC_2024",
                                   grepl("Delannoy-BrunoO_2021", First.Author) ~ "DB,O_2021",
                                   grepl("Delannoy-BrunoO_2022", First.Author) ~ "DB,O_2022",
                                   grepl("DeMartino", First.Author) ~ "DM,P_2022",
                                   grepl("Devarakonda", First.Author) ~ "D,SLS_2024",
                                   grepl("Roager", First.Author) ~ "R,HM_2019",
                                   grepl("Walsh", First.Author) ~ "W,LH_2023",
                                   grepl("Wastyk", First.Author) ~ "W,HC_2021",
                                   .default = gsub("[a-z-]+", ",", First.Author)))

N_samples_df %>% filter(grepl("Fructans", Fiber.Type)) %>% pull(N_Samples) %>% sum()
N_samples_df %>% filter(grepl("RS2", Fiber.Type)) %>% pull(N_Samples) %>% sum()
N_samples_df %>% filter(grepl("RS2", Fiber.Type)) %>% pull(N_Samples) %>% sum()

ggplot(N_samples_df, aes(x = Fiber.Type, y = N_Samples, group = First.Author, fill = Notes)) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  geom_text(aes(label = Study_Compact), position = position_stack(vjust = 0.5), color = "white", size = 5.5) +
  scale_fill_manual(values = c("Shotgun" = "grey", "16s" = "black")) +
  ylab("Number of Samples") +
  theme_publication() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none")

ggsave("../figures/Stacked BarPlot.pdf", device = cairo_pdf,
       width = 17, height = 8)


N_individuals_df <- fiber_data %>%
  separate_longer_delim(Fiber.Type, delim = " + ") %>% 
  group_by(Fiber.Type, First.Author) %>%
  summarise(N = sum(N), .groups = "drop") %>% 
  mutate(Study_Compact = case_when(grepl("Baxter", First.Author) ~ "B,NT_2019",
                                   grepl("Benitez", First.Author) ~ "BP,A_2019",
                                   grepl("Dahl", First.Author) ~ "D,WJ_2015",
                                   grepl("DeehanEC_2020", First.Author) ~ "D,EC_2020",
                                   grepl("DeehanEC_2024", First.Author) ~ "D,EC_2024",
                                   grepl("Delannoy-BrunoO_2021", First.Author) ~ "DB,O_2021",
                                   grepl("Delannoy-BrunoO_2022", First.Author) ~ "DB,O_2022",
                                   grepl("DeMartino", First.Author) ~ "DM,P_2022",
                                   grepl("Devarakonda", First.Author) ~ "D,SLS_2024",
                                   grepl("Roager", First.Author) ~ "R,HM_2019",
                                   grepl("Walsh", First.Author) ~ "W,LH_2023",
                                   grepl("Wastyk", First.Author) ~ "W,HC_2021",
                                   .default = gsub("[a-z-]+", ",", First.Author))) %>% 
  merge(fiber_data %>% select(First.Author, Notes))

ggplot(N_individuals_df, aes(x = Fiber.Type, y = N, group = First.Author, fill = Notes)) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  geom_text(aes(label = Study_Compact), position = position_stack(vjust = 0.5), color = "white", size = 5.5) +
  scale_fill_manual(values = c("Shotgun" = "grey", "16s" = "black")) +
  ylab("Number of Samples") +
  theme_publication() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        legend.position = "none")
