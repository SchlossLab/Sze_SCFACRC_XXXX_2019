### Analysis of Relative Abundance of previous Top 10 important variables
### from previous RF models
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse"))


# Previous research top 10

# Adenoma: Lachnospiraceae, Clostridiales, Blautia, Bacteroides, Odoribacter, Ruminococcaceae, Flavonifractor, 
# Roseburia, Escherichia/Shigella, Clostridium XIVa

# Advanced Adenoma: Clostridiales, Clostridium XIVa, Roseburia, Odoribacter, Bacteroides, Ruminococcaceae, Streptococcus, 
# Blautia, Lachnospiraceae, Ruminococcus

# Carcinoma: Porphyromonas, Parvimonas, Fusobacterium, Gemella, Prevotella, Streptococcus, Coprobacillus, Pasteurellaceae, 
# Collinsella, Bilophila


shared_data <- read_tsv("data/process/final.0.03.subsample.shared")

taxonomic_data <- read_tsv("data/process/final.taxonomy") %>% 
  mutate(Taxonomy = str_replace_all(Taxonomy, "\\((\\d{2,3})\\)", ""), 
         Taxonomy = str_replace_all(Taxonomy, "_unclassified", ""), 
         Taxonomy = str_replace_all(Taxonomy, "2", "")) %>% 
  separate(Taxonomy, c("kingdom", "phyla", "class", "order", "family", "genus", "species"), sep = ";")

combined_meta_data <- read_csv("data/raw/metadata/metaI_final.csv") %>% 
  mutate(Dx_Bin = ifelse(Dx_Bin == "High Risk Normal" | Dx_Bin == "Normal", invisible("control"), 
                         ifelse(Dx_Bin == "Adenoma", invisible("adenoma"), 
                                ifelse(Dx_Bin == "adv Adenoma", invisible("adv_adenoma"), 
                                       ifelse(Dx_Bin == "Cancer", invisible("cancer"), invisible(Dx_Bin)))))) %>% 
  select(sample, Dx_Bin, dx) %>% 
  bind_rows(read_csv("data/raw/metadata/good_metaf_final.csv") %>% 
              mutate(sample = initial) %>% 
              select(sample, Dx_Bin, dx))
  


##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################













##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################
