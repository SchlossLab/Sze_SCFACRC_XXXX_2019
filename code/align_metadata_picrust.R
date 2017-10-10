### Align Metadata with shared file 
### Preps a metadata file for creation of a biom file for picrust
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse"))

## Read in needed data tables

metaI <- read_csv("data/raw/metadata/metaI_final.csv") %>% 
  select(sample, Dx_Bin, dx, fit_result) %>% 
  rename(Group = sample) %>% 
  mutate()

metaF <- read_csv("data/raw/metadata/good_metaf_final.csv") %>% 
  gather(key = "time_point", value = "Group", initial, followUp) %>% 
  select(Group, Dx_Bin, dx, fit_result)

shared <- read_tsv("data/process/final.shared")

# Merge the metadata together into a single table
combined_meta <- metaI %>% 
  filter(!(Group %in% metaF$Group)) %>% 
  bind_rows(metaF) %>% 
  mutate(Dx_Bin = stringr::str_replace(Dx_Bin, "adv Adenoma", "adv_adenoma"), 
         Dx_Bin = stringr::str_replace(Dx_Bin, "High Risk Normal", "Normal"))

# reorder combined meta file to match the shared
combined_meta <- combined_meta %>% slice(match(shared$Group, Group))


# Write out the table to the process dir
write_tsv(combined_meta, "data/process/picrust_metadata")









