### Compare the number of significant correlations that are in the top 10 variables 
### Get a count of whether SCFAs replace OTUs associated with them
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "gridExtra"))

#setup variables that will be used
scfas <- c("acetate", "butyrate", "propionate")

# Read in first data table
sig_reg <- read_csv("data/process/tables/significant_reg_otu_comp_summary.csv") %>% 
  left_join(read_tsv("data/process/final.taxonomy") %>% 
              select(-Size) %>% 
              mutate(Taxonomy = str_replace_all(Taxonomy, "\\(\\d*\\)", "")) %>% 
              separate(Taxonomy, c("kingdom", "phyla", "class", "order", "family", "genus", "species"), sep = ";") %>% 
              select(OTU, family), by = c("otu" = "OTU")) %>% 
  select(otu, family, genus, dx, scfa, estimate, p.value, bh) %>% 
  mutate_at(vars(p.value:bh), function(x) format(x, scientific = T, digits = 3)) %>% 
  mutate(estimate = round(estimate, digits = 3), 
         family = str_replace_all(family, "_unclassified", ""), 
         family = str_replace_all(family, "_", " "), 
         genus = str_replace_all(family, "_unclassified", ""), 
         genus = str_replace_all(family, "_", " "))

# Read in second data table
taxonomy <- read_tsv("data/process/final.taxonomy") %>% 
  select(-Size) %>% 
  mutate(Taxonomy = str_replace_all(Taxonomy, "\\(\\d*\\)", "")) %>% 
  separate(Taxonomy, c("kingdom", "phyla", "class", "order", "family", "genus", "species"), sep = ";")

all_data <- read_csv("data/process/tables/adn_full_MDA_Summary.csv") %>% 
  mutate(model = "adn_full") %>% 
  filter(Variable != "isobutyrate") %>% 
  slice(1:10) %>% 
  bind_rows(read_csv("data/process/tables/crc_full_MDA_Summary.csv") %>% 
              mutate(model = "crc_full") %>% 
              filter(Variable != "isobutyrate") %>% 
              slice(1:10), 
            read_csv("data/process/tables/adn_otu_only_MDA_Summary.csv") %>% 
              mutate(model = "adn_otu") %>% 
              slice(1:10), 
            read_csv("data/process/tables/crc_otu_only_MDA_Summary.csv") %>% 
              mutate(model = "crc_otu") %>% 
              slice(1:10)) %>% 
  left_join(select(taxonomy, OTU, family), by = c("Variable" = "OTU")) %>% 
  mutate(family = case_when(
    Variable %in% c("acetate", "butyrate", "propionate") ~ str_to_title(Variable), 
    TRUE ~ family), 
    family = str_replace(family, "_unclassified", ""), 
    family = str_replace_all(family, "_", " "), 
    Variable = case_when(
      Variable %in% c("acetate", "butyrate", "propionate") ~ str_to_title(Variable), 
      TRUE ~ Variable))

# Set up filter variables
adn_otus_scfas <- all_data %>% filter(model == "adn_full") %>% pull(Variable)
crc_otus_scfas <- all_data %>% filter(model == "crc_full") %>% pull(Variable)
adn_otus_only <- all_data %>% filter(model == "adn_otu") %>% pull(Variable)
crc_otus_only <- all_data %>% filter(model == "crc_otu") %>% pull(Variable)

# Set up read out table
read_out_data <- sig_reg %>% filter(otu %in% adn_otus_scfas) %>% mutate(model = "adn_full") %>% 
  bind_rows(sig_reg %>% filter(otu %in% adn_otus_only) %>% mutate(model = "adn_otu"), 
            sig_reg %>% filter(otu %in% crc_otus_scfas) %>% mutate(model = "crc_full"), 
            sig_reg %>% filter(otu %in% crc_otus_only) %>% mutate(model = "crc_otu"))


# Write out the csv file
write_csv(read_out_data, "data/process/tables/sig_association_and_RF_top10.csv")



