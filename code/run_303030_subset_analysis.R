### Integration G.Hannigan OPF into workflow
### Use the 30/30/30 to identify if SCFAs are important
### Marc Sze



# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse"))

# Load in needed data
meta_data <- read_tsv("data/raw/metadata/303030_subset/NexteraXT003Map.tsv", col_names = F) %>% 
  select(X2, X22, X28, X30) %>% 
  rename(hannigan_name = X2, sample_id = X22, sample_type = X28, disease = X30)

opf_data <- read_tsv("data/process/ClusteredOpfAbund.tsv") %>% 
  rename(gene_name = V1, hannigan_name = V2, cor_counts = sum)
















