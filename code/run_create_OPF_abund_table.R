# Create a OPF abundance table
# Need files clu.tsv and orf_abund.tsv
# Within workflow pathing for the files are
  # data/raw/mmseq2_opf_run/clu.tsv
  # data/raw/diamond_analysis/orf_abund.tsv


# Load libraries
library(tidyverse)

# Command line input 
input <- commandArgs(TRUE)
opf_variable_name <- input[1]
orf_variable_name <- input[2]


# Load in needed data
opf_metadata <- read_tsv(opf_variable_name, col_names = F) %>% 
  rename(opf_cluster = X1, seq_name = X2) %>% select(-X3)

orf_count_data <- read_tsv(orf_variable_name, col_names = F) %>% 
  rename(total_counts = X1, seq_name = X2, sample_id = X3)

# Create a tidy OPF abundance table
combined_data <- orf_count_data %>% 
  inner_join(opf_metadata, by = "seq_name") %>% 
  group_by(sample_id, opf_cluster) %>% 
  summarise(total_counts = sum(total_counts))

# Write out the data
write_tsv(combined_data, "data/process/opf_shared.tsv")

