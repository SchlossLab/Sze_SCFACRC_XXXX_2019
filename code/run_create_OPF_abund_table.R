# Create a OPF abundance table
# Need files contig_length_table.tsv, clu.tsv, and orf_abund.tsv
# Within workflow pathing for the files are
  # data/raw/mmseq2_opf_run/clu.tsv
  # data/raw/diamond_analysis/orf_abund.tsv
  # data/raw/contig_length_table.tsv


# Load libraries
library(tidyverse)

# Command line input 
input <- commandArgs(TRUE)
opf_variable_name <- input[1]
orf_variable_name <- input[2]
#contig_length_name <- input[3]


# Load in needed data
opf_metadata <- read_tsv(opf_variable_name, col_names = F) %>% 
  rename(opf_cluster = X1, seq_name = X2) %>% select(-X3)

# test <- opf_metadata %>% separate(opf_cluster, c(paste("X", 1:10, sep = "")), sep = "_") %>% 
#   select(X1, X2, X3, seq_name) %>% 
#   unite(contig, X1, X2, X3, sep = "_")

#opf_metadata <- opf_metadata %>% inner_join(test, by = "seq_name")

# rm(test)

orf_count_data <- read_tsv(orf_variable_name, col_names = F) %>% 
  rename(total_counts = X1, seq_name = X2, sample_id = X3)


#contig_length_data <- read_tsv(contig_length_name, col_names = F) %>% 
#  rename(seq_name = X1, contig_length = X2)

# Create a tidy OPF abundance table
combined_data <- orf_count_data %>% 
  inner_join(opf_metadata, by = "seq_name") %>% 
#  inner_join(contig_length_data, by = c("contig" = "seq_name")) %>% 
  group_by(sample_id, opf_cluster) %>% 
  summarise(sum_counts = sum(total_counts), 
            total_opf_members = length(total_counts)) %>%
  mutate(cor_total_counts = sum_counts/total_opf_members)

# Write out the data
write_tsv(combined_data, "data/process/opf_shared.tsv")

