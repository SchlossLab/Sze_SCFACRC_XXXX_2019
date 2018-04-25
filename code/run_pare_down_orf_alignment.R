### Pare down orf alignment file to only those for OPFs
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse"))

# Command line input 
input <- commandArgs(TRUE)
aligment_file_path <- input[1]
opf_file_path <- input[2]
identity_filter <- as.numeric(input[3])

# Load in needed data and remove those whose sequence ID is less than 50%
temp_align_data <- read_tsv(aligment_file_path, col_names = F) %>% 
  rename(gene_id = X1, seq_name = X2, seq_id = X3, align_length = X4, mismatch_num = X5, 
         gap_num = X6, domain_start = X7, domain_end = X8, seq_start = X9, seq_end = X10, 
         e_value = X11, bit_score = X12) %>% 
  filter(seq_id > identity_filter)

# Load in the opf data
opf_shared_data <- read_tsv(opf_file_path)

# generate combined table with only relevent IDs
combined_data <- as.data.frame(opf_shared_data) %>% 
  left_join(as.data.frame(temp_align_data), by = c("opf_cluster" = "seq_name")) %>% 
  select(sample_id, opf_cluster, sum_counts, total_opf_members, cor_total_counts, gene_id)


# write out the data
write_tsv("data/process/final_opf_shared.tsv")





