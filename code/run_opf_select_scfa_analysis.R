### Selective search and test of specific OPFs
### Generating the relative abundances
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "optparse"))

# Command line input 
option_list = list(
  make_option(c("-mf", "--match_file"), type="character", default="data/process/select_scfa_opf_matches.tsv", 
              help="KEGG to contig match file name", metavar="character"),
  make_option(c("-of", "--opf_file"), type="character", default="data/process/opf_shared.tsv", 
              help="opf shared file [default= %default]", metavar="character"), 
  make_option(c("-md", "--metadata"), type="character", default="data/process/sra_meta_conversion.txt", 
              help="metadata of the samples [default= %default]", metavar="character"), 
  make_option(c("-o", "--output_file"), type="character", default="data/process/select_scfa_opf_kruskal_summary.csv", 
              help="Output kruskal summary name [default= %default]", metavar="character"), 
  make_option(c("-eo", "--extra_output"), type="character", default="data/process/select_scfa_opf_data.csv", 
              help="Output kruskal summary name [default= %default]", metavar="character"))

opt_parser = OptionParser(option_list=option_list)
command_line_input = parse_args(opt_parser)

# Unzip needed files
system(paste("gunzip ", command_line_input[["opf_file"]], ".gz", sep = ""), intern = TRUE, ignore.stderr = TRUE)

# Load in needed data
meta_data <- read_tsv(command_line_input[["metadata"]]) %>% 
  filter(samp_mat_process != "Virome") %>% 
  select(Run, Sample_Name, DiseaseClass) %>% 
  rename(sample_id = Run, hannigan_name = Sample_Name, disease = DiseaseClass)

# Load in OPF data and combine with meta_data
opf_data <- read_tsv(command_line_input[["opf_file"]]) %>% 
  left_join(meta_data, by = "sample_id")
  
# Load in the match data
match_data <- read_tsv(command_line_input[["match_file"]], col_names = F) %>% 
  rename(opf_cluster = X1, kegg_id = X2)

# Get a vector of all the opfs to keep
opfs_to_keep <- match_data %>% pull(opf_cluster)

# Filter out uneeded opf clusters
filtered_df <- opf_data %>% 
  filter(opf_cluster %in% opfs_to_keep) %>% 
  left_join(match_data, by = "opf_cluster")

# Collapse down similar kegg ids
collapsed_df <- filtered_df %>% 
  group_by(sample_id, kegg_id) %>% 
  summarise(combined_count = sum(cor_total_counts)) %>% 
  left_join(meta_data, by = "sample_id") %>% 
  ungroup()


# Run a kruskal test to find differences between specific genes
kruskal_test <- collapsed_df %>% 
  mutate(disease = factor(disease, 
                          levels = c("Healthy", "Adenoma", "Cancer"), 
                          labels = c("Normal", "Adenoma", "Cancer"))) %>% 
  group_by(kegg_id) %>% 
  nest() %>% 
  mutate(kruskal_analysis = map(data, ~kruskal.test(combined_count ~ disease, data = .x)), 
         summary_data = map(kruskal_analysis, broom::tidy)) %>% 
  select(kegg_id, summary_data) %>% 
  unnest(summary_data) %>% 
  select(kegg_id, method, parameter, statistic, p.value) %>% 
  mutate(bh = p.adjust(p.value, method = "BH")) %>% 
  arrange(p.value)

# Write out the results
write_csv(kruskal_test, command_line_input[["output_file"]])
write_csv(collapsed_df, command_line_input[["extra_output"]])

# Zip needed files
system(paste("gzip ", command_line_input[["opf_file"]], sep = ""), intern = TRUE, ignore.stderr = TRUE)



