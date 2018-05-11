### Integration G.Hannigan OGU into workflow
### Generating the relative abundances
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "optparse"))

# Command line input 
option_list = list(
  make_option(c("-cf", "--cluster_file"), type="character", default="data/process/clustering_gt1000.csv", 
              help="Contig Cluster file (CONCOCT output)", metavar="character"),
  make_option(c("-lt", "--length_table"), type="character", default="data/process/contig_length_table.tsv", 
              help="Length of each contig [default= %default]", metavar="character"), 
  make_option(c("-iaf", "--initial_abundance_file"), type="character", default="data/process/total_contig_1to10kb_rel_abund.tsv", 
              help="Initial abundance file [default= %default]", metavar="character"), 
  make_option(c("-o", "--output_file"), type="character", default="data/process/cluster_contig_final_abund.csv", 
              help="Output file name to write to [default= %default]", metavar="character"))

opt_parser = OptionParser(option_list=option_list)
command_line_input = parse_args(opt_parser)

# Unzip needed files
system(paste("gunzip ", command_line_input[["initial_abundance_file"]], ".gz", sep = ""), intern = TRUE, ignore.stderr = TRUE)

# Read contig clustering data
cluster_data <- read_csv(command_line_input[["cluster_file"]], col_names = F) %>% 
  rename(contig = X1, cluster_group = X2)

# Read in contig length data
contig_length_data <- read_tsv(command_line_input[["length_table"]], col_names = F) %>% 
  rename(contig = X1, contig_length = X2)

# Read in raw adundance data
raw_abund_data <- read_tsv(command_line_input[["initial_abundance_file"]]) %>% 
  filter(sample_name != "sample_name")

# Create a corrected abund column
corr_adund_data <- raw_abund_data %>% 
  mutate(corr_count = (10^6 * count) / contig_length)

# Add the cluster data to the table
corr_adund_data <- corr_adund_data %>% 
  left_join(cluster_data, by = "contig")

# Summarize data down to clusters to get the corrected cluster counts
final_abund_data <- corr_adund_data %>% 
  group_by(sample_name, cluster_group) %>% 
  summarise(final_abund = sum(corr_count) / length(corr_count))

# Write out the corrected abundnace table
write_csv(final_abund_data, command_line_input[["output_file"]])

# Gzip back up used file
system(paste("gzip ", command_line_input[["initial_abundance_file"]], sep = ""), intern = TRUE, ignore.stderr = TRUE)



