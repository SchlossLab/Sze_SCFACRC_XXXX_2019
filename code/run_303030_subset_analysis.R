### Integration G.Hannigan OPF into workflow
### Use the 30/30/30 to identify if SCFAs are important
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "caret"))

# Load in needed data
meta_data <- read_tsv("data/process/sra_meta_conversion.txt") %>% 
  filter(samp_mat_process != "Virome") %>% 
  select(Run, Sample_Name, DiseaseClass) %>% 
  rename(sample_id = Run, hannigan_name = Sample_Name, disease = DiseaseClass)
  
opf_data <- read_tsv("data/process/opf_shared.tsv")

# Get clusters that are not singletons
greater_than_1_count <- opf_data %>% 
  group_by(opf_cluster) %>% 
  summarise(total_counts = sum(cor_total_counts)) %>% 
  filter(total_counts > 1) %>% 
  pull(opf_cluster)

# Remove these from the data file
opf_data <- opf_data %>% 
  filter(opf_cluster %in% greater_than_1_count)

# Get counts of clusters that are due to more than 1 sample 
in_more_than_1_individual <- opf_data %>% 
  group_by(opf_cluster) %>% 
  summarise(sample_counts = n_distinct(sample_id)) %>% 
  filter(sample_counts > 1) %>% 
  pull(opf_cluster)
  
# Remove these from the data file
opf_data <- opf_data %>% 
  filter(opf_cluster %in% in_more_than_1_individual)

# Spread the data to make it amenable for near zero variance testing
test <- opf_data %>% 
  select(-sum_counts, -total_opf_members) %>% 
  spread(opf_cluster, cor_total_counts, fill = 0) %>% 
  left_join(select(meta_data, sample_id, disease), by = "sample_id") %>% 
  select(sample_id, disease, everything())
 
# Generate the near zero variance list 
nzv <- nearZeroVar(as.data.frame(select(test, -sample_id)))
# Generate the final data frame holder
filtered_df <- as.data.frame(select(test, -sample_id))
# Filter out the genes with near zero variance from the final data frame
filtered_df <- filtered_df[, -nzv]

# Create Correlation network
descrCor <- cor(select(filtered_df, -disease))
# Find and remove variables with a high correlation to each other
highlyCorDescr <- findCorrelation(descrCor, cutoff = .90)

filtered_df <- filtered_df[, -highlyCorDescr] %>% 
  mutate(sample_id = test$sample_id, 
         disease = filtered_df$disease) %>% 
  select(sample_id, disease, everything())

filtered_df <- filtered_df %>% 
  gather("opf_cluster", "corr_count", 3:length(colnames(filtered_df)))


write_tsv(filtered_df, "data/process/reduced_opf_shared.tsv")



