### Integration G.Hannigan OPF into workflow
### Use the 30/30/30 to identify if SCFAs are important
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "caret", "FSA"))

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

# Remove ambiguous clusters
filtered_df <- filtered_df %>% select(-Prodigal_Seq_11_1, -Prodigal_Seq_20_1, -SRR56_1, -SRR566_1)

# gather the columns back together
filtered_df <- filtered_df %>% 
  mutate(sample_id = test$sample_id) %>% 
  select(sample_id, disease, everything())

filtered_df <- filtered_df %>% 
  gather("opf_cluster", "corr_counts", 3:length(colnames(filtered_df)))

# Save memory and remove uneeded data
rm(greater_than_1_count, in_more_than_1_individual, test, opf_data)

# Run a kruskal test to find differences between specific genes
kruskal_test <- filtered_df %>% 
  mutate(disease = factor(disease, 
                          levels = c("Healthy", "Adenoma", "Cancer"), 
                          labels = c("Normal", "Adenoma", "Cancer"))) %>% 
  group_by(opf_cluster) %>% 
  nest() %>% 
  mutate(kruskal_analysis = map(data, ~kruskal.test(corr_counts ~ disease, data = .x)), 
         summary_data = map(kruskal_analysis, broom::tidy)) %>% 
  select(opf_cluster, summary_data) %>% 
  unnest(summary_data) %>% 
  select(opf_cluster, method, parameter, statistic, p.value) %>% 
  mutate(bh = p.adjust(p.value, method = "BH")) %>% 
  arrange(p.value)

# Run a dunn's test post-hoc to find which components are most significant 
# Used a P-value cutoff of 0.05 since nothing was significant after multiple comparison correction
opf_clusters_to_keep <- kruskal_test %>% 
  filter(p.value < 0.05) %>% 
  pull(opf_cluster)

dunn_test <- filtered_df %>% 
  mutate(disease = factor(disease, 
                          levels = c("Healthy", "Adenoma", "Cancer"), 
                          labels = c("Normal", "Adenoma", "Cancer"))) %>% 
  filter(opf_cluster %in% opf_clusters_to_keep) %>% 
  group_by(opf_cluster) %>% 
  nest() %>% 
  mutate(dunn_analysis = map(data, ~dunnTest(corr_counts ~ disease, data = .x, method = "bh")[["res"]])) 

dunn_testing <- dunn_test %>% select(opf_cluster, dunn_analysis) %>% unnest()

# Find all groups that have a consistent direction
all_neg <- dunn_testing %>% group_by(opf_cluster) %>% 
  count(fulfills_req = Z < 0) %>% 
  filter(n == 3 & fulfills_req == TRUE) %>% 
  pull(opf_cluster)

all_pos <- dunn_testing %>% group_by(opf_cluster) %>% 
  count(fulfills_req = Z > 0) %>% 
  filter(n == 3 & fulfills_req == TRUE) %>% 
  pull(opf_cluster)
  
direction_consistent_dunn <- dunn_testing %>% 
  filter(opf_cluster %in% c(all_neg, all_pos))

# Pull only those that have at least 2 out of 3 significant comparisons
sig_comparisons <- direction_consistent_dunn %>% 
  group_by(opf_cluster) %>% 
  count(fulfills_req = P.adj < 0.05) %>% 
  filter(n >= 2 & fulfills_req == TRUE) %>% 
  pull(opf_cluster)

sig_dir_consistent_dunn <- direction_consistent_dunn %>% 
  filter(opf_cluster %in% sig_comparisons)

# Find only those that have carcinoma and adenoma different versus normal
specific_sig_comparisons <- sig_dir_consistent_dunn %>% 
  group_by(opf_cluster) %>% 
  filter(Comparison == "Adenoma - Normal" | Comparison == "Cancer - Normal") %>% 
  count(fulfills_req = P.adj < 0.05) %>% 
  filter(n == 2 & fulfills_req == TRUE) %>% 
  pull(opf_cluster)

specific_sig_dir_cons_dunn <- sig_dir_consistent_dunn %>% 
  filter(opf_cluster %in% specific_sig_comparisons)

# Arrange by most significant groups
most_sig_compare <- specific_sig_dir_cons_dunn %>% 
  filter(Comparison != "Adenoma - Cancer") %>% 
  group_by(opf_cluster) %>% 
  summarise(most_sig_mean = mean(P.adj)) %>% 
  arrange(most_sig_mean) %>% 
  pull(opf_cluster)

ordered_sp_sig_dir_cons_dunn <- specific_sig_dir_cons_dunn %>% 
  mutate(opf_cluster = factor(opf_cluster, 
                              levels = most_sig_compare, 
                              labels = most_sig_compare)) %>% 
  arrange(opf_cluster)


# Write out the data frames of interest
write_csv(kruskal_test, "data/process/opf_kruskal_analysis.csv")

write_csv(ordered_sp_sig_dir_cons_dunn, "data/process/opf_dunn_select_results.csv")

write_csv(dunn_testing, "data/process/opf_all_dunn_results.csv")





