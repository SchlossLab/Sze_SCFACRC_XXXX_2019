### Investigate how well the picrust predicted metagenomes work
### compare most differentially expressed genes by disease (kruskal - wallis)
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "biomformat", "caret", "FSA"))


# Load in needed metadata and rename sample column
picrust_meta <- read_tsv("data/process/picrust_metadata") %>% 
  mutate(Group = as.character(Group))

# Load in more specific follow up data
metaF <- read_csv("data/raw/metadata/good_metaf_final.csv")

# Load in predicted metagenome data
all_biom_data <- read_biom("data/process/predicted_metagenomes.biom")

# Converts the needed biom data into a matrix
full_pi_data <- t(as(biom_data(all_biom_data), "matrix")) %>% 
  as.data.frame() %>% 
  mutate(sample_id = rownames(.)) %>% 
  select(sample_id, everything())

# Get the kegg mappings for each coded sample
kegg_table <- observation_metadata(all_biom_data) %>% 
  mutate(kegg_id = rownames(.)) %>% 
  select(kegg_id, everything())

# Generate the near zero variance
nzv <- nearZeroVar(full_pi_data)

# Remove these from the data set
nzv_full_pi_data <- full_pi_data[, -nzv]


# Use only cross-sectional data and add the needed meta data
cross_sectional_data <- nzv_full_pi_data %>% 
  filter(!(sample_id %in% c(metaF$initial, metaF$followUp))) %>% 
  left_join(picrust_meta, by = c("sample_id" = "Group")) %>% 
  select(sample_id, Dx_Bin, dx, everything()) %>% 
  select(-fit_result) %>% 
  gather("kegg_ortholog", "rel_abund", 4:length(colnames(.)))
  

# Run the Kruskal-wallis analysis
kruskal_test <- cross_sectional_data %>% 
  group_by(kegg_ortholog) %>% 
  mutate(dx = factor(dx, 
                     levels = c("normal", "adenoma", "cancer"), 
                     labels = c("Normal", "Adenoma", "Cancer"))) %>% 
  nest() %>% 
  mutate(kruskal_analysis = map(data, ~kruskal.test(rel_abund ~ dx, data = .x)), 
         summary_data = map(kruskal_analysis, broom::tidy)) %>% 
  select(kegg_ortholog, summary_data) %>% 
  unnest(summary_data) %>% 
  select(kegg_ortholog, method, parameter, statistic, p.value) %>% 
  mutate(bh = p.adjust(p.value, method = "BH")) %>% 
  arrange(p.value)


# Run a dunn's test post-hoc to find which components are most significant 
# Used an BH cutoff of 0.05 
picrust_genes <- kruskal_test %>% 
  filter(bh < 0.05) %>% 
  pull(kegg_ortholog)

dunn_test <- cross_sectional_data %>% 
  mutate(dx = factor(dx, 
                     levels = c("normal", "adenoma", "cancer"), 
                     labels = c("Normal", "Adenoma", "Cancer"))) %>% 
  filter(kegg_ortholog %in% picrust_genes) %>% 
  group_by(kegg_ortholog) %>% 
  nest() %>% 
  mutate(dunn_analysis = map(data, ~dunnTest(rel_abund ~ dx, data = .x, method = "bh")[["res"]])) 

dunn_testing <- dunn_test %>% select(kegg_ortholog, dunn_analysis) %>% unnest()

# Save the resulting tables for both the kruskal-wallis and dunn's post hoc test
write_csv(kruskal_test, "data/process/picrust_kruskal_summary.csv")
write_csv(dunn_testing, "data/process/picrust_dunns_post_hoc_summary.csv")











































