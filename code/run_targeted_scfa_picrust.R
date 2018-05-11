### Investigate how well the picrust predicted metagenomes work
### compare most differentially expressed genes by disease (kruskal - wallis)
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "biomformat", "caret", "FSA"))

library(KEGGREST)


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

# Pull specific KEGG IDs based on pathway and literature search
# Target Butyrate, Propionate, and Acetate pathway genes 
# Enzymes at the very end of the pathways

butyrate_genes <- c("K00929", "K01034", "K01035")
propionate_genes <- c("K01908", "K01895", "K00932", "K19697")
acetate_genes <- c("K18372", "K00467", "K00156", "K00925", "K01512", "K01067", "K18118", "K01026", 
                   "K01905", "K22224", "K01895", "K00128", "K14085", "K00149", "K00138")

# Select the specific data from the cleaned data set
selected_data <- cross_sectional_data %>% 
  filter(kegg_ortholog %in% c(butyrate_genes, propionate_genes, acetate_genes))

# Run the Kruskal-wallis analysis
kruskal_test <- selected_data %>% 
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

# Write out the results
write_csv(kruskal_test, "data/process/specific_scfa_kruskal_picrust_summary.csv")




