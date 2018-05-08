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
# Used an pvalue cutoff of 0.05 to be more lenient 
picrust_genes <- kruskal_test %>% 
  filter(p.value < 0.05) %>% 
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


# Find consistent directions
all_neg <- dunn_testing %>% group_by(kegg_ortholog) %>% 
  count(fulfills_req = Z < 0) %>% 
  filter(n == 3 & fulfills_req == TRUE) %>% 
  pull(kegg_ortholog)

all_pos <- dunn_testing %>% group_by(kegg_ortholog) %>% 
  count(fulfills_req = Z > 0) %>% 
  filter(n == 3 & fulfills_req == TRUE) %>% 
  pull(kegg_ortholog)

direction_consistent_dunn <- dunn_testing %>% 
  filter(kegg_ortholog %in% c(all_neg, all_pos))

# Pull only those that have at least 2 out of 3 significant comparisons
sig_comparisons <- direction_consistent_dunn %>% 
  group_by(kegg_ortholog) %>% 
  count(fulfills_req = P.adj < 0.05) %>% 
  filter(n >= 2 & fulfills_req == TRUE) %>% 
  pull(kegg_ortholog)

sig_dir_consistent_dunn <- direction_consistent_dunn %>% 
  filter(kegg_ortholog %in% sig_comparisons)


# Find only those that have carcinoma and adenoma different versus normal
specific_sig_comparisons <- sig_dir_consistent_dunn %>% 
  group_by(kegg_ortholog) %>% 
  filter(Comparison == "Adenoma - Normal" | Comparison == "Cancer - Normal") %>% 
  count(fulfills_req = P.adj < 0.05) %>% 
  filter(n == 2 & fulfills_req == TRUE) %>% 
  pull(kegg_ortholog)


specific_sig_dir_cons_dunn <- sig_dir_consistent_dunn %>% 
  filter(kegg_ortholog %in% specific_sig_comparisons)


# Save the resulting tables for both the kruskal-wallis and dunn's post hoc test
write_csv(kruskal_test, "data/process/picrust_kruskal_summary.csv")
write_csv(dunn_testing, "data/process/picrust_dunns_post_hoc_summary.csv")


# Identify the Kegg orthologs

##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################


# Function to do the actual data pulling
pull_respective_kegg_data <- function(ID_of_int){
  
  tempKEGG <- try(keggGet(ID_of_int))
  
  tempEntry <- data_frame(kegg_id = try(unname(tempKEGG[[1]]$ENTRY)), 
          gene = try(tempKEGG[[1]]$NAME), 
          full_gene = try(paste(unname(tempKEGG[[1]]$PATHWAY), "_", 
                                collapse = "", sep = "")), 
          ortholog = try(paste(names(tempKEGG[[1]]$PATHWAY), "_", collapse = "", sep = "")))
  
  if(length(tempEntry) <= 1){
    
    tempEntry <- data_frame(kegg_id = ID_of_int, gene = NA, full_gene = NA, ortholog = NA)
  }
  
  print(paste("Completed processing KEGG ID: ", ID_of_int, sep = ""))
  
  
  
  return(tempEntry)
}



##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################


kegg_ids_of_int <- dunn_testing %>% 
  group_by(kegg_ortholog) %>% 
  nest() %>% 
  mutate(test = map(kegg_ortholog, function(x) pull_respective_kegg_data(x))) %>% 
  select(kegg_ortholog, test) %>% 
  unnest(test)
  # Note K00534 is the only gene that is up in both adenoma and carcinoma based in the same direction

pull_respective_kegg_data("K03851")

unlist(kegg_ids_of_int$test)





































