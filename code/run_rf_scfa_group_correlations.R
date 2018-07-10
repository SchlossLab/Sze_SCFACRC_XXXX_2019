### Check for correlations between the top 10 OTUs of each model
### Group SCFA analysis
# Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "Hmisc"))

#setup variables that will be used
scfas <- c("acetate", "butyrate", "isobutyrate", "propionate")

# Load in the needed data
taxonomy <- read_tsv("data/process/final.taxonomy") %>% 
  select(-Size) %>% 
  mutate(Taxonomy = str_replace_all(Taxonomy, "\\(\\d*\\)", "")) %>% 
  separate(Taxonomy, c("kingdom", "phyla", "class", "order", "family", "genus", "species"), sep = ";")

shared <- read_tsv("data/process/final.0.03.subsample.shared")

metadata <- read_csv("data/raw/metadata/metaI_final.csv")


# setup variable groups
data_models <- list(final_data_normal = c(), final_data_adenoma = c(), final_data_carcinoma = c())

imp_vars_models <- list(final_data_normal = c(), final_data_adenoma = c(), final_data_carcinoma = c())


# Load in the data
for(i in names(data_models)){
  
  data_models[[i]] <- sapply(scfas, function(x) read_csv(paste("data/process/tables/", i, "_",
                                                               x, "_classification_RF_group_summary.csv", sep = "")) %>% 
                               summarise(median_train_auc = median(train_AUC), 
                                         max_train_auc = max(train_AUC), 
                                         min_train_auc = min(train_AUC), 
                                         median_test_auc = median(test_AUC), 
                                         max_test_auc = max(test_AUC), 
                                         min_test_auc = min(test_AUC), 
                                         most_common_mtry = names(table(mtry))[1]) %>% 
                               mutate(scfa_model = x), simplify = F) %>% 
    bind_rows() %>% 
    mutate(model_group = i)
  
  imp_vars_models[[i]] <- sapply(scfas, function(x)
    read_csv(paste("data/process/tables/", i, "_", x, 
                   "imp_otus_classification_RF_group_summary.csv", sep = "")) %>% 
      group_by(otu) %>% 
      summarise(median_mda = median(Overall), 
                max_mda = max(Overall), 
                min_mda = min(Overall)) %>% 
      arrange(desc(median_mda)) %>% 
      inner_join(select(taxonomy, OTU, genus), by = c("otu" = "OTU")) %>% 
      ungroup() %>% 
      mutate(scfa_model = x) %>% 
      slice(1:10), simplify = F) %>% 
    bind_rows() %>% 
    mutate(model_group = i)
  
}


test_model_summary <- data_models %>% bind_rows() %>% 
  filter(scfa_model != "isobutyrate")

test_imp_vars_summary <- imp_vars_models %>% bind_rows() %>% 
  filter(scfa_model != "isobutyrate") %>% 
  mutate(genus = str_replace_all(genus, "_unclassified", ""), 
         genus = str_replace_all(genus, "_", " "), 
         rankings = rep(c(1:10), 9))


# reduce the shared to only the cross-sectional data
cross_section_shared <- shared %>% 
  filter(Group %in%  pull(metadata, sample))

# Get the top ten OTUs by scfa
acetate_otus <- imp_vars_models %>% 
  map(function(x) 
  filter(x, scfa_model == "acetate") %>% 
  pull(otu)) %>% 
  bind_rows() %>% 
  gather("group_name", "otu", final_data_normal:final_data_carcinoma) %>% 
  pull(otu)

butyrate_otus <- imp_vars_models %>% 
  map(function(x) 
    filter(x, scfa_model == "butyrate") %>% 
      pull(otu)) %>% 
  bind_rows() %>% 
  gather("group_name", "otu", final_data_normal:final_data_carcinoma) %>% 
  pull(otu)

propionate_otus <- imp_vars_models %>% 
  map(function(x) 
    filter(x, scfa_model == "propionate") %>% 
      pull(otu)) %>% 
  bind_rows() %>% 
  gather("group_name", "otu", final_data_normal:final_data_carcinoma) %>% 
  pull(otu)


# reduce the OTUs to only the ones that are in the top 10
acetate_shared <- cross_section_shared %>% 
  select(one_of(c("Group", acetate_otus)))

butyrate_shared <- cross_section_shared %>% 
  select(one_of(c("Group", butyrate_otus)))

propionate_shared <- cross_section_shared %>% 
  select(one_of(c("Group", propionate_otus)))


# Run the correlation for each OTU against each other (perhaps a correlation matrix)

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    first_otu = rownames(cormat)[row(cormat)[ut]],
    second_otu = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# Code from http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software


test <- rcorr(as.matrix(select(acetate_shared, -Group)), type = "spearman")

acetate_table <- flattenCorrMatrix(test$r, test$P) %>% 
  arrange(p) %>% 
  mutate(var1 = case_when(
    first_otu %in% acetate_otus[1:10] ~ "normal", 
    first_otu %in% acetate_otus[11:20] ~ "adenoma", 
    first_otu %in% acetate_otus[21:30] ~ "carcinoma", 
    TRUE ~ "uh oh"), 
    var2 = case_when(
      second_otu %in% acetate_otus[1:10] ~ "normal", 
      second_otu %in% acetate_otus[11:20] ~ "adenoma", 
      second_otu %in% acetate_otus[21:30] ~ "carcinoma", 
      TRUE ~ "uh oh")) %>% 
  left_join(select(taxonomy, OTU, genus), by = c("first_otu" = "OTU")) %>% 
  rename(first_genus = genus) %>% 
  mutate(first_genus = str_replace_all(first_genus, "_unclassified", "")) %>% 
  left_join(select(taxonomy, OTU, genus), by = c("second_otu" = "OTU")) %>% 
  rename(second_genus = genus) %>% 
  mutate(second_genus = str_replace_all(second_genus, "_unclassified", ""))

test <- rcorr(as.matrix(select(butyrate_shared, -Group)), type = "spearman")

butyrate_table <- flattenCorrMatrix(test$r, test$P) %>% 
  arrange(p) %>% 
  mutate(var1 = case_when(
    first_otu %in% butyrate_otus[1:10] ~ "normal", 
    first_otu %in% butyrate_otus[11:20] ~ "adenoma", 
    first_otu %in% butyrate_otus[21:30] ~ "carcinoma", 
    TRUE ~ "uh oh"), 
    var2 = case_when(
      second_otu %in% butyrate_otus[1:10] ~ "normal", 
      second_otu %in% butyrate_otus[11:20] ~ "adenoma", 
      second_otu %in% butyrate_otus[21:30] ~ "carcinoma", 
      TRUE ~ "uh oh")) %>% 
  left_join(select(taxonomy, OTU, genus), by = c("first_otu" = "OTU")) %>% 
  rename(first_genus = genus) %>% 
  mutate(first_genus = str_replace_all(first_genus, "_unclassified", "")) %>% 
  left_join(select(taxonomy, OTU, genus), by = c("second_otu" = "OTU")) %>% 
  rename(second_genus = genus) %>% 
  mutate(second_genus = str_replace_all(second_genus, "_unclassified", ""))


test <- rcorr(as.matrix(select(propionate_shared, -Group)), type = "spearman")

propionate_table <- flattenCorrMatrix(test$r, test$P) %>% 
  arrange(p) %>% 
  mutate(var1 = case_when(
    first_otu %in% propionate_otus[1:10] ~ "normal", 
    first_otu %in% propionate_otus[11:20] ~ "adenoma", 
    first_otu %in% propionate_otus[21:30] ~ "carcinoma", 
    TRUE ~ "uh oh"), 
    var2 = case_when(
      second_otu %in% propionate_otus[1:10] ~ "normal", 
      second_otu %in% propionate_otus[11:20] ~ "adenoma", 
      second_otu %in% propionate_otus[21:30] ~ "carcinoma", 
      TRUE ~ "uh oh")) %>% 
  left_join(select(taxonomy, OTU, genus), by = c("first_otu" = "OTU")) %>% 
  rename(first_genus = genus) %>% 
  mutate(first_genus = str_replace_all(first_genus, "_unclassified", "")) %>% 
  left_join(select(taxonomy, OTU, genus), by = c("second_otu" = "OTU")) %>% 
  rename(second_genus = genus) %>% 
  mutate(second_genus = str_replace_all(second_genus, "_unclassified", ""))


