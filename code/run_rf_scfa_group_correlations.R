### Check for correlations between the top 10 OTUs of each model
### Group SCFA analysis
# Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "gridExtra"))

#setup variables that will be used
scfas <- c("acetate", "butyrate", "isobutyrate", "propionate")

# Load in the taxonomy
taxonomy <- read_tsv("data/process/final.taxonomy") %>% 
  select(-Size) %>% 
  mutate(Taxonomy = str_replace_all(Taxonomy, "\\(\\d*\\)", "")) %>% 
  separate(Taxonomy, c("kingdom", "phyla", "class", "order", "family", "genus", "species"), sep = ";")

shared <- read_tsv("data/process/final.0.03.subsample.shared")


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

# reduce the OTUs to only the ones that are in the top 10

# Run the correlation for each OTU against each other (perhaps a correlation matrix)





