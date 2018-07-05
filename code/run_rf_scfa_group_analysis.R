### The analysis of the RF results from gruop prediction of SCFAs
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


# setup variable groups
data_models <- list(final_data_normal = c(), final_data_adenoma = c(), final_data_carcinoma = c())

imp_vars_models <- list(final_data_normal = c(), final_data_adenoma = c(), final_data_carcinoma = c())


# Load in the data
for(i in names(data_groups)){
  
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


# Use the decscending importance graph (simialr to PCR paper graph) to highlight how 
# important members are changing based on disease state

test_imp_vars_summary %>% 
  filter(scfa_model == "acetate") %>% 
  ggplot(aes(rankings, median_mda, color = genus)) + 
  geom_pointrange(aes(ymin = min_mda, ymax = max_mda)) + 
  facet_wrap(~model_group) + theme_bw() + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        legend.position = "bottom", 
        legend.text = element_text(size = 8), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())







