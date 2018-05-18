### Create graph of SCFA classification AUCs
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs("tidyverse")


# Read in the data
classification_model_data <- read_csv("data/process/tables/acetate_classification_RF_summary.csv") %>% 
  mutate(scfa = "acetate") %>% 
  bind_rows(
    read_csv("data/process/tables/butyrate_classification_RF_summary.csv") %>% 
      mutate(scfa = "butyrate"), 
    read_csv("data/process/tables/propionate_classification_RF_summary.csv") %>% 
      mutate(scfa = "propionate"))

classification_summary <- classification_model_data %>% 
  group_by(scfa) %>% 
  summarise(median_train_auc = median(train_AUC), 
            max_train_auc = max(train_AUC), 
            min_train_auc = min(train_AUC), 
            median_test_auc = median(test_AUC), 
            max_test_auc = max(test_AUC), 
            min_test_auc = min(test_AUC), 
            median_sens = median(Sens), max_sens = max(Sens), min_sens = min(Sens), 
            median_spec = median(Spec), max_spec = max(Spec), min_spec = min(Spec)) %>% 
  gather("groupings", "auc", median_train_auc:min_test_auc)


classification_imp_otus <- read_csv("data/process/tables/acetateimp_otus_classification_RF_summary.csv") %>% 
  mutate(scfa = "acetate") %>% 
  bind_rows(
    read_csv("data/process/tables/butyrateimp_otus_classification_RF_summary.csv") %>% 
      mutate(scfa = "butyrate"), 
    read_csv("data/process/tables/propionateimp_otus_classification_RF_summary.csv") %>% 
      mutate(scfa = "propionate")) %>% 
  group_by(scfa, otu) %>% 
  summarise(median_mda = median(Overall), 
            max_mda = max(Overall), 
            min_mda = min(Overall)) %>% 
  arrange(desc(median_mda), .by_group = T) %>% 
  slice(1:10)

# Read in taxnomic data
tax <- read_tsv("data/process/final.taxonomy") %>% 
  mutate(Taxonomy = str_replace_all(Taxonomy, "\\(\\d*\\)", "")) %>% 
  separate(Taxonomy, c("domain", "phyla", "class", "order", "family", "genus", "species"), sep = ";") %>% 
  mutate(genus = str_replace_all(genus, "_unclassified", ""), 
         genus = str_replace_all(genus, "_", " "))

combined_imp_class_summary <- classification_imp_otus %>% 
  left_join(select(tax, OTU, genus), by = c("otu" = "OTU"))


### Plot the classification AUC data###

classification_summary %>% 
  select(scfa, groupings, auc) %>% 
  separate(groupings, c("range_data", "groups"), sep = "_") %>% 
  spread(range_data, auc) %>% 
  mutate(groups = factor(groups, 
                         levels = c("train", "test"), 
                         labels = c("Train", "Test")), 
         scfa = factor(scfa, 
                       levels = c("acetate", "butyrate", "propionate"), 
                       labels = c("Acetate", "Butyrate", "Propionate"))) %>% 
  ggplot(aes(scfa, median, color = groups, group = groups)) + 
  geom_pointrange(aes(ymin = min, ymax = max), position = position_dodge(width = 0.5), size = 1) + 
  geom_vline(xintercept=seq(1.5, length(unique(classification_summary$scfa))-0.5, 1), 
             lwd=1, colour="gray", alpha = 0.6) + 
  theme_bw() + coord_cartesian(ylim = c(0, 0.8)) + 
  labs(x = "", y = "AUC") + 
  scale_color_manual(name = "Group", 
                     values = c("black", "darkgray")) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))


combined_imp_class_summary %>% 
  ungroup() %>% 
  mutate(scfa = factor(scfa, 
                       levels = c("acetate", "butyrate", "propionate"), 
                       labels = c("Acetate", "Butyrate", "Propionate"))) %>% 
  ggplot(aes(scfa, median_mda, color = genus, group = otu)) + 
  geom_pointrange(aes(ymin = min_mda, ymax = max_mda), position = position_dodge(width = 1), size = 1) + 
  geom_vline(xintercept=seq(1.5, length(unique(classification_imp_otus$scfa))-0.5, 1), 
             lwd=1, colour="gray", alpha = 0.6) + 
  theme_bw() + labs(x = "", y = "Median MDA") + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom", 
        legend.title = element_blank(),
        legend.text = element_text(face = "italic"), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))






