### Create graph of SCFA classification AUCs
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "gridExtra"))

#setup variables that will be used
scfas <- c("acetate", "butyrate", "propionate")


# Read in classification data
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


test <- classification_model_data %>% 
  gather("groupings", "auc", train_AUC:test_AUC) %>% 
  separate(groupings, c("group", "metric"), sep = "_") %>% 
  select(-metric)


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
  left_join(select(tax, OTU, family), by = c("otu" = "OTU")) %>% 
  mutate(family = str_replace_all(family, "_unclassified", ""))



# Read in the data
# The data is obtained from the model tracking of the run_rf_scfa_predictions.R script
# Put in a temporary run number until the analysis is re-run with it in the data frame.
class_data <- map(scfas, function(x) 
  read_csv(paste("data/process/tables/", x, "_classification_RF_train_probs_summary.csv", sep = "")) %>% 
    mutate(scfa = x, 
           model = "class_train")) %>% bind_rows() %>% 
  bind_rows(map(scfas, function(x) 
    read_csv(paste("data/process/tables/", x, "_classification_RF_test_probs_summary.csv", sep = "")) %>% 
      mutate(run = rep(1:100, each = 85), 
             scfa = x, 
             model = "class_test")) %>% 
      bind_rows() %>% 
      rename(pred = tempPredictions, obs = high_low, sample_id = Group)) %>% 
  mutate(correct_class = case_when(
    pred == obs ~ "yes", 
    TRUE ~ "no")) %>% 
  group_by(model, scfa, run, dx) %>% 
  summarise(yes = table(correct_class)[2]/(table(correct_class)[2] + table(correct_class)[1]), 
            no = table(correct_class)[1]/(table(correct_class)[2] + table(correct_class)[1])) %>% 
  ungroup() %>% 
  group_by(model, scfa, dx) %>% 
  summarise(median_yes = median(yes, na.rm = T), 
            median_no = median(no, na.rm = T), 
            min_yes = min(yes, na.rm = T), 
            max_yes = max(yes, na.rm = T))


### Plot the classification AUC data ###
auc_classification <- classification_model_data %>% 
  gather("groupings", "auc", train_AUC:test_AUC) %>% 
  separate(groupings, c("group", "metric"), sep = "_") %>% 
  select(-metric) %>% 
  mutate(group = factor(group, 
                        levels = c("train", "test"), 
                        labels = c("Train", "Test")), 
         scfa = factor(scfa, 
                       levels = c("acetate", "butyrate", "propionate"), 
                       labels = c("Acetate", "Butyrate", "Propionate"))) %>% 
  ggplot(aes(scfa, auc, fill = group)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_vline(xintercept=seq(1.5, length(unique(classification_summary$scfa))-0.5, 1), 
             lwd=1, colour="gray", alpha = 0.6) + 
  theme_bw() + coord_cartesian(ylim = c(0, 1.0)) + 
  labs(x = "", y = "AUC") + ggtitle("A") + 
  scale_fill_manual(name = "", 
                    values = c("white", "darkgray")) + 
  theme(plot.title = element_text(face="bold", hjust = -0.20, size = 20), 
        legend.position = "bottom", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))

# Plot the top 10 OTUs in the model graph
acetate_imp_graph <- combined_imp_class_summary %>% 
  ungroup() %>% 
  filter(scfa == "acetate") %>% 
  mutate(otu = factor(otu, 
                      levels = rev(otu), 
                      labels = rev(otu)), 
         family = case_when(
           family %in% c("Firmicutes", "Pasteurellaceae", "Verrucomicrobiaceae") ~ "Other", 
           TRUE ~ family), 
         family = factor(family, 
                         levels = c("Clostridiales", "Lachnospiraceae", "Ruminococcaceae", "Other"), 
                         labels = c("Clostridiales", "Lachnospiraceae", "Ruminococcaceae", "Other"))) %>% 
  ggplot(aes(otu, median_mda, color = family)) + 
  geom_pointrange(aes(ymin = min_mda, ymax = max_mda), position = position_dodge(width = 1), size = 1) + 
  coord_flip(ylim = c(-5, 25)) + theme_bw() + ggtitle("C") + 
  annotate("text", label = paste("Acetate"), x = 4.7, y = 16.5, size = 3.5) +
  labs(x = "", y = "Mean Decrease in Accuracy") + 
  scale_color_manual(name = "", values = c('#9B30FF', '#63B8FF', '#FF7F00', '#6C7B8B')) + 
  theme(plot.title = element_text(face="bold", hjust = -0.35, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = c(0.7, 0.2), 
        legend.background = element_rect(color = "black"), 
        legend.title = element_blank(),
        legend.text = element_text(face = "italic", size = 6.5), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))

butyrate_imp_graph <- combined_imp_class_summary %>% 
  ungroup() %>% 
  filter(scfa == "butyrate") %>% 
  mutate(otu = factor(otu, 
                      levels = rev(otu), 
                      labels = rev(otu)), 
         family = case_when(
           family %in% c("Firmicutes", "Pasteurellaceae", "Verrucomicrobiaceae") ~ "Other", 
           TRUE ~ family), 
         family = factor(family, 
                         levels = c("Clostridiales", "Lachnospiraceae", "Ruminococcaceae", "Other"), 
                         labels = c("Clostridiales", "Lachnospiraceae", "Ruminococcaceae", "Other"))) %>% 
  ggplot(aes(otu, median_mda, color = family)) + 
  geom_pointrange(aes(ymin = min_mda, ymax = max_mda), position = position_dodge(width = 1), size = 1) + 
  coord_flip(ylim = c(-5, 25)) + theme_bw() + 
  annotate("text", label = paste("Butyrate"), x = 4.7, y = 16.5, size = 3.5) + 
  labs(x = "", y = "Mean Decrease in Accuracy") + ggtitle("D") + 
  scale_color_manual(name = "", values = c('#9B30FF', '#63B8FF', '#FF7F00', '#6C7B8B')) + 
  theme(plot.title = element_text(face="bold", hjust = -0.35, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = c(0.7, 0.2), 
        legend.background = element_rect(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "italic", size = 6.5), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))

propionate_imp_graph <- combined_imp_class_summary %>% 
  ungroup() %>% 
  filter(scfa == "propionate") %>% 
  mutate(otu = factor(otu, 
                      levels = rev(otu), 
                      labels = rev(otu)), 
         family = case_when(
           family %in% c("Firmicutes", "Pasteurellaceae", "Verrucomicrobiaceae") ~ "Other", 
           TRUE ~ family), 
         family = factor(family, 
                         levels = c("Clostridiales", "Lachnospiraceae", "Ruminococcaceae", "Other"), 
                         labels = c("Clostridiales", "Lachnospiraceae", "Ruminococcaceae", "Other"))) %>% 
  ggplot(aes(otu, median_mda, color = family)) + 
  geom_pointrange(aes(ymin = min_mda, ymax = max_mda), position = position_dodge(width = 1), size = 1) + 
  coord_flip(ylim = c(-5, 25)) + theme_bw() + 
  annotate("text", label = paste("Propionate"), x = 4.7, y = 16.5, size = 3.5) + 
  labs(x = "", y = "Mean Decrease in Accuracy") + ggtitle("E") + 
  scale_color_manual(name = "", values = c('#9B30FF', '#63B8FF', '#FF7F00', '#6C7B8B')) + 
  theme(plot.title = element_text(face="bold", hjust = -0.35, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = c(0.7, 0.2), 
        legend.background = element_rect(color = "black"), 
        legend.title = element_blank(),
        legend.text = element_text(face = "italic", size = 6.5), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))


# imp_otus_classification <- combined_imp_class_summary %>% 
#   ungroup() %>% 
#   mutate(scfa = factor(scfa, 
#                        levels = c("acetate", "butyrate", "propionate"), 
#                        labels = c("Acetate", "Butyrate", "Propionate"))) %>% 
#   ggplot(aes(scfa, median_mda, color = family, group = otu)) + 
#   geom_pointrange(aes(ymin = min_mda, ymax = max_mda), position = position_dodge(width = 1), size = 1) + 
#   geom_vline(xintercept=seq(1.5, length(unique(classification_imp_otus$scfa))-0.5, 1), 
#              lwd=1, colour="gray", alpha = 0.6) + 
#   theme_bw() + labs(x = "", y = "Median MDA") + ggtitle("C") + 
#   scale_color_manual(name = "", values = c("#B0171F", "#FF6EB4", "#9B30FF", "#4169E1", 
#                                            "#63B8FF", "#FFD700", "#00FF00", "#FF7F00")) + 
#   theme(plot.title = element_text(face="bold", hjust = -0.04, size = 20), 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         legend.position = "bottom", 
#         legend.title = element_blank(),
#         legend.text = element_text(face = "italic"), 
#         axis.text.y = element_text(size = 10), 
#         axis.title = element_text(size = 12), 
#         axis.text.x = element_text(size = 10, face = "bold"))



### Create the classification graph
test_data_graph <- class_data %>% 
  ungroup() %>% 
  mutate(model = factor(model, 
                        levels = c("class_train", "class_test"), 
                        labels = c("Training (Classification)", "Testing (Classification)")), 
         dx = factor(dx, 
                     levels = c("normal", "adenoma", "cancer"), 
                     labels = c("Normal", "Adenoma", "Cancer")), 
         scfa = factor(scfa, 
                       levels = c("acetate", "butyrate", "propionate"), 
                       labels = c("Acetate", "Butyrate", "Propionate"))) %>% 
  ggplot(aes(scfa, median_yes, color = dx, group = dx)) + 
  geom_pointrange(aes(ymin = min_yes, ymax = max_yes), position = position_dodge(width = 0.6)) +  
  facet_wrap(~model) + theme_bw() + ggtitle("B") + 
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "", y = "Correct Classification Probability") + 
  scale_color_manual(name = "", values = c('#228B22', '#FFD700', '#DC143C')) + 
  theme(plot.title = element_text(face="bold", hjust = -0.09, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom", 
        legend.text = element_text(size = 10, face = "bold"), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))


### Create a merged graph
prediction_plot <- grid.arrange(auc_classification, test_data_graph, 
                                acetate_imp_graph, butyrate_imp_graph, propionate_imp_graph, 
                                layout_matrix = rbind(c(1, 2, 2), c(3, 4, 5)))

# Write out to specific directory
ggsave("results/figures/FigureS1.pdf", prediction_plot, width = 11, height = 8)







