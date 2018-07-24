### Graph the aggregated data
### Visualize how well do the predictions do by group (normal, adenoma, cancer)?
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse"))

#setup variables that will be used
scfas <- c("acetate", "butyrate", "propionate")

# Read in taxnomic data
tax <- read_tsv("data/process/final.taxonomy") %>% 
  mutate(Taxonomy = str_replace_all(Taxonomy, "\\(\\d*\\)", "")) %>% 
  separate(Taxonomy, c("domain", "phyla", "class", "order", "family", "genus", "species"), sep = ";") %>% 
  mutate(genus = str_replace_all(genus, "_unclassified", ""), 
         genus = str_replace_all(genus, "_", " "))

# Read in the data
# The data is obtained from the model tracking of the run_rf_scfa_predictions.R script
# Put in a temporary run number until the analysis is re-run with it in the data frame.
reg_data <- map(scfas, function(x) 
  read_csv(paste("data/process/tables/", x, "_regression_RF_train_conc_summary.csv", sep = "")) %>% 
    mutate(scfa = x, 
           model = "reg_train")) %>% 
  bind_rows() %>% 
  bind_rows(map(scfas, function(x) 
    read_csv(paste("data/process/tables/", x, "_regression_RF_test_conc_summary.csv", sep = "")) %>% 
      mutate(run = rep(1:100, each = 85), 
             scfa = x, 
             model = "reg_test")) %>% 
      bind_rows() %>% 
      rename(obs = mmol_kg, pred = tempPredictions, sample_id = Group)) %>% 
  mutate(difference = obs - pred) %>% 
  group_by(model, scfa, run, dx) %>% 
  summarise(median_diff = median(difference), 
            min_dff = min(difference), 
            max_diff = max(difference)) %>% 
  ungroup() %>% 
  group_by(model, scfa, dx) %>% 
  summarise(median_yes = median(median_diff, na.rm = T), 
            min_yes = min(median_diff, na.rm = T), 
            max_yes = max(median_diff, na.rm = T))

# log_reg_data <- map(scfas, function(x) 
#   read_csv(paste("data/process/tables/log_", x, "_regression_RF_train_conc_summary.csv", sep = "")) %>% 
#     mutate(scfa = x, 
#            model = "reg_train")) %>% 
#   bind_rows() %>% 
#   bind_rows(map(scfas, function(x) 
#     read_csv(paste("data/process/tables/log_", x, "_regression_RF_test_conc_summary.csv", sep = "")) %>% 
#       mutate(run = rep(1:100, each = 85), 
#              scfa = x, 
#              model = "reg_test")) %>% 
#       bind_rows() %>% 
#       rename(obs = mmol_kg, pred = tempPredictions, sample_id = Group)) %>% 
#   mutate(difference = obs - pred) %>% 
#   group_by(model, scfa, run, dx) %>% 
#   summarise(median_diff = median(difference), 
#             min_dff = min(difference), 
#             max_diff = max(difference)) %>% 
#   ungroup()

# Read in data associated with the training models
regression_model_data <- read_csv("data/process/tables/acetate_regression_RF_summary.csv") %>%
  mutate(scfa = "acetate") %>%
  bind_rows(
    read_csv("data/process/tables/butyrate_regression_RF_summary.csv") %>%
      mutate(scfa = "butyrate"),
    read_csv("data/process/tables/propionate_regression_RF_summary.csv") %>%
      mutate(scfa = "propionate")) %>% 
  gather("model", "r2", train_r2:test_r2)

# Create the summary data
# Not needed anymore since switch to boxplot
# regression_summary <- regression_model_data %>%
#   group_by(scfa) %>%
#   summarise(median_train_r2 = median(train_r2),
#             max_train_r2 = max(train_r2),
#             min_train_r2 = min(train_r2),
#             median_test_r2 = median(test_r2),
#             max_test_r2 = max(test_r2),
#             min_test_r2 = min(test_r2),
#             median_sens = median(RMSE), max_sens = max(RMSE), min_sens = min(RMSE),
#             median_spec = median(MAE), max_spec = max(MAE), min_spec = min(MAE)) %>%
#   gather("groupings", "r2", median_train_r2:min_test_r2)
#   


#### Read in the model training important OTUs
regression_imp_otus <- read_csv("data/process/tables/acetateimp_otus_regression_RF_summary.csv") %>%
  mutate(scfa = "acetate") %>%
  bind_rows(
    read_csv("data/process/tables/butyrateimp_otus_regression_RF_summary.csv") %>%
      mutate(scfa = "butyrate"),
    read_csv("data/process/tables/propionateimp_otus_regression_RF_summary.csv") %>%
      mutate(scfa = "propionate")) %>%
  group_by(scfa, otu) %>%
  summarise(median_mda = median(Overall),
            max_mda = max(Overall),
            min_mda = min(Overall)) %>%
  arrange(desc(median_mda), .by_group = T) %>%
  slice(1:10)

combined_imp_reg_summary <- regression_imp_otus %>% 
  left_join(select(tax, OTU, family), by = c("otu" = "OTU")) %>% 
  mutate(family = str_replace_all(family, "_unclassified", ""))

### Plot the regression data ###
r2_classification <- regression_model_data %>%
  select(scfa, model, r2) %>%
  mutate(model = factor(model,
                         levels = c("train_r2", "test_r2"),
                         labels = c("Train", "Test")),
         scfa = factor(scfa,
                       levels = c("acetate", "butyrate", "propionate"),
                       labels = c("Acetate", "Butyrate", "Propionate"))) %>%
  ggplot(aes(scfa, r2, fill = model)) +
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_vline(xintercept=seq(1.5, length(unique(regression_model_data$scfa))-0.5, 1),
             lwd=1, colour="gray", alpha = 0.6) +
  theme_bw() + coord_cartesian(ylim = c(0, 0.4)) +
  labs(x = "", y = "R Squared") + ggtitle("A") + 
  scale_fill_manual(name = "", 
                    values = c("white", "darkgray")) + 
  theme(plot.title = element_text(face="bold", hjust = -0.18, size = 20),
        legend.position = "bottom", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"))

acetate_imp_graph <- combined_imp_reg_summary %>% 
  ungroup() %>% 
    filter(scfa == "acetate") %>% 
  mutate(otu = factor(otu, 
                      levels = rev(otu), 
                      labels = rev(otu)), 
    family = case_when(
           family %in% c("Firmicutes", "Pasteurellaceae") ~ "Other", 
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

butyrate_imp_graph <- combined_imp_reg_summary %>% 
  ungroup() %>% 
  filter(scfa == "butyrate") %>% 
  mutate(otu = factor(otu, 
                      levels = rev(otu), 
                      labels = rev(otu)), 
         family = case_when(
           family %in% c("Firmicutes", "Pasteurellaceae") ~ "Other", 
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

propionate_imp_graph <- combined_imp_reg_summary %>% 
  ungroup() %>% 
  filter(scfa == "propionate") %>% 
  mutate(otu = factor(otu, 
                      levels = rev(otu), 
                      labels = rev(otu)), 
         family = case_when(
           family %in% c("Firmicutes", "Pasteurellaceae") ~ "Other", 
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

  # ggplot(aes(scfa, median_mda, color = family, group = otu)) + 
  # geom_pointrange(aes(ymin = min_mda, ymax = max_mda), position = position_dodge(width = 1), size = 1) + 
  # geom_vline(xintercept=seq(1.5, length(unique(regression_imp_otus$scfa))-0.5, 1), 
  #            lwd=1, colour="gray", alpha = 0.6) + 
  # theme_bw() + labs(x = "", y = "Median MDA") + ggtitle("C") + 
  # scale_color_manual(name = "", values = c('#9B30FF', '#63B8FF', '#FF7F00', '#6C7B8B')) + 
  # theme(plot.title = element_text(face="bold", hjust = -0.04, size = 20), 
  #       panel.grid.major = element_blank(), 
  #       panel.grid.minor = element_blank(), 
  #       legend.position = "bottom", 
  #       legend.title = element_blank(),
  #       legend.text = element_text(face = "italic"), 
  #       axis.text.y = element_text(size = 10), 
  #       axis.title = element_text(size = 12), 
  #       axis.text.x = element_text(size = 10, face = "bold"))

### Create the regression graph
test_data_graph <- reg_data %>% 
  ungroup %>% 
  mutate(model = factor(model, 
                        levels = c("reg_train", "reg_test"), 
                        labels = c("Training (Regression)", "Testing (Regression)")), 
         dx = factor(dx, 
                     levels = c("normal", "adenoma", "cancer"), 
                     labels = c("Normal", "Adenoma", "Cancer")), 
         scfa = factor(scfa, 
                       levels = c("acetate", "butyrate", "propionate"), 
                       labels = c("Acetate", "Butyrate", "Propionate"))) %>% 
  ggplot(aes(scfa, median_yes, color = dx, group = dx)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") + 
  geom_pointrange(aes(ymin = min_yes, ymax = max_yes), position = position_dodge(width = 0.6)) + 
  facet_wrap(~model) + theme_bw() + ggtitle("B") + 
  coord_cartesian(ylim = c(-20, 20)) + 
  labs(x = "", y = "Difference from Actual Concentration") + 
  scale_color_manual(name = "", values = c('#228B22', '#FFD700', '#DC143C')) + 
  annotate("text", label = paste("Higher in Observed"), x = 2.9, y = 18, size = 3.5) + 
  annotate("text", label = paste("Higher in Predicted"), x = 2.9, y = -18, size = 3.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom", 
        legend.text = element_text(size = 10, face = "bold"), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))


### Create the log regression graph
# log_reg_data %>% 
#   ungroup %>% 
#   mutate(model = factor(model, 
#                         levels = c("reg_train", "reg_test"), 
#                         labels = c("Training (Regression)", "Testing (Regression)")), 
#          dx = factor(dx, 
#                      levels = c("normal", "adenoma", "cancer"), 
#                      labels = c("Normal", "Adenoma", "Cancer")), 
#          scfa = factor(scfa, 
#                        levels = c("acetate", "butyrate", "propionate"), 
#                        labels = c("Acetate", "Butyrate", "Propionate"))) %>% 
#   ggplot(aes(scfa, median_diff, color = dx, group = dx)) + 
#   geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") + 
#   geom_point(position = position_dodge(width = 0.6)) + 
#   stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
#                colour = "black", geom = "crossbar", size = 0.5, width = 0.4, 
#                position = position_dodge(width = 0.6)) + 
#   facet_wrap(~model) + theme_bw() + 
#   coord_cartesian(ylim = c(-20, 20)) + 
#   labs(x = "", y = expression(Difference~from~Actual~Log["10"]~Concentration)) + 
#   scale_color_manual(name = "", values = c('#228B22', '#FFD700', '#DC143C')) + 
#   theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         legend.position = "bottom", 
#         legend.text = element_text(size = 10, face = "bold"), 
#         axis.text.y = element_text(size = 10), 
#         axis.title = element_text(size = 12), 
#         axis.text.x = element_text(size = 10, face = "bold"))


### Create a merged graph
prediction_plot <- grid.arrange(r2_classification, test_data_graph, 
                                acetate_imp_graph, butyrate_imp_graph, propionate_imp_graph, 
                                layout_matrix = rbind(c(1, 2, 2), c(3, 4, 5)))

# Write out to specific directory
ggsave("results/figures/Figure4.pdf", prediction_plot, width = 11, height = 8)
