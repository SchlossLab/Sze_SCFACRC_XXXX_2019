### Graph AUC distributions 
### With and without scfas added
### Marc Sze



# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "gridExtra"))

# Load in summary data tables that will be used
adn_summary <- read_csv("data/process/tables/adn_full_AUC_model_summary.csv") %>% 
  mutate(model = "full") %>% 
  bind_rows(read_csv("data/process/tables/adn_otu_only_AUC_model_summary.csv") %>% 
              mutate(model = "otu"))

crc_summary <- read_csv("data/process/tables/crc_full_AUC_model_summary.csv") %>% 
  mutate(model = "full") %>% 
  bind_rows(read_csv("data/process/tables/crc_otu_only_AUC_model_summary.csv") %>% 
              mutate(model = "otu"))



##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

adenoma <- adn_summary %>% 
  mutate(model = factor(model, 
                        levels = c("full", "otu"), 
                        labels = c("SCFA + OTU\nModel", "OTU\nModel"))) %>% 
  ggplot(aes(model, test_auc, color = model, group = model)) + 
  geom_jitter(width = 0.2, size = 2.5, show.legend = F) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.75, width = 0.5, alpha = 0.5) + 
  coord_cartesian(ylim = c(0, 1)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  labs(x = "", y = "Area Under the Curve") + 
  theme_bw() + ggtitle("A") + 
  scale_color_manual(values = c('#EEEE00', '#FFC125')) + 
  annotate("text", label = paste("Adenoma"), x = 1.5, y = 1, size = 4) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))


carcinoma <- crc_summary %>% 
  mutate(model = factor(model, 
                        levels = c("full", "otu"), 
                        labels = c("SCFA + OTU\nModel", "OTU\nModel"))) %>% 
  ggplot(aes(model, test_auc, color = model, group = model)) + 
  geom_jitter(width = 0.2, size = 2.5, show.legend = F) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.75, width = 0.5, alpha = 0.5) + 
  coord_cartesian(ylim = c(0, 1)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  labs(x = "", y = "Area Under the Curve") + 
  theme_bw() + ggtitle("B") + 
  scale_color_manual(values = c('#FFC0CB', '#B0171F')) + 
  annotate("text", label = paste("Carcinoma"), x = 1.5, y = 1, size = 4) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))

# Combine the graphs together

combined_plot <- grid.arrange(adenoma, carcinoma, ncol = 2, nrow = 1)

# Write out to specific directory
ggsave("results/figures/rf_otu_and_scfa_aucs_graph.pdf", combined_plot, width = 6, height = 6)
