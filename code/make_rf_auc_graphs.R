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

# Load needed libraries
loadLibs(c("tidyverse", "gridExtra"))

taxonomy <- read_tsv("data/process/final.taxonomy") %>% 
  select(-Size) %>% 
  mutate(Taxonomy = str_replace_all(Taxonomy, "\\(\\d*\\)", "")) %>% 
  separate(Taxonomy, c("kingdom", "phyla", "class", "order", "family", "genus", "species"), sep = ";")

all_data <- read_csv("data/process/tables/adn_full_MDA_Summary.csv") %>% 
  mutate(model = "adn_full") %>% 
  filter(Variable != "isobutyrate") %>% 
  slice(1:10) %>% 
  bind_rows(read_csv("data/process/tables/crc_full_MDA_Summary.csv") %>% 
              mutate(model = "crc_full") %>% 
              filter(Variable != "isobutyrate") %>% 
              slice(1:10), 
            read_csv("data/process/tables/adn_otu_only_MDA_Summary.csv") %>% 
              mutate(model = "adn_otu") %>% 
              slice(1:10), 
            read_csv("data/process/tables/crc_otu_only_MDA_Summary.csv") %>% 
              mutate(model = "crc_otu") %>% 
              slice(1:10)) %>% 
  left_join(select(taxonomy, OTU, family), by = c("Variable" = "OTU")) %>% 
  mutate(family = case_when(
    Variable %in% c("acetate", "butyrate", "propionate") ~ str_to_title(Variable), 
    TRUE ~ family), 
    family = str_replace(family, "_unclassified", ""), 
    family = str_replace_all(family, "_", " "), 
    Variable = case_when(
      Variable %in% c("acetate", "butyrate", "propionate") ~ str_to_title(Variable), 
      TRUE ~ Variable))


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
  theme_bw() + ggtitle("D") + 
  scale_color_manual(values = c('#FFC0CB', '#B0171F')) + 
  annotate("text", label = paste("Carcinoma"), x = 1.5, y = 1, size = 4) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))


all_adn <- all_data %>% 
  filter(model == "adn_full") %>% 
  mutate(Variable = factor(Variable, 
                           levels = rev(Variable), 
                           labels = rev(Variable)), 
         family = str_replace(family, "Enterobacteriaceae", "Other"), 
         family = factor(family, 
                         levels = c("Acetate", "Butyrate", "Clostridiales", "Lachnospiraceae", 
                                    "Porphyromonadaceae", "Ruminococcaceae", "Other"), 
                         labels = c("Acetate", "Butyrate", "Clostridiales", "Lachnospiraceae", 
                                    "Porphyromonadaceae", "Ruminococcaceae", "Other"))) %>% 
  ggplot(aes(Variable, median_mda, color = family)) + 
  geom_pointrange(aes(ymin = iqr25, ymax = iqr75), position = position_dodge(width = 1), size = 1) + 
  coord_flip(ylim = c(0, 3)) + theme_bw() + ggtitle("B") + 
  annotate("text", label = paste("SCFA + OTU"), x = 6, y = 2, size = 3.5) + 
  labs(x = "", y = "Mean Decrease in Accuracy") + 
  scale_color_manual(name = "", values = c('#FF3E96', '#8B475D', '#9B30FF', 
                                           '#63B8FF', '#008080', '#FF7F00', '#6C7B8B'), 
                     guide = guide_legend(nrow = 4)) + 
  theme(plot.title = element_text(face="bold", hjust = -0.2, size = 20), 
        legend.position = c(0.67, 0.3), 
        legend.background = element_rect(color = "black"), 
        legend.text = element_text(size = 6.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"))


only_otu_adn <- all_data %>% 
  filter(model == "adn_otu") %>% 
  mutate(Variable = factor(Variable, 
                           levels = rev(Variable), 
                           labels = rev(Variable)), 
         family = str_replace(family, "Bacteroidaceae", "Other"), 
         family = factor(family, 
                         levels = c("Clostridiales", "Lachnospiraceae", 
                                    "Ruminococcaceae", "Other"), 
                         labels = c("Clostridiales", "Lachnospiraceae", 
                                    "Ruminococcaceae", "Other"))) %>% 
  ggplot(aes(Variable, median_mda, color = family)) + 
  geom_pointrange(aes(ymin = iqr25, ymax = iqr75), position = position_dodge(width = 1), size = 1) + 
  coord_flip(ylim = c(0, 3)) + theme_bw() + ggtitle("C") + 
  annotate("text", label = paste("OTU"), x = 5.2, y = 2, size = 3.6) + 
  labs(x = "", y = "Mean Decrease in Accuracy") + 
  scale_color_manual(name = "", values = c("#9B30FF", "#63B8FF", "#FF7F00", "#6C7B8B"), 
                     guide = guide_legend(nrow = 2)) + 
  theme(plot.title = element_text(face="bold", hjust = -0.2, size = 20), 
        legend.position = c(0.67, 0.3), 
        legend.background = element_rect(color = "black"), 
        legend.text = element_text(size = 6.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"))



all_crc <- all_data %>% 
  filter(model == "crc_full") %>% 
  mutate(Variable = factor(Variable, 
                           levels = rev(Variable), 
                           labels = rev(Variable)), 
         family = case_when(
           family %in% c("Bacillales Incertae Sedis XI", "Bacteroidaceae", 
                         "Coriobacteriaceae", "Streptococcaceae") ~ "Other", 
           TRUE ~ family), 
         family = factor(family, 
                         levels = c("Fusobacteriaceae", "Lachnospiraceae", "Porphyromonadaceae", 
                                    "Ruminococcaceae", "Other"), 
                         labels = c("Fusobacteriaceae", "Lachnospiraceae", "Porphyromonadaceae", 
                                    "Ruminococcaceae", "Other"))) %>% 
  ggplot(aes(Variable, median_mda, color = family)) + 
  geom_pointrange(aes(ymin = iqr25, ymax = iqr75), position = position_dodge(width = 1), size = 1) + 
  coord_flip(ylim = c(0, 5)) + theme_bw() + ggtitle("E") + 
  annotate("text", label = paste("SCFA + OTU"), x = 5.7, y = 3.5, size = 3.5) + 
  labs(x = "", y = "Mean Decrease in Accuracy") + 
  scale_color_manual(name = "", values = c("#000000", "#63B8FF", "#008080", "#FF7F00", "#6C7B8B"),
                     guide = guide_legend(nrow = 3)) + 
  theme(plot.title = element_text(face="bold", hjust = -0.2, size = 20), 
        legend.position = c(0.67, 0.3), 
        legend.background = element_rect(color = "black"), 
        legend.text = element_text(size = 6.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"))



only_otu_crc <- all_data %>% 
  filter(model == "crc_otu") %>% 
  mutate(Variable = factor(Variable, 
                           levels = rev(Variable), 
                           labels = rev(Variable)), 
         family = case_when(
           family %in% c("Bacillales Incertae Sedis XI", "Bacteroidaceae", 
                         "Coriobacteriaceae", "Desulfovibrionaceae", 
                         "Prevotellaceae", "Streptococcaceae") ~ "Other", 
           TRUE ~ family), 
         family = factor(family, 
                         levels = c("Lachnospiraceae", "Porphyromonadaceae", "Other"), 
                         labels = c("Lachnospiraceae", "Porphyromonadaceae", "Other"))) %>% 
  ggplot(aes(Variable, median_mda, color = family)) + 
  geom_pointrange(aes(ymin = iqr25, ymax = iqr75), position = position_dodge(width = 1), size = 1) + 
  coord_flip(ylim = c(0, 5)) + theme_bw() + ggtitle("F") + 
  annotate("text", label = paste("OTU"), x = 5.2, y = 3.4, size = 3.5) + 
  labs(x = "", y = "Mean Decrease in Accuracy") + 
  scale_color_manual(name = "", values = c("#63B8FF", "#008080", "#6C7B8B"),
                     guide = guide_legend(nrow = 2)) + 
  theme(plot.title = element_text(face="bold", hjust = -0.2, size = 20), 
        legend.position = c(0.67, 0.3), 
        legend.background = element_rect(color = "black"), 
        legend.text = element_text(size = 6.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"))

# Combine the graphs together

prediction_plot <- grid.arrange(adenoma, all_adn, only_otu_adn, 
                                carcinoma, all_crc, only_otu_crc, 
                                layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6)))

# Write out to specific directory
ggsave("results/figures/Figure2.pdf", prediction_plot, width = 15, height = 8)
