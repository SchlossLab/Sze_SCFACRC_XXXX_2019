

# Load in needed functions and libraries
source('code/functions.R')

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
  coord_flip(ylim = c(0, 3)) + theme_bw() + 
  labs(x = "", y = "Mean Decrease in Accuracy") + 
  scale_color_manual(name = "", values = c('#FF3E96', '#8B475D', '#9B30FF', 
                                           '#63B8FF', '#008080', '#FF7F00', '#6C7B8B'), 
                     guide = guide_legend(nrow = 4)) + 
  theme(legend.position = c(0.7, 0.3), 
        legend.background = element_rect(color = "black"), 
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
  coord_flip(ylim = c(0, 3)) + theme_bw() + 
  labs(x = "", y = "Mean Decrease in Accuracy") + 
  scale_color_manual(name = "", values = c("#9B30FF", "#63B8FF", "#FF7F00", "#6C7B8B"), 
                     guide = guide_legend(nrow = 2)) + 
  theme(legend.position = c(0.7, 0.3), 
        legend.background = element_rect(color = "black"), 
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
  coord_flip(ylim = c(0, 5)) + theme_bw() + 
  labs(x = "", y = "Mean Decrease in Accuracy") + 
  scale_color_manual(name = "", values = c("#000000", "#63B8FF", "#008080", "#FF7F00", "#6C7B8B"),
                     guide = guide_legend(nrow = 3)) + 
  theme(legend.position = c(0.7, 0.3), 
        legend.background = element_rect(color = "black"), 
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
  coord_flip(ylim = c(0, 5)) + theme_bw() + 
  labs(x = "", y = "Mean Decrease in Accuracy") + 
  scale_color_manual(name = "", values = c("#63B8FF", "#008080", "#6C7B8B"),
                     guide = guide_legend(nrow = 2)) + 
  theme(legend.position = c(0.7, 0.3), 
        legend.background = element_rect(color = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(size = 10, face = "bold"))


### Create a merged graph
prediction_plot <- grid.arrange(all_adn, only_otu_adn, all_crc, only_otu_crc, 
                                layout_matrix = rbind(c(1, 2), c(3, 4)))

# Write out to specific directory
ggsave("results/figures/test.pdf", prediction_plot, width = 14, height = 8)

