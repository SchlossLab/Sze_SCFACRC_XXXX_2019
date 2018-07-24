

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
    family = str_replace_all(family, "_", " "))

all_adn <- all_data %>% 
  filter(model == "adn_full") %>% 
  mutate(Variable = factor(Variable, 
                           levels = rev(Variable), 
                           labels = rev(Variable))) %>% 
  ggplot(aes(Variable, median_mda, color = family)) + 
  geom_pointrange(aes(ymin = iqr25, ymax = iqr75), position = position_dodge(width = 1), size = 1) + 
  coord_flip(ylim = c(0, 6)) + theme_bw() + 
  scale_color_manual(name = "", values = c('#FF3E96', '#8B475D', '#8B008B', '#63B8FF', '#008080', '#FF7F00')) + 
  theme(legend.position = "bottom")

#values = c('#9B30FF', '#63B8FF', '#008080', '#8FBC8F', '#FF7F00', '#6C7B8B')
#
#values = c("#FF6EB4", "#9B30FF", "#4B0082", "#4169E1", 
#           "#63B8FF", "#FFD700", "#00FF00", "#FF7F00")
  
only_otu_adn <- all_data %>% 
  filter(model == "adn_otu") %>% 
  mutate(Variable = factor(Variable, 
                           levels = rev(Variable), 
                           labels = rev(Variable))) %>% 
  ggplot(aes(Variable, median_mda, color = family)) + 
  geom_pointrange(aes(ymin = iqr25, ymax = iqr75), position = position_dodge(width = 1), size = 1) + 
  coord_flip(ylim = c(0, 6)) + theme_bw() + 
  scale_color_manual(name = "", values = c("#8B7E66", "#9B30FF", "#63B8FF", "#FF7F00")) + 
  theme(legend.position = "bottom")


all_crc <- all_data %>% 
  filter(model == "crc_full") %>% 
  mutate(Variable = factor(Variable, 
                           levels = rev(Variable), 
                           labels = rev(Variable))) %>% 
  ggplot(aes(Variable, median_mda, color = family)) + 
  geom_pointrange(aes(ymin = iqr25, ymax = iqr75), position = position_dodge(width = 1), size = 1) + 
  coord_flip(ylim = c(0, 5)) + theme_bw() + 
  theme(legend.position = "bottom")


only_otu_crc <- all_data %>% 
  filter(model == "crc_otu") %>% 
  mutate(Variable = factor(Variable, 
                           levels = rev(Variable), 
                           labels = rev(Variable))) %>% 
  ggplot(aes(Variable, median_mda, color = family)) + 
  geom_pointrange(aes(ymin = iqr25, ymax = iqr75), position = position_dodge(width = 1), size = 1) + 
  coord_flip(ylim = c(0, 5)) + theme_bw() + 
  theme(legend.position = "bottom")


### Create a merged graph
prediction_plot <- grid.arrange(only_otu_adn, all_adn, all_crc, only_otu_crc, 
                                layout_matrix = rbind(c(1, 2), c(3, 4)))

# Write out to specific directory
ggsave("results/figures/test.pdf", prediction_plot, width = 14, height = 8)

