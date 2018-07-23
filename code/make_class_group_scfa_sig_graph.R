### Graph the significant correlations and differences -- classification
### Group SCFA analysis
# Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))

taxonomy <- read_tsv("data/process/final.taxonomy") %>% 
  select(-Size) %>% 
  mutate(Taxonomy = str_replace_all(Taxonomy, "\\(\\d*\\)", "")) %>% 
  separate(Taxonomy, c("kingdom", "phyla", "class", "order", "family", "genus", "species"), sep = ";")

class_data <- read_csv("data/process/tables/significant_class_otu_comp_summary.csv") %>% 
  left_join(select(taxonomy, OTU, family), by = c("otu" = "OTU")) %>% 
  mutate(family = str_replace_all(family, "_unclassified", ""), 
         family = str_replace_all(family, "_", " "))

graph_data <- class_data %>% 
  mutate(family = case_when(
    family %in% c("Acidaminococcaceae", "Actinomycetaceae", "Bacteria",
                  "Pseudomonadaceae", "Streptococcaceae", "Sutterellaceae",
                  "Proteobacteria", "Firmicutes", "Bacteroidaceae",
                  "Clostridia", "Erysipelotrichaceae", "Oxalobacteraceae",
                  "Pasteurellaceae") ~ "Other",
    TRUE ~ family)) %>%
  group_by(scfa, dx, family) %>% 
  count()



number_class_plot <- graph_data %>% 
  ungroup() %>% 
  mutate(dx = factor(dx, 
                     levels = c("normal", "adenoma", "cancer"), 
                     labels = c("Normal", "Adenoma", "Carcinoma")), 
         family = factor(family, 
                         levels = c('Clostridiales', 'Lachnospiraceae', 'Porphyromonadaceae', 
                                    'Rikenellaceae', 'Ruminococcaceae', 'Other'), 
                         labels = c('Clostridiales', 'Lachnospiraceae', 'Porphyromonadaceae', 
                                    'Rikenellaceae', 'Ruminococcaceae', 'Other')), 
         scfa = factor(scfa, 
                       levels = c("acetate", "butyrate", "propionate"), 
                       labels = c("Acetate", "Butyrate", "Propionate"))) %>% 
  ggplot(aes(dx, n, fill = family)) + 
  geom_bar(stat = "identity") + theme_bw() + 
  # geom_bar(stat = "identity", position = position_dodge(1), show.legend = F) + 
  facet_grid(~scfa) +
  labs(x = "", y = "Number of Significant OTUs") + 
  scale_fill_manual(name = "", values = c('#9B30FF', '#63B8FF', '#008080', '#8FBC8F', '#FF7F00', '#6C7B8B')) + 
  theme(legend.position = "bottom") + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom", 
        legend.text = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))


ggsave("results/figures/FigureS2.pdf", number_class_plot, width = 8, height = 8)

