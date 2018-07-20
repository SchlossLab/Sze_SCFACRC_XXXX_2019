### Graph the significant correlations and differences 
### Group SCFA analysis
# Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))


reg_data <- read_csv("data/process/tables/significant_reg_otu_comp_summary.csv") %>% 
  mutate(genus = str_replace_all(genus, "_unclassified", ""), 
         genus = str_replace_all(genus, "_", " "))


graph_data <- reg_data %>% 
  group_by(scfa, dx, genus) %>% 
  count()

graph_data %>% 
  ungroup() %>% 
  mutate(dx = factor(dx, 
                     levels = c("normal", "adenoma", "cancer"), 
                     labels = c("Normal", "Adenoma", "Carcinoma"))) %>% 
  ggplot(aes(dx, n, group = genus, fill = genus)) + 
  geom_bar(stat = "identity") + 
 # geom_bar(stat = "identity", position = position_dodge(1), show.legend = F) + 
  facet_grid(~scfa) + 
  labs(x = "", y = "Number of Significant OTUs") + 
  scale_fill_viridis(discrete = T, option = "magma") + 
  theme(legend.position = "bottom")

