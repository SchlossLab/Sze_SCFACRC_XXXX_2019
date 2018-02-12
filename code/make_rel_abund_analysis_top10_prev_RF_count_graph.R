### Graph looking at number of OTUs with P-value < 0.05 for each respective top 10 genera
### from previous RF models
### Marc Sze



# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse"))


lowest_txID <- list(adn = c("Lachnospiraceae", "Clostridiales", "Blautia", "Bacteroides", 
                            "Odoribacter", "Ruminococcaceae", "Flavonifractor", 
                            "Roseburia", "Escherichia/Shigella", "Clostridium_XIVa"), 
                    adv_adn = c("Clostridiales", "Clostridium_XIVa", "Roseburia", "Odoribacter", "Bacteroides", 
                                "Ruminococcaceae", "Streptococcus", "Blautia", "Lachnospiraceae", "Ruminococcus"), 
                    crc = c("Porphyromonas", "Parvimonas", "Fusobacterium", "Gemella", "Prevotella", "Streptococcus", 
                            "Coprobacillus", "Pasteurellaceae", "Collinsella", "Bilophila"))


# generate vector to make sapply run
tumors <- c("adn", "adv_adn", "crc")

data_tables <- sapply(tumors, 
                      function(x) read_csv(paste("data/process/tables/", x, 
                                                 "_16S_top10_RF_taxa_testing.csv", sep = "")), simplify = F)


# Filter based on pvalues
pvalue_filtered_tables <- sapply(tumors, 
                                 function(x) filter(data_tables[[x]], pvalue < 0.05), simplify = F)

##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################

# Function to tally up the number of OTUs with pvalues less than 0.05
get_pvalue_counts <- function(tumor_type, dataList){
  
  tempData <- dataList[[tumor_type]] %>% 
    mutate(temp_count = 1) %>% 
    group_by(genus) %>% 
    summarise(total_otus = sum(temp_count)) %>% 
    arrange(total_otus)
  
  return(tempData)
  
}





##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

# Generate the counts
pvalue_taxa_counts <- sapply(tumors, 
                             function(x) get_pvalue_counts(x, pvalue_filtered_tables), simplify = F)

# Generate the labels
labels_to_use <- sapply(tumors, 
                        function(x) pvalue_taxa_counts[[x]]$genus, simplify = F)

# Generate the adenoma graph
adenoma_graph <- pvalue_taxa_counts[["adn"]] %>% 
  mutate(genus = factor(genus, levels = labels_to_use[["adn"]], 
    labels = labels_to_use[["adn"]])) %>% 
  ggplot(aes(genus, total_otus)) + 
  geom_bar(stat = "identity") + coord_flip(ylim = c(0, 10)) + 
  labs(x = "Lowest Taxa ID", y = "OTUs with a P-value < 0.05") + 
  theme_bw() + 
  annotate("text", label = paste("Control vs. Adenoma"), x = 0.75, y = 8, size = 4) + 
  ggtitle("A") + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))

# Generate the advanced adenoma graph
adv_adenoma_graph <- pvalue_taxa_counts[["adv_adn"]] %>% 
  mutate(genus = factor(genus, levels = labels_to_use[["adv_adn"]], 
                        labels = labels_to_use[["adv_adn"]])) %>% 
  ggplot(aes(genus, total_otus)) + 
  geom_bar(stat = "identity") + coord_flip(ylim = c(0, 10)) + 
  labs(x = "Lowest Taxa ID", y = "OTUs with a P-value < 0.05") + 
  theme_bw() + 
  annotate("text", label = paste("Control vs. Adv. Adenoma"), x = 1, y = 7.5, size = 4) + 
  ggtitle("B") + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))

# Generate the carcinoma graph
carcinoma_graph <- pvalue_taxa_counts[["crc"]] %>% 
  mutate(genus = factor(genus, levels = labels_to_use[["crc"]], 
                        labels = labels_to_use[["crc"]])) %>% 
  ggplot(aes(genus, total_otus)) + 
  geom_bar(stat = "identity") + coord_flip(ylim = c(0, 10)) + 
  labs(x = "Lowest Taxa ID", y = "OTUs with a P-value < 0.05") + 
  theme_bw() + 
  annotate("text", label = paste("Control vs. Carcinoma"), x = 0.85, y = 8, size = 4) + 
  ggtitle("C") + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))


# Combine the graphs together

sig_occurence_plot <- grid.arrange(adenoma_graph, adv_adenoma_graph, carcinoma_graph, nrow = 1, ncol = 3)

# Write out to specific directory
ggsave("results/figures/top10_RF_sig_occurence.tiff", sig_occurence_plot, width = 15, height = 6)






