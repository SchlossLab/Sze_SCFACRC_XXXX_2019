### Graph Important Pathways and genes in picrust models
### Pull top pathways 
### Marc Sze



# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "scales"))


# Read in needed data
model_data <- list(
  adn = read_csv("data/process/tables/picrust_adn_top10percent_imp_model_summary.csv"), 
  crc = read_csv("data/process/tables/picrust_crc_top10percent_imp_model_summary.csv"))

# Create vector of models used
tumors <- c("adn", "crc")


##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################


# Function to generate pathway counts
get_scfa_pathway_counts <- function(tumor_type, dataList){
  
  tempData <- dataList[[tumor_type]] %>% 
    select(propionate_pathway, butyrate_pathway, acetate_pathway, isobutyrate_pathway, other_pathway) %>% 
    gather(pathway, occurence, everything()) %>% 
    group_by(pathway) %>% 
    summarise(total_count = sum(occurence)) %>% 
    mutate(tumor = tumor_type, 
           percent_total = total_count/sum(total_count))
  
  return(tempData)
}



  







##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

total_scfa_occurence_data <- sapply(tumors, 
                                    function(x) get_scfa_pathway_counts(x, model_data), simplify = F) %>% 
  bind_rows() %>% 
  mutate(pathway = factor(pathway, 
                          levels = c("acetate_pathway", "butyrate_pathway", "isobutyrate_pathway", 
                                     "propionate_pathway", "other_pathway"), 
                          labels = c("Acetate", "Butyrate", "Isobutyrate", 
                                     "Propionate", "Other")), 
         tumor = factor(tumor, 
                        levels = c("adn", "crc"), labels = c("Adenoma", "Carcinoma")))


# Generate the adenoma graph
picrust_top10_model_graph <- total_scfa_occurence_data %>%  
  ggplot(aes(tumor, percent_total, group = pathway, fill = pathway)) + 
  geom_bar(stat = "identity") + 
  labs(x = "", y = "Pathway Proportion in Top 10% of Model") + 
  theme_bw() + 
  scale_fill_manual(name = "Pathway", 
                    values = c('#000080', '#836FFF', '#00F5FF', '#1E90FF', '#808080')) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))
  
  
  
# Write out to specific directory
ggsave("results/figures/picrust_RF_model_pathway_proportions.tiff", picrust_top10_model_graph, width = 7, height = 6)















