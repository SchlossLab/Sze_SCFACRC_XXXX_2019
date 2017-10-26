### Identify Important Pathways in picrust model
### Pull top pathways 
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse"))

# Load in needed kegg pathway data
kegg_pathway <- read_csv("data/process/tables/kegg_id_key.csv")

# Load in model information
adn_model_imp <- read_csv("data/process/tables/adenoma_imp_otus_classification_RF_summary.csv")
crc_model_imp <- read_csv("data/process/tables/cancer_imp_otus_classification_RF_summary.csv")


##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################

# Function to summarize all the Keggs by median and quantiles and order highest to lowest
group_stats <- function(dataTable){
  
  tempData <- get(dataTable)
  
  tempData <- tempData %>% 
    group_by(kegg_id) %>% 
    summarise(median_MDA = median(Overall), 
              quantile_25 = quantile(Overall, probs = 0.25), 
              quantile_75 = quantile(Overall, probs = 0.75)) %>% 
    arrange(desc(median_MDA))
  
  return(tempData)
}


# Function to ID the top X pathways
top_model_pathways <- function(dataTable, mapFile, topNumber){
  
  tempData <- dataTable %>% 
    slice(1:topNumber) %>% 
    left_join(mapFile, by = "kegg_id")
  
  return(tempData)
}





##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

# Get the summary imp data
summary_imp_lists <- sapply(c("adn_model_imp", "crc_model_imp"), 
                            function(x) group_stats(x), simplify = F)
  
top_ten_pathways <- sapply(c("adn_model_imp", "crc_model_imp"), 
               function(x) top_model_pathways(summary_imp_lists[[x]], kegg_pathway, 10), simplify = F)
  
