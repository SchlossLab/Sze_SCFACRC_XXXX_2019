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
  
  tempData <- dataTable %>% 
    group_by(kegg_id) %>% 
    summarise(median_MDA = median(Overall), 
              quantile_25 = quantile(Overall, probs = 0.25), 
              quantile_75 = quantile(Overall, probs = 0.75)) %>% 
    arrange(desc(median_MDA))
  
  return(tempData)
}




##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

test <- group_stats(adn_model_imp)



