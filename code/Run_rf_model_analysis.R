### Test with t-test if AUC distributions are different
### With and without scfas added
### Marc Sze



# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse"))

# Load in summary data tables that will be used
full_adn_summary <- read_csv("data/process/tables/adn_full_AUC_model_summary.csv")
full_crc_summary <- read_csv("data/process/tables/crc_full_AUC_model_summary.csv")

otu_adn_summary <- read_csv("data/process/tables/adn_otu_only_AUC_model_summary.csv")
otu_crc_summary <- read_csv("data/process/tables/crc_otu_only_AUC_model_summary.csv")



##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################

# Function to create the output summary data
make_comparison <- function(data_vector1, data_vector2){
  
  result_vector <- c(
    pvalue = t.test(data_vector1, data_vector2, 
                     alternative = "two.sided")$p.value, 
    avg_vec1 = mean(data_vector1, na.rm = T), 
    avg_vec2 = mean(data_vector2, na.rm = T), 
    sd_vec1 = sd(data_vector1, na.rm = T), 
    sd_vec2 = sd(data_vector2, na.rm = T))
  
  return(result_vector)
  
  
}



##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

summary_results <- rbind(make_comparison(full_adn_summary$test_auc, otu_adn_summary$test_auc), 
              make_comparison(full_crc_summary$test_auc, otu_crc_summary$test_auc)) %>% 
  as.data.frame() %>% 
  mutate(bh = p.adjust(pvalue, method = "BH"), 
         model = c("adn", "crc")) %>% 
  rename(avg_full = avg_vec1, avg_otu = avg_vec2, sd_full = sd_vec1, sd_otu = sd_vec2) %>% 
  select(model, avg_full, sd_full, avg_otu, sd_otu, pvalue, bh)
  








