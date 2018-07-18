### Graph the aggregated data
### Visualize how well do the predictions do by group (normal, adenoma, cancer)?
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse"))

#setup variables that will be used
scfas <- c("acetate", "butyrate", "propionate")


# Read in the data
# The data is obtained from the model tracking of the run_rf_scfa_predictions.R script
# Put in a temporary run number until the analysis is re-run with it in the data frame.
model_data <- list(
  class_train = sapply(scfas, function(x) 
    read_csv(paste("data/process/tables/", x, "_classification_RF_train_probs_summary.csv", sep = "")), simplify = F), 
  class_test = sapply(scfas, function(x) 
    read_csv(paste("data/process/tables/", x, "_classification_RF_test_probs_summary.csv", sep = "")) %>% 
      mutate(run = rep(1:100, each = 85)), simplify = F), 
  reg_train = sapply(scfas, function(x) 
    read_csv(paste("data/process/tables/", x, "_regression_RF_train_conc_summary.csv", sep = "")), simplify = F),
  reg_test = sapply(scfas, function(x) 
    read_csv(paste("data/process/tables/", x, "_regression_RF_test_conc_summary.csv", sep = "")) %>% 
      mutate(run = rep(1:100, each = 85)), simplify = F),
  log_reg_train = sapply(scfas, function(x) 
    read_csv(paste("data/process/tables/log_", x, "_regression_RF_train_conc_summary.csv", sep = "")), simplify = F),
  log_reg_test = sapply(scfas, function(x) 
    read_csv(paste("data/process/tables/log_", x, "_regression_RF_test_conc_summary.csv", sep = "")) %>% 
      mutate(run = rep(1:100, each = 85)), simplify = F))


# Create summary for chi-square test for proportion of correct and incorrect samples
create_count_table <- function(datatable, train_data = T){
  # datatable is each respective data frame within the model data list
  # train_data denotes whether information was from training (T) or testing (F)
  
  # creates the class column based on train or test
  if(train_data == T){
    
    tempData <- datatable %>% 
      mutate(correct_class = case_when(
        pred == obs ~ "yes", 
        TRUE ~ "no")) 
  } else{
    
    tempData <- datatable %>% 
      mutate(correct_class = case_when(
        tempPredictions == high_low ~ "yes", 
        TRUE ~ "no")) 
  }
  # Uses the class column and creates yes and no columns
  tempData <- tempData %>% 
    group_by(run, dx) %>% 
    summarise(yes = table(correct_class)[2], 
              no = table(correct_class)[1]) %>% 
    ungroup() %>% 
    group_by(dx) %>% 
    summarise(yes = round(mean(yes, na.rm = T)), 
              no = round(mean(no, na.rm = T))) %>% 
    as.data.frame()
  # Creates rownames from dx 
  rownames(tempData) <- tempData$dx
  # removes the dx column resulting in a confusion matrix that can be used with chi square analysis
  tempData <- as.matrix(tempData[, -1])
  
  # Return this to the local environment
  return(tempData)
  
}

# Create needed data for the analysis of the regression models
get_reg_diff <- function(datatable, train_data = T){
  # datatable is each respective data frame within the model data list
  # train_data denotes whether information was from training (T) or testing (F)
  
  # creates the difference column based on train or test
  if(train_data == T){
    
    tempData <- datatable %>% 
      mutate(difference = obs - pred)
    
  } else{
    
    tempData <- datatable %>% 
      mutate(difference = mmol_kg - tempPredictions)
  }
  # generates summary data for each run 
  tempData <- tempData %>% 
    group_by(run, dx) %>% 
    summarise(median_diff = median(difference), 
              min_dff = min(difference), 
              max_diff = max(difference)) %>% 
    ungroup()
  # returns the summary table to the local environment
  return(tempData)
}





