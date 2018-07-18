### Aggregate and summarize the model data
### How well do the predictions do by group (normal, adenoma, cancer)?
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


# Assign SCFA label to the respective comparison
add_scfa_label <- function(dataList, list_names){
  
  for(i in list_names){
    
    temp_list <- dataList[[i]]
    
    scfa_names <- names(temp_list)
    
    if(str_detect(i, "class")){
      
      modified_list <- sapply(scfa_names, function(x) c(temp_list[[x]], scfa = x), simplify = F)
      
    } else{
      
      modified_list <- sapply(scfa_names, function(x) temp_list[[x]] %>% 
                                tbl_df() %>% mutate(scfa = x), simplify = F)
    }
    
    
    dataList[[i]] <- modified_list
  }
  
  return(dataList)
  
}




# Set up the names of the list
data_sets_names <- names(model_data)

# Set up empty list to store the information
test_results <- list(
  class_train = c(), class_test = c(), reg_train = c(), reg_test = c(), 
  log_reg_train = c(), log_reg_test = c())

# Initial for loop to run through the analysis for each model train and test (everything in model_data)
for(i in data_sets_names){
  # Set the model type for accurate analysis
  if(str_detect(i, "train")){
    
    model_type = T
    
  } else{
    
    model_type = F
  }
  # differentiate whether it is a regression or classification model
  if(str_detect(i, "reg")){
    # creates the regression data frames needed
    temp_tables <- map(model_data[[i]], function(x) get_reg_diff(x, train_data = model_type))
    # Runs the initial kruskal wallis test
    analysis_summary <- map(temp_tables, function(x)
      x %>% nest() %>%
        mutate(k_test = map(data, ~kruskal.test(median_diff ~ factor(dx), data = .x)),
               k_summary = map(k_test, broom::tidy)) %>%
        select(k_summary) %>%
        unnest())
    # Cycles through each of the kruskal wallis test run
    for(j in names(analysis_summary)){
      # Pull the p-value for that specific SCFA
      temp_pvalue <- analysis_summary[[j]] %>% 
        select(p.value) %>% pull()
      # Check if it is below 0.05
      if(temp_pvalue < 0.05){
        # Run a dunn's analysis and replce the kruskal results with the dunns results
        analysis_summary[[j]] <- as.data.frame.list(
          dunn.test::dunn.test(temp_tables[[j]]$median_diff, factor(temp_tables[[j]]$dx), method = "bh")) %>% 
          mutate(model = i)
      }
    }
    # Add it to the overall storage list
    test_results[[i]] <- analysis_summary

      
  } else{
    # create data needed for the classification analysis
    temp_tables <- map(model_data[[i]], function(x) create_count_table(x, train_data = model_type))
    # Run a chi square test and add it to the overall storage list
    test_results[[i]] <- map(temp_tables, function(x) c(statistic = unname(chisq.test(x)$statistic), 
                                                        pvalue = chisq.test((x))$p.value, 
                                                        model = i))
    
  }
  
  # Keeps track of progress on stdout
  print(paste("Finished analysis of ", i, "...", sep = ""))
}

test <- add_scfa_label(test_results, data_sets_names)







# Assess probability differences


