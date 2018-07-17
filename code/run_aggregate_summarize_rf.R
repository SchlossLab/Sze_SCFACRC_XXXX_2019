### Aggregate and summarize the model data
### How well do the predictions do by group (normal, adenoma, cancer)?
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse"))

#setup variables that will be used
scfas <- c("acetate", "butyrate", "isobutyrate", "propionate")


# Read in the data
model_data <- list(
  class_train = sapply(scfas, function(x) 
    read_csv(paste("data/process/tables/", x, "_classification_RF_train_probs_summary.csv", sep = "")), simplify = F), 
  class_test = sapply(scfas, function(x) 
    read_csv(paste("data/process/tables/", x, "_classification_RF_test_probs_summary.csv", sep = "")), simplify = F), 
  reg_train = sapply(scfas, function(x) 
    read_csv(paste("data/process/tables/", x, "_regression_RF_train_conc_summary.csv", sep = "")), simplify = F),
  reg_test = sapply(scfas, function(x) 
    read_csv(paste("data/process/tables/", x, "_regression_RF_test_conc_summary.csv", sep = "")), simplify = F),
  log_reg_train = sapply(scfas, function(x) 
    read_csv(paste("data/process/tables/log_", x, "_regression_RF_train_conc_summary.csv", sep = "")), simplify = F),
  log_reg_test = sapply(scfas, function(x) 
    read_csv(paste("data/process/tables/log_", x, "_regression_RF_test_conc_summary.csv", sep = "")), simplify = F))


# Create summary for chi-square test for proportion of correct and incorrect samples
create_count_table <- function(datatable, train_data = T){
  
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
  
  tempData <- tempData %>% 
    group_by(run, dx) %>% 
    summarise(yes = table(correct_class)[2], 
              no = table(correct_class)[1]) %>% 
    ungroup() %>% 
    group_by(dx) %>% 
    summarise(yes = round(mean(yes)), 
              no = round(mean(no))) %>% 
    as.data.frame()
  
  rownames(tempData) <- tempData$dx

  tempData <- as.matrix(tempData[, -1])
  
  
  return(tempData)
  
}

# Create needed data for the analysis
get_reg_diff <- function(datatable, train_data = T){
  
  if(train_data == T){
    
    tempData <- datatable %>% 
      mutate(difference = obs - pred)
    
  } else{
    
    tempData <- datatable %>% 
      mutate(difference = mmol_kg - tempPredictions)
  }
  
  tempData <- tempData %>% 
    group_by(run, dx) %>% 
    summarise(median_diff = median(difference), 
              min_dff = min(difference), 
              max_diff = max(difference)) %>% 
    ungroup()
  
  return(tempData)
}


# Set up the names of the list
data_sets_names <- names(model_data)

# Set up empty list to store the information
test_results <- list(
  class_train = c(), class_test = c(), reg_train = c(), reg_test = c(), 
  log_reg_train = c(), log_reg_test = c())


test <- model_data$reg_test$acetate %>% 
  mutate(difference = pred - obs) %>% 
  group_by(run, dx) %>% 
  summarise(median_diff = median(difference), 
            min_dff = min(difference), 
            max_diff = max(difference)) %>% 
  ungroup()


kruskal_summary <- test %>% 
  nest() %>% 
  mutate(k_test = map(data, ~kruskal.test(median_diff ~ factor(dx), data = .x)), 
         k_summary = map(k_test, broom::tidy)) %>% 
  select(k_summary) %>% 
  unnest()


dunn_summary <- as.data.frame.list(dunn.test::dunn.test(test$median_diff, factor(test$dx), method = "bh"))
  


for(i in data_sets_names){
  
  if(str_detect(i, "train")){
    
    model_type = T
    
  } else{
    
    model_type = F
  }
  
  if(str_detect(i, "reg")){
    
    temp_tables <- map(model_data[[i]], function(x) get_reg_diff(x, train_data = model_type))
    
    analysis_summary <- map(temp_tables, function(x)
      x %>% nest() %>%
        mutate(k_test = map(data, ~kruskal.test(median_diff ~ factor(dx), data = .x)),
               k_summary = map(k_test, broom::tidy)) %>%
        select(k_summary) %>%
        unnest())
    
    for(j in names(analysis_summary)){
      
      temp_pvalue <- analysis_summary[[j]] %>% 
        select(p.value) %>% pull()
      
      str(temp_pvalue)
      
      if(temp_pvalue < 0.05){
        
        analysis_summary[[j]] <- as.data.frame.list(
          dunn.test::dunn.test(temp_tables[[j]]$median_diff, factor(temp_tables[[j]]$dx), method = "bh"))
      }
    }
      
    test_results[[i]] <- analysis_summary

      
  } else{
    
    temp_tables <- map(model_data[[i]], function(x) create_count_table(x, train_data = model_type))
    
    test_results[[i]] <- map(temp_tables, function(x) c(statistic = unname(chisq.test(x)$statistic), 
                                                        pvalue = chisq.test((x))$p.value))
    
  }
  
  
  
  print(paste("Finished analysis of ", i, "...", sep = ""))
}








# Assess probability differences


