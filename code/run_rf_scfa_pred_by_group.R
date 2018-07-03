### Investigate how important OTUs based on MDA change for each group
### Using RF models do OTUs change if in non-tumor, adenoma, or carcinoma samples
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "caret", "pROC", "randomForest"))

#setup variables that will be used
scfas <- c("acetate", "butyrate", "isobutyrate", "propionate")

# Samples that had to little volume to actually be measured
not_measured <- c("3367653", "2041650", "2043650", "3027650")

##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################

# Function to load in the respective scfa data
upload_scfa_data <- function(scfa_name, path_to_file, ending){
  # scfa_name is a variable that stores the scfa of interest e.g. "acetate"
  # path_to_file stores the directory path where the files are
  # ending is the common ending of the file of interest
  
  # the command that reads in the files
  tempData <- read_csv(paste(path_to_file, scfa_name, ending, sep = "")) %>% 
    rename(Group = study_id)
  # sends the data to the global work environment
  return(tempData)
}

# Function to remove specific data from sets of interest
pare_down_data <- function(ids_to_remove, scfa_list, otu_data, meta_data){
  
  # Remove data from tables first
  temp_otuData <- otu_data %>% filter(!(Group %in% ids_to_remove))
  
  nzv <- nearZeroVar(temp_otuData)
  
  nzv_rm_tempOTU <- temp_otuData[, -nzv]
  
  temp_meta_data <- meta_data %>% filter(!(Group %in% ids_to_remove))
  
  # Remove data from lists
  tempList <- lapply(scfa_list, 
                     function(x) filter(x, !(Group %in% ids_to_remove)) %>% 
                       inner_join(temp_meta_data, by = "Group") %>% 
                       inner_join(nzv_rm_tempOTU, by = "Group")) 
  
  return(tempList)
}

# Function to generate the high low groupings for each respective SCFA
get_high_low <- function(all_data_list, scfa_median_list, scfas){
  
  tempList <- sapply(scfas, 
                     function(x) all_data_list[[x]] %>% 
                       mutate(high_low = ifelse(mmol_kg <= scfa_median_list[[x]], 
                                                invisible("low"), invisible("high"))) %>% 
                       select(Group, mmol_kg, high_low, Dx_Bin, dx, fit_result, everything()), 
                     simplify = F)
  
  return(tempList)
  
}

# Function to create relevent RF data frames
split_dataframe <- function(i, dataList){
  
  tempData <- dataList[[i]]
  
  temp_regression <- tempData %>% select(mmol_kg, contains("Otu"))
  
  temp_high_low <- tempData %>% select(high_low, contains("Otu"))
  
  finalList <- list(rf_regression = temp_regression, 
                    rf_groups = temp_high_low)
  return(finalList)
}


# Function to make an initial 80/20 split within the data
eighty_twenty_split <- function(i, dataList){
  
  tempData <- dataList[[i]]
  
  sampling_vector <- c(1:length(rownames(tempData)))
  
  trainValues <- sample(sampling_vector, round(length(sampling_vector)*0.8))
  
  #tempTrain <- tempData %>% sample_frac(0.8, replace = FALSE)
  
  #tempTest <- tempData %>% filter(!(rownames(tempTrain) %in% rownames(.)))
  
  finalList <- list(
    training_data = tempData[trainValues, ], 
    test_data = tempData[-trainValues, ])
  
  
  return(finalList)
  
}


# Function that will run and create the needed models
make_rf_model <- function(i, run_marker, train_data_name, 
                          cat_column, method_to_use, 
                          metric_to_use, dataList){
  # i is the scfa of interest
  # run_marker is the model iteration that has currently completed
  # train_data_name is a vector for the data set to be used
  # dataList is the data table to be used for model training
  # cat_column is the grouping column to be used
  # method_to_use are "rf" for classification and "rf" for regression rf
  # metric_to_use are "ROC" for classification and "Rsquared" for regression rf
  
  tempData <- dataList[[i]][[train_data_name]] %>% 
    select(high_low, contains("Otu"))
  
  #Create Overall specifications for model tuning
  # number controls fold of cross validation
  # Repeats control the number of times to run it
  if(metric_to_use == "ROC"){
    
    fitControl <- trainControl(## 10-fold CV
      method = "cv",
      number = 10,
      p = 0.8, 
      classProbs = TRUE, 
      summaryFunction = twoClassSummary, 
      savePredictions = "final")
  } else{
    
    fitControl <- trainControl(## 10-fold CV
      method = "cv",
      number = 10,
      p = 0.8, 
      savePredictions = "final")
  }
  
  
  #Train the model
  training_model <- 
    train(formula(paste(cat_column, " ~ .", sep = "")), data = tempData, 
          method = method_to_use, 
          ntree = 500, 
          trControl = fitControl,
          metric = metric_to_use,
          importance = TRUE, 
          verbose = FALSE)
  
  #Print out tracking message
  print(paste("Completed ", run_marker, " RF model for ", 
              i, " ", cat_column,  " using cv", sep = ""))
  
  # Return the model object
  return(training_model)
}

# Create a function to run the prediction on the test data
run_prediction <- function(i, model_data, train_data_name, control_type, 
                           var_of_int, dataList){
  
  tempData <- dataList[[i]][[train_data_name]] %>% 
    select(high_low, contains("Otu"))
  
  tempModel <- model_data[[i]]
  
  
  
  tempPredictions <- predict(tempModel, tempData, type = control_type)
  
  tempPredictions <- cbind(tempPredictions, actual = tempData[, var_of_int])
  
  
  return(tempPredictions)
}


# Generate ROC curves for the test data
get_test_roc <- function(i, var_of_int, pred_name, pred_list, classif = T){
  
  if(classif == T){
    
    actual_values <- factor(pred_list[[i]][, var_of_int])
    pred_values <- factor(pred_list[[i]][, pred_name], ordered = TRUE)
    
  } else{
    
    actual_values <- pred_list[[i]][, var_of_int]
    pred_values <- pred_list[[i]][, pred_name]
  }
  
  tempPred_ROC <- roc(actual_values ~ pred_values)
  
  return(tempPred_ROC)
}


# Function to aggregate test data
add_model_summary_data <- function(i, model_list, test_data, data_Table, classif = T){
  
  tempTrainResults <- model_list[[i]]$results
  
  if(classif == T){
    
    tempROC_test <- test_data[[i]]
    
    tempdata <- tempTrainResults %>% 
      filter(ROC == max(ROC)) %>% 
      rename(train_AUC = ROC) %>% 
      mutate(test_AUC = tempROC_test$auc[1]) %>% 
      select(mtry, train_AUC, test_AUC, everything())
    
    data_Table[[i]] <- data_Table[[i]] %>% bind_rows(tempdata)
    
  } else{
    
    temp_test <- summary(lm(test_data[[i]]$mmol_kg ~ test_data[[i]]$tempPredictions))
    
    tempdata <- tempTrainResults %>% 
      filter(Rsquared == max(Rsquared)) %>% 
      rename(train_r2 = Rsquared) %>% 
      mutate(test_r2 = temp_test$adj.r.squared) %>% 
      select(mtry, train_r2, test_r2, everything())
    
    data_Table[[i]] <- data_Table[[i]] %>% bind_rows(tempdata)
  }
  
  return(data_Table[[i]])
}

#Function to store the importance information for each run
grab_importance <- function(i, run_number, type_of_model, modelList, dataHolder){
  
  tempModel <- modelList[[i]]
  
  if(type_of_model == "classification"){
    
    test <- as.data.frame.list(varImp(tempModel, scale = FALSE)) %>% 
      mutate(otu = rownames(.), Overall = importance.high, 
             run = run_number) %>% 
      select(-model, -calledFrom, -importance.high, -importance.low) %>% 
      select(Overall, otu, run)
  } else{
    
    test <- as.data.frame.list(varImp(tempModel, scale = FALSE)) %>% 
      mutate(otu = rownames(.), run = run_number) %>% 
      select(-model, -calledFrom)
  }
  
  dataHolder[[i]] <- dataHolder[[i]] %>% bind_rows(test)
  
  return(dataHolder[[i]])
}



##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

# Load in needed metadata and rename sample column
picrust_meta <- read_tsv("data/process/picrust_metadata") %>% 
  mutate(Group = as.character(Group))

# Load in more specific follow up data
metaF <- read_csv("data/raw/metadata/good_metaf_final.csv")

# Load in shared data
shared <- read_tsv("data/process/final.0.03.subsample.shared") %>% 
  mutate(Group = as.character(Group)) %>% 
  select(-label, -numOtus) 


# Generate IDs that are initial and follow up
initial_IDs <- as.character(metaF$initial)
follow_IDs <- as.character(metaF$followUp)

# read in all the respective scfa data tables
scfa_data <- sapply(scfas, 
                    function(x) upload_scfa_data(x, "data/process/tables/", "_final_data.csv"), 
                    simplify = F)

# Combine all data into a single list of lists and remove follow samples
all_data <- pare_down_data(c(initial_IDs, follow_IDs), scfa_data, shared, picrust_meta)


# Get medians for each scfa measured
scfa_medians <- lapply(all_data, function(x) median(x$mmol_kg))

# get final data tables to be used in the RF analysis
final_data <- get_high_low(all_data, scfa_medians, scfas)


# Split into normal, adenoma, and cancer groups
split_final_data <- list(
  final_data_normal = sapply(scfas, 
                             function(x) final_data[[x]] %>% 
                               filter(dx == "normal"), simplify = F), 
  final_data_adenoma = sapply(scfas, 
                              function(x) final_data[[x]] %>% 
                                filter(dx == "adenoma"), simplify = F), 
  final_data_carcinoma = sapply(scfas, 
                                function(x) final_data[[x]] %>% 
                                  filter(dx == "cancer"), simplify = F))

# Set up initial store for summary data
class_summary_model_data <- list(
  final_data_normal = list(acetate = data_frame(), butyrate = data_frame(), 
                           isobutyrate = data_frame(), propionate = data_frame()), 
  final_data_adenoma = list(acetate = data_frame(), butyrate = data_frame(), 
                            isobutyrate = data_frame(), propionate = data_frame()), 
    final_data_carcinoma = list(acetate = data_frame(), butyrate = data_frame(), 
                              isobutyrate = data_frame(), propionate = data_frame()))
  
  
class_important_vars <- list(
  final_data_normal = list(acetate = data_frame(), butyrate = data_frame(), 
                           isobutyrate = data_frame(), propionate = data_frame()), 
  final_data_adenoma = list(acetate = data_frame(), butyrate = data_frame(), 
                            isobutyrate = data_frame(), propionate = data_frame()), 
  final_data_carcinoma = list(acetate = data_frame(), butyrate = data_frame(), 
                              isobutyrate = data_frame(), propionate = data_frame()))


# Run a 100 different splits for each group and SCFA on classification
for(i in names(split_final_data)){
  
  print(paste("Getting model data from ", i, " only....", sep = ""))
  
  for(j in 1:1){
    
    # Generate an 80/20 data split
    rf_data <- sapply(scfas, function(x) eighty_twenty_split(x, split_final_data[[i]]), simplify = F)
    
    class_test <- sapply(scfas, function(x)
      make_rf_model(x, j, "training_data", "high_low", "rf", "ROC", rf_data),
      simplify = F)

    class_pred <- sapply(scfas, function(x)
      run_prediction(x, class_test, "test_data", "raw", "high_low", rf_data),
      simplify = F)

    ROC_data <- sapply(scfas, function(x)
      get_test_roc(x, "high_low", "tempPredictions", class_pred, classif = T), simplify = F)

    class_summary_model_data[[i]] <- sapply(scfas, function(x)
      add_model_summary_data(x, class_test, ROC_data, class_summary_model_data[[i]], classif = T), simplify = F)

    class_important_vars[[i]] <- sapply(scfas, function(x)
      grab_importance(x, j, "classification", class_test, class_important_vars[[i]]), simplify = F)
  }
  
  print(paste("Completed getting all model data for ", i, " only.", sep = ""))
}






