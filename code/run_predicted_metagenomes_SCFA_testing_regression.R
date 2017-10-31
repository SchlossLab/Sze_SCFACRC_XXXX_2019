### Investigate how well the picrust predicted metagenomes predict SCFA levels
### compare most differentially expressed genes by disease
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "biomformat", "caret", "randomForest", "pROC"))

# Load in needed metadata and rename sample column
picrust_meta <- read_tsv("data/process/picrust_metadata") %>% 
  mutate(Group = as.character(Group))

# Load in more specific follow up data
metaF <- read_csv("data/raw/metadata/good_metaf_final.csv")

# Load in SCFA data
scfa_data <- sapply(c("acetate", "butyrate", "isobutyrate", "propionate"), 
                    function(x) read_csv(paste("data/process/tables/", x, "_final_data.csv", sep = "")), 
                    simplify = F)


# Load in predicted metagenome data
all_biom_data <- read_biom("data/process/predicted_metagenomes.biom")

# Converts the needed biom data into a matrix
full_pi_data <- t(as(biom_data(all_biom_data), "matrix")) %>% 
  as.data.frame() %>% 
  mutate(sample_id = rownames(.)) %>% 
  select(sample_id, everything())

# Generate the near zero variance
nzv <- nearZeroVar(full_pi_data)

# Remove these from the data set
nzv_full_pi_data <- full_pi_data[, -nzv]

##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################

# Function to create data groups to be used in the RF model
make_rf_data <- function(scfa_name, dataTable, scfa_data, 
                         exclude_samples, column_of_int){
  
  tempData <- scfa_data[[scfa_name]] %>% 
    filter(!(study_id %in% exclude_samples)) %>% 
    select(study_id, contains(column_of_int)) %>% 
    inner_join(dataTable, by = c("study_id" = "sample_id")) %>% 
    select(-study_id)
  
  return(tempData)
}


# Function to generate the high low groupings for each respective SCFA
get_high_low <- function(all_data_list, scfa_median_list, scfas){
  
  tempList <- sapply(scfas, 
                     function(x) all_data_list[[x]] %>% 
                       mutate(high_low = ifelse(mmol_kg <= scfa_median_list[[x]], 
                                                invisible("low"), invisible("high"))) %>% 
                       select(study_id, mmol_kg, high_low, everything()), 
                     simplify = F)
  
  return(tempList)
  
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
  
  tempData <- dataList[[i]][[train_data_name]]
  
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



#Function to store the importance information for each run
grab_importance <- function(i, run_number, modelList, dataHolder){
  
  tempModel <- modelList[[i]]
  
  if (i ==  "adenoma"){
    
    test <- as.data.frame.list(varImp(tempModel, scale = FALSE)) %>% 
      mutate(kegg_id = rownames(.), Overall = importance.adenoma, 
             run = run_number)
  } else{
    
    test <- as.data.frame.list(varImp(tempModel, scale = FALSE)) %>% 
      mutate(kegg_id = rownames(.), Overall = importance.cancer, 
             run = run_number)
  }
  
  dataHolder[[i]] <- dataHolder[[i]] %>% bind_rows((test %>% select(Overall, kegg_id, run)))
  

  return(dataHolder[[i]])
}

# Create a function to run the prediction on the test data
run_prediction <- function(i, model_data, train_data_name, control_type, 
                           var_of_int, dataList){
  
  tempData <- dataList[[i]][[train_data_name]]
  tempModel <- model_data[[i]]
  
  tempPredictions <- predict(tempModel, tempData, type = control_type)
  
  tempPredictions <- cbind(tempPredictions, actual = tempData[, var_of_int])
  
  
  return(tempPredictions)
}


# Generate ROC curves for the test data
get_test_roc <- function(i, var_of_int, pred_name, pred_list){
  
  actual_values <- pred_list[[i]][, var_of_int]
  pred_values <- factor(pred_list[[i]][, pred_name], ordered = TRUE)
  
  tempPred_ROC <- roc(actual_values ~ pred_values)
  
  return(tempPred_ROC)
}


# Function to aggregate test data
add_model_summary_data <- function(i, model_list, test_data, data_Table){
  
  tempTrainResults <- model_list[[i]]$results
  
  tempROC_test <- test_data[[i]]
  
  tempdata <- tempTrainResults %>% 
    filter(ROC == max(ROC)) %>% 
    rename(train_AUC = ROC) %>% 
    mutate(test_AUC = tempROC_test$auc[1]) %>% 
    select(mtry, train_AUC, test_AUC, everything())
  
  data_Table[[i]] <- data_Table[[i]] %>% bind_rows(tempdata)
  
  return(data_Table[[i]])
}


##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

# Set up initial store for summary data
class_important_vars <- list(acetate = data_frame(), butyrate = data_frame(), 
                             isobutyrate = data_frame(), propionate = data_frame())

class_summary_data <- list(acetate = data_frame(), butyrate = data_frame(), 
                           isobutyrate = data_frame(), propionate = data_frame())

regress_important_vars <- list(acetate = data_frame(), butyrate = data_frame(), 
                               isobutyrate = data_frame(), propionate = data_frame())

regress_summary_data <- list(acetate = data_frame(), butyrate = data_frame(), 
                             isobutyrate = data_frame(), propionate = data_frame())

# Set up model type variables
scfas <- c("acetate", "butyrate", "isobutyrate", "propionate")

# Get scfas medians
scfa_median_vector <- sapply(scfas, function(x) median(scfa_data[[x]]$mmol_kg))

# Generate High/Low variable vector
all_scfa_data <- get_high_low(scfa_data, scfa_median_vector, scfas)

# Create RF data tables to be used
combined_data_amount <- sapply(scfas, 
                  function(x) make_rf_data(x, nzv_full_pi_data, all_scfa_data, 
                                           c(metaF$initial, metaF$followUp), "mmol_kg"), simplify = F)


# run through a 100 different 80/20 splits for modeling
for(j in 1:100){
  
  # Create test and train data for both conditions
  final_rf_data <- sapply(scfas, 
                          function(x) eighty_twenty_split(x, combined_data_hl), simplify = F)
  
  # Create the two models
  rf_models <- sapply(scfas, 
                      function(x) make_rf_model(x, j, "training_data", "disease", "rf", "ROC", final_rf_data),
                      simplify = F)
  
  # grab the importance variables MDA
  model_important_vars <- sapply(scfas, 
                                 function(x) 
                                   grab_importance(x, j, rf_models, class_important_vars), simplify = F)
  
  # Run the prediction
  rf_test_data <- sapply(scfas, 
                         function(x) 
                           run_prediction(x, rf_models, "test_data", "raw", "disease", final_rf_data), 
                         simplify = F)
  
  # Generate the ROC curves
  test_roc <- sapply(scfas, 
                     function(x) get_test_roc(x, "disease", "tempPredictions", rf_test_data), simplify = F)
  
  # Generate the summary data
  class_summary_data <- sapply(scfas, function(x) 
    add_model_summary_data(x, rf_models, test_roc, class_summary_data), 
    simplify = F)
  
}

# Write out the necessary data files
sapply(scfas, function(x) 
  write_csv(class_summary_data[[x]], 
            paste("data/process/tables/", x, "_classification_RF_summary.csv", sep = "")))

sapply(scfas, function(x) 
  write_csv(class_important_vars[[x]], 
            paste("data/process/tables/", x, "_imp_otus_classification_RF_summary.csv", sep = "")))

