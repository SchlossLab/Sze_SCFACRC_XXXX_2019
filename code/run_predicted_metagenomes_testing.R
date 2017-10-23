### Investigate how well the picrust predicted metagenomes work
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

# Load in predicted metagenome data
all_biom_data <- read_biom("data/process/predicted_metagenomes.biom")

# Converts the needed biom data into a matrix
full_pi_data <- t(as(biom_data(all_biom_data), "matrix")) %>% 
  as.data.frame() %>% 
  mutate(sample_id = rownames(.)) %>% 
  select(sample_id, everything())

# Get the kegg mappings for each coded sample
kegg_table <- observation_metadata(all_biom_data) %>% 
  mutate(kegg_id = rownames(.)) %>% 
  select(kegg_id, everything())

# Generate the near zero variance
nzv <- nearZeroVar(full_pi_data)

# Remove these from the data set
nzv_full_pi_data <- full_pi_data[, -nzv]

##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################

# Function to create data groups to be used in the RF model
make_rf_data <- function(dataTable, metadata, exclude_samples, sep_name){
  
  tempData <- metadata %>% 
    filter(!(Group %in% exclude_samples), dx != sep_name) %>% 
    select(Group, dx) %>% 
    inner_join(dataTable, by = c("Group" = "sample_id")) %>% 
    select(-Group) %>% 
    rename(disease = dx)
  
  return(tempData)
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
model_important_vars <- list(adenoma = data_frame(), cancer = data_frame())

model_summary_data <- list(adenoma = data_frame(), cancer = data_frame())

# Set up model type variables
models <- c("adenoma", "cancer")
# Create RF data tables to be used
rf_data <- sapply(models, 
                  function(x) make_rf_data(nzv_full_pi_data, picrust_meta, 
                                           c(metaF$initial, metaF$followUp), x), simplify = F)
# Correctly assign names
names(rf_data) <- c("cancer", "adenoma")

# run through a 100 different 80/20 splits for modeling
for(j in 1:100){
  
  # Create test and train data for both conditions
  final_rf_data <- sapply(models, 
                          function(x) eighty_twenty_split(x, rf_data), simplify = F)
  
  # Create the two models
  rf_models <- sapply(models, 
                      function(x) make_rf_model(x, j, "training_data", "disease", "rf", "ROC", final_rf_data),
                      simplify = F)
  
  # grab the importance variables MDA
  model_important_vars <- sapply(models, 
                                 function(x) 
                                   grab_importance(x, j, rf_models, model_important_vars), simplify = F)
  
  # Run the prediction
  rf_test_data <- sapply(models, 
                         function(x) 
                           run_prediction(x, rf_models, "test_data", "raw", "disease", final_rf_data), 
                         simplify = F)
  
  # Generate the ROC curves
  test_roc <- sapply(models, 
                     function(x) get_test_roc(x, "disease", "tempPredictions", rf_test_data), simplify = F)
  
  # Generate the summary data
  model_summary_data <- sapply(models, function(x) 
    add_model_summary_data(x, rf_models, test_roc, model_summary_data), 
    simplify = F)
  
}

# Write out the necessary data files
sapply(models, function(x) 
  write_csv(model_summary_data[[x]], 
            paste("data/process/tables/", x, "_classification_RF_summary.csv", sep = "")))

sapply(models, function(x) 
  write_csv(model_important_vars[[x]], 
            paste("data/process/tables/", x, "_imp_otus_classification_RF_summary.csv", sep = "")))

# Write out the relevant Kegg ID key
write_csv(kegg_table, "data/process/tables/kegg_id_key.csv")

