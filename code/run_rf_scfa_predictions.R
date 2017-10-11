### Investigate how well one can predict SCFAs
### Using RF models can we ID what predicts SCFAs
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "caret"))

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
eighty_twenty_split <- function(i, data_of_int, dataList){
  
  tempData <- dataList[[i]][[data_of_int]]
  
  totalLength <- length(rownames(tempData))
  
  tempTrain <- tempData %>% sample_frac(0.8, replace = FALSE)
  
  tempTest <- tempData %>% filter(!(rownames(.) %in% rownames(tempTrain)))
  
  finalList <- list(
    training_data = tempTrain, 
    test_data = tempTest)
  
  
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
  fitControl <- trainControl(## 10-fold CV
    method = "cv",
    number = 10,
    p = 0.8, 
    classProbs = TRUE, 
    summaryFunction = twoClassSummary, 
    savePredictions = "final")

  
  #Train the model
  training_model <- 
    train(formula(paste(cat_column, " ~ .", sep = "")), data = tempData, 
          method = method_to_use, 
          ntree = 500, 
          trControl = fitControl,
          metric = metric_to_use, 
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
  
  tempData <- dataList[[i]][[train_data_name]]
  tempModel <- model_data[[i]]
  
  
  
  tempPredictions <- predict(tempModel, tempData, type = control_type)
  
  tempPredictions <- cbind(tempPredictions, actual = tempData[, var_of_int])
  
  
  return(tempPredictions)
}


# Generate ROC curves for the test data
get_test_roc <- function(i, train_data_name, var_of_int, 
                         pred_list, dataList){
  
  pred_roc[[paste("run_", i, sep = "")]] <- 
    roc(test_data_list[[paste("run_", i, sep="")]]$lesion ~ 
          factor(run_predictions[[paste("run_", i, sep = "")]], ordered = TRUE))
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

#split the data for the downstream nzv that needs to occur
final_sp_data <- sapply(scfas, 
                        function(x) split_dataframe(x, final_data), simplify = F)
# Generate an 80/20 data split
rf_data <- sapply(scfas, function(x) eighty_twenty_split(x, "rf_groups", final_sp_data), simplify = F)

class_test <- sapply(scfas, function(x) 
                       make_rf_model(x, 1, "training_data", "high_low", "rf", "ROC", rf_data), 
                     simplify = F)

class_pred <- sapply(scfas, function(x) 
                       run_prediction(x, class_test, "test_data", "prob", "high_low", rf_data), 
                     simplify = F)


# Generate an 80/20 data split for regression
reg_rf_data <- sapply(scfas, 
                  function(x) eighty_twenty_split(x, "rf_regression", final_sp_data), simplify = F)

regression_test <- sapply(scfas, 
                          function(x) make_rf_model(x, 1, "training_data", 
                                 "mmol_kg", "rf", "Rsquared", reg_rf_data), simplify = F)

reg_pred <- sapply(scfas, function(x) 
                     run_prediction(x, regression_test, "test_data", "raw", "mmol_kg", reg_rf_data), 
                   simplify = F)





