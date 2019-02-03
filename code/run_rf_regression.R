# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "caret", "pROC", "randomForest"))

#setup variables that will be used
scfas <- c("acetate", "butyrate", "isobutyrate", "propionate")

# Samples that had to little volume to actually be measured
not_measured <- c("3367653", "2041650", "2043650", "3027650")

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
  
  temp_regression <- tempData %>% select(mmol_kg, Group, dx, contains("Otu"))
  
  temp_high_low <- tempData %>% select(high_low, Group, dx, contains("Otu"))
  
  finalList <- list(rf_regression = temp_regression, 
                    rf_groups = temp_high_low)
  return(finalList)
}


# Function to make an initial 80/20 split within the data
eighty_twenty_split <- function(i, data_of_int, dataList){
  
  tempData <- dataList[[i]][[data_of_int]]
  
  sampling_vector <- c(1:length(rownames(tempData)))
  
  trainValues <- sample(sampling_vector, round(length(sampling_vector)*0.8))
  
  #tempTrain <- tempData %>% sample_frac(0.8, replace = FALSE)
  
  #tempTest <- tempData %>% filter(!(rownames(tempTrain) %in% rownames(.)))
  
  finalList <- list(
    training_data = tempData[trainValues, ], 
    test_data = tempData[-trainValues, ])
  
  
  return(finalList)
  
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

#split the data for each SCFA 
final_sp_data <- sapply(scfas, 
                        function(x) split_dataframe(x, final_data), simplify = F)

# We will only analyze acetate downstream
final_acetate_data <- final_sp_data$acetate$rf_regression %>% 
  select(-Group, -dx)

pipeline <- function(dataset, model){
  # Create vectors to save cv and test AUC values for every data-split
  results_total <-  data.frame()
  test_aucs <- c()
  cv_aucs <- c()
  # Stratified data partitioning %80 training - %20 testing
  inTraining <- createDataPartition(dataset$mmol_kg, p = .80, list = FALSE)
  training <- dataset[ inTraining,]
  testing  <- dataset[-inTraining,]
  # Scale all features between 0-1
  preProcValues <- preProcess(training, method = "range")
  trainTransformed <- predict(preProcValues, training)
  testTransformed <- predict(preProcValues, testing)

  # Cross-validation method
  cv <- trainControl(method="repeatedcv",
                   repeats = 10,
                   number=10,
                   returnResamp="final",
                   classProbs=FALSE,
                   indexFinal=NULL,
                   savePredictions = TRUE)
  # Hyper-parameter tuning budget
  grid <-  expand.grid(mtry = c(80,500, 1000, 1500))
  # Train model for 1 data-split but with 10fold-10repeat CV
  trained_model <-  train(mmol_kg ~ .,
                          data=trainTransformed,
                          method = "rf",
                          trControl = cv,
                          metric = "RMSE",
                          tuneGrid = grid,
                          ntree=1000)
  # Mean AUC value over repeats of the best mtry parameter during training
  cv_RMSE <- getTrainPerf(trained_model)$TrainRMSE
  # Training results for all the mtry parameters
  model_results <- trained_model$results
  # Predictions with the best mtry parameter for 1 data-split
  predictions <- predict(trained_model, testTransformed, type = "raw")
  realval_predictions <- cbind(predictions, actual = testTransformed$mmol_kg)
  # RMSE values for test data from 1 datasplit
  test_RMSE <- RMSE(predictions,testTransformed$mmol_kg)
  
  results <- list(cv_RMSE, test_RMSE, model_results)
  return(results)
}

######################################################################
#------------------------- DEFINE FUNCTION -------------------#
######################################################################
get_AUCs <- function(models, split_number){
  for(ml in models){
    
    # Save results of the modeling pipeline as a list
    results <- pipeline(final_acetate_data, ml) 
    
    # ------------------------------------------------------------------ 
    # Create a matrix with cv_aucs and test_aucs from 100 data splits
    RMSE_results <- matrix(c(results[[1]], results[[2]]), ncol=2) 
    # Convert to dataframe and add a column noting the model name
    RMSE_results_dataframe <- data.frame(RMSE_results) %>% 
      rename(cv_RMSE=X1, test_RMSE=X2) %>% 
      mutate(model=ml) %>% 
      write.csv(file=paste0("data/best_hp_results_", ml,"_", split_number, ".csv"), row.names=F)
    # ------------------------------------------------------------------   
    
    # ------------------------------------------------------------------   
    # Save all tunes from 100 data splits and corresponding AUCs
    all_results <- results[3]
    # Convert to dataframe and add a column noting the model name
    dataframe <- data.frame(all_results) %>% 
      mutate(model=ml) %>% 
      write.csv(file=paste0("data/all_hp_results_", ml,"_", split_number, ".csv"), row.names=F)
    # ------------------------------------------------------------------ 
  }
}

######################## RUN PIPELINE #############################
# Get the cv and test AUCs for 100 data-splits
start_time <- Sys.time()

input <- commandArgs(trailingOnly=TRUE) # recieve input from model
# Get variables from command line
seed <- as.numeric(input[1])
model <- input[2]

set.seed(seed)
get_AUCs(model, input[1]) # model will be "Random_Forest"

end_time <- Sys.time()
print(end_time - start_time)
###################################################################

