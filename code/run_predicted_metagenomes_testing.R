### Investigate how well the picrust predicted metagenomes work
### compare most differentially expressed genes by disease
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "biomformat", "caret", "randomForest"))

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

# Set up initial store for summary data
model_important_vars <- list(adenoma = data_frame(), cancer = data_frame())

# Set up model type variables
models <- c("adenoma", "cancer")
# Create RF data tables to be used
rf_data <- sapply(models, 
                  function(x) make_rf_data(nzv_full_pi_data, picrust_meta, 
                                           c(metaF$initial, metaF$followUp), x), simplify = F)
# Correctly assign names
names(rf_data) <- c("cancer", "adenoma")


# Create test and train data for both conditions
final_rf_data <- sapply(models, 
                        function(x) eighty_twenty_split(x, rf_data), simplify = F)

# Create the two models
rf_models <- sapply(models, 
       function(x) make_rf_model(x, 1, "training_data", "disease", "rf", "ROC", final_rf_data),
       simplify = F)

