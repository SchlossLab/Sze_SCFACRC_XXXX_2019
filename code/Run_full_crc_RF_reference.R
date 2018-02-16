### Build the best adenoma model possible
### Running the model
## Marc Sze


#Load needed libraries
source('code/functions.R')

loadLibs(c("tidyverse", "caret","scales"))

# Set i variable from command line
input_cmds <- commandArgs(TRUE)
i <- as.numeric(input_cmds[1])


# Read in the needed data
test_data <- as.data.frame(read_csv("data/process/tables/crc_full_test_data.csv"))
eighty_twenty_splits <- as.data.frame(read_csv("data/process/tables/crc_full_eighty_twenty_splits.csv"))

#################################################################################
#                                                                               #
#                                                                               #
#               Model Training and Parameter Tuning                             #
#                                                                               #
#################################################################################

#Get train data
train_test_data <- test_data[eighty_twenty_splits[, i], ]

# Get test data
test_test_data <- test_data[-eighty_twenty_splits[, i], ]

#Create Overall specifications for model tuning
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated twenty times
  repeats = 20, 
  p = 0.8, 
  classProbs = TRUE, 
  summaryFunction = twoClassSummary)


#################################################################################
#                                                                               #
#                                                                               #
#               Model Training and Parameter Tuning                             #
#                                                                               #
#################################################################################

#Set up lists to store the data
#test_tune_list <- list()
#test_predictions <- list()

#Train the model
#train_name <- paste("data_split", i, sep = "")

#set.seed(3457)
#test_tune_list[[paste("data_split", i, sep = "")]] <- assign(train_name, 
#                                                             train(dx ~ ., data = train_test_data, 
#                                                                   method = "rf", 
#                                                                   ntree = 500, 
#                                                                   trControl = fitControl, 
#                                                                   metric = "ROC", 
#                                                                   verbose = FALSE))

#test_predictions[[paste("data_split", i, sep = "")]] <- 
#  predict(test_tune_list[[paste("data_split", i, sep = "")]], 
#          test_test_data)


# Save image with data and relevant parameters
#save.image(paste("exploratory/adn_RF_model_", i, ".RData", sep=""))