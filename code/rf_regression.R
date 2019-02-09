# Load needed libraries
library("tidyverse")
library("caret")
library("pROC")
library("randomForest")

# Function to run Begum's pipeline
pipeline <- function(dataset){

	# Create vectors to save cv and test AUC values for every data-split
  results_total <-  data.frame()
  test_aucs <- c()
  cv_aucs <- c()

	# Stratified data partitioning %80 training - %20 testing
  inTraining <- createDataPartition(dataset$regress, p = 0.80, list = FALSE)
  training <- dataset[ inTraining,]
  testing  <- dataset[-inTraining,]

	# remove columns that have no variance within the training set. These are likely to be all zero
	# and will not enter into the model
	zero_variance <- training %>%
		gather(-regress, key="feature", value="value") %>%
		group_by(feature) %>%
		summarize(v=var(value)) %>%
		filter(v == 0) %>%
		pull(feature)

	training <- training %>% select(-zero_variance)
	testing <- testing %>% select(-zero_variance)

	n_features <- ncol(training) - 1

	# powers <- 1

	if(n_features < 19){ mtry <- 1:6
	} else { mtry <- floor(seq(1, n_features/3, length=6)) }

	mtry <- mtry[mtry <= n_features]

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
  grid <-  expand.grid(mtry = mtry)

  # Train model for 1 data-split but with 10fold-10repeat CV
  trained_model <-  train(regress ~ .,
                          data=trainTransformed,
                          method = "rf",
                          trControl = cv,
                          metric = "RMSE",
                          tuneGrid = grid,
                          ntree=1000)

  # RMSE value over repeats of the best mtry parameter during training
  cv_RMSE <- getTrainPerf(trained_model)$TrainRMSE

  # Training results for all the mtry parameters
  model_results <- trained_model$results

  # Predictions with the best mtry parameter for 1 data-split
  predictions <- predict(trained_model, testTransformed, type = "raw")
  realval_predictions <- cbind(predictions, actual = testTransformed$regress)

  # RMSE values for test data from 1 datasplit
  test_RMSE <- RMSE(predictions,testTransformed$regress)

  results <- list(cv_RMSE, test_RMSE, model_results)
  return(results)
}

# Function to save the RMSE values and save them as .csv
get_RMSE <- function(dataset, split_number, dir){

  # Save results of the modeling pipeline as a list
  results <- pipeline(dataset)

  # ------------------------------------------------------------------
  # Create a matrix with cv_aucs and test_aucs from 100 data splits
  RMSE_results <- matrix(c(results[[1]], results[[2]]), ncol=2)

  # Convert to dataframe and add a column noting the model name
  RMSE_results_dataframe <- data.frame(RMSE_results) %>%
    rename(cv_RMSE=X1, test_RMSE=X2) %>%
    write_csv(path=paste0(dir, "/optimum_mtry.", split_number, ".csv"))
  # ------------------------------------------------------------------

  # ------------------------------------------------------------------
  # Save all tunes from 100 data splits and corresponding AUCs
  all_results <- results[3]

	# Convert to dataframe and add a column noting the model name
  dataframe <- data.frame(all_results) %>%
    write_csv(path=paste0(dir, "/all_mtry.", split_number, ".csv"))
  # ------------------------------------------------------------------

}

##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

get_data <- function(path){

	parse_path <- unlist(str_split(unlist(str_split(path, '/'))[3], "_"))
	target_scfa <- parse_path[1]
	feature_sources <- parse_path[-1]


	# Read in metadata
	meta <- read_csv('data/raw/metadata/cross_section.csv',
									col_types=cols(sample=col_character())) %>%
	  select(sample, dx, fit_result)

	# Read in SCFAs spread columns
	scfa <- read_tsv('data/scfa/scfa_composite.tsv', col_types=cols(study_id=col_character())) %>%
		spread(key=scfa, value=mmol_kg)

	data <- inner_join(meta, scfa, by=c("sample"="study_id"))

	if("otu" %in% feature_sources){
	# Read in OTU table and remove label and numOtus columns
		data <- read_tsv('data/mothur/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.shared', col_types=cols(Group=col_character())) %>%
		  select(-label, -numOtus) %>%
			inner_join(data,., by=c("sample"="Group"))
	}

	if(setequal(feature_sources, c("otu", "fit"))){
		data <- data %>% select(target_scfa, fit_result, starts_with("Otu"))
	} else if(setequal(feature_sources, "fit")){
		data <- data %>% select(target_scfa, fit_result)
	} else if(setequal(feature_sources, "otu")){
		data <- data %>% select(target_scfa, starts_with("Otu"))
	}

	data %>%
		rename(regress = target_scfa) %>%
	  drop_na()
}

######################## RUN PIPELINE #############################
# Get the cv and test AUCs for 100 data-splits
start_time <- Sys.time()

input <- commandArgs(trailingOnly=TRUE) # recieve input from model
# Get variables from command line
seed <- as.numeric(input[1])
path <- input[2]

if(!dir.exists(path)){
	dir.create(path, recursive=TRUE)
}

set.seed(seed)
data <- get_data(path)
get_RMSE(data, seed, path) # model will be "Random_Forest"

end_time <- Sys.time()
print(end_time - start_time)
###################################################################
