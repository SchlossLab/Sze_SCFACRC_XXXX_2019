# Load needed libraries
library("tidyverse")
library("caret")
library("pROC")
library("randomForest")
library("data.table")

# Function to run Begum's pipeline
pipeline <- function(dataset){

# Create vectors to save cv and test AUC values for every data-split
results_total <-  data.frame()
test_aucs <- c()
cv_aucs <- c()

# We are doing the pre-processing to the full dataset and then splitting 80-20
# Scale all features between 0-1
preProcValues <- preProcess(dataset, method = "range")
dataTransformed <- predict(preProcValues, dataset)

# remove columns that only appear within one or fewer samples of the training set. These are
# likely to be all zero and will not enter into the model
frequent <- names(which(apply(dataTransformed[, -1] > 0, 2, sum) > 1))

dataTransformed <- dataTransformed %>% select(regress, frequent)

# Do the 80-20 data-split
# Stratified data partitioning %80 training - %20 testing

inTraining <- createDataPartition(dataTransformed$regress, p = .80, list = FALSE)
trainTransformed <- dataTransformed[ inTraining,]
testTransformed  <- dataTransformed[-inTraining,]

	n_features <- ncol(trainTransformed) - 1
	if(n_features > 20000) n_features <- 20000

	if(n_features < 19){ mtry <- 1:6
	} else { mtry <- floor(seq(1, n_features/3, length=6)) }

	mtry <- mtry[mtry <= n_features]

# cv index to make sure the internal 5-folds are stratified for diagnosis classes and also resampled 100 times.
# 100 repeat internally is necessary to get robust readings of hyperparameter setting performance
  folds <- 5
	cvIndex <- createMultiFolds(factor(trainTransformed$regress), folds, times=100) #returnTrain = T default for multifolds
  cv <- trainControl(method="repeatedcv",
                     number= folds,
                     index = cvIndex,
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
	cv_best <- getTrainPerf(trained_model)

  # Training results for all the mtry parameters
  cv_results <- trained_model$results

  # Predictions with the best mtry parameter for 1 data-split
  predictions <- predict(trained_model, testTransformed, type = "raw")

  # quality values for test data from 1 datasplit
	test_results <- postResample(testTransformed$regress, predictions)

  results <- list(cv_best=cv_best, cv_all=cv_results, test=test_results)
  return(results)
}

# Function to save the RMSE values and save them as .csv
get_RMSE_R2_MAE <- function(dataset, split_number, dir){

  # Save results of the modeling pipeline as a list
  results <- pipeline(dataset)

  # ------------------------------------------------------------------
  # Create a matrix with cv_aucs and test_aucs from 100 data splits
	colnames(results$cv_best) <- str_replace(colnames(results$cv_best), "Train", "train_")
	names(results$test) <- paste0("test_", names(results$test))


  # Convert to dataframe and add a column noting the model name
  enframe(unlist(c(results$cv_best[1, -4], results$test))) %>%
		spread(name, value) %>%
    write_csv(path=paste0(dir, "/optimum_mtry.", split_number, ".csv"))

  # ------------------------------------------------------------------

  # ------------------------------------------------------------------
  # Save all tunes from 100 data splits and corresponding AUCs
  results$cv_all %>%
    write_csv(path=paste0(dir, "/all_mtry.", split_number, ".csv"))
  # ------------------------------------------------------------------

}

##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

read_shared <- function(shared_file_name, min_samples=2, min_abundance=1){

	data <- fread(shared_file_name, header=T, colClasses=c(Group="character"))[, -c(1,3)]

	if(min_abundance != 1){
		abundant <- names(which(apply(data[,-1], 2, sum) >= min_abundance))
		data <- select(data, Group, abundant)
	}

	frequent <- names(which(apply(data[, -1] > 0, 2, sum) >= min_samples))

	select(data, Group, frequent)

}

get_data <- function(path){

	parse_path <- unlist(str_split(unlist(str_split(path, '/'))[3], "_"))
	target_scfa <- parse_path[1]
	feature_sources <- parse_path[-1]

	metagenomics <- c("opf", "kegg")
	tax_levels <- c("kingdom", "phylum", "class", "order", "family", "genus")
	picrust2 <- c("pc2ko", "pc2ec", "pc2pathways")

	# Read in metadata
	if(any(feature_sources %in% metagenomics)){
		data <- read_tsv('data/metadata/zackular_metadata.tsv') %>%
						mutate(sample=str_replace(sample, "(\\d{1,2})", " \\1")) %>%
						separate(sample, into=c("disease", "subject")) %>%
						mutate(disease=str_replace(disease, "^(.).*", "\\1"),
										dx=tolower(dx)) %>%
						unite(sample, disease, subject, sep="") %>%
						select(sample, fit_result, dx)
	} else {
		data <- read_csv('data/metadata/cross_section.csv', col_types=cols(sample=col_character()))
	}

	if("fit" %in% feature_sources){
		data <- data %>% select(sample, fit_result)
	} else {
		data <- data %>% select(sample)
	}

	if("asv" %in% feature_sources){

		data <- read_shared('data/asv/crc.asv.shared') %>%
			inner_join(data, ., by=c("sample"="Group"))

	}

	# Read in OTU table and remove label and numOtus columns
	if("otu" %in% feature_sources){

		data <- read_shared('data/mothur/crc.otu.shared') %>%
			inner_join(data, ., by=c("sample"="Group"))

	}

	if(any(feature_sources %in% tax_levels)){
		taxon <- feature_sources[which(feature_sources %in% tax_levels)]
		shared_taxon_file <- paste0('data/phylotype/crc.', taxon, ".shared")

		data <- read_shared(shared_taxon_file) %>%
			inner_join(data, ., by=c("sample"="Group"))

	}

	if(any("picrust1" %in% feature_sources)){

		data <- read_shared("data/picrust1/crc.picrust1.shared") %>%
			inner_join(data, ., by=c("sample"="Group"))

	}

	if(any(feature_sources %in% picrust2)){
		pc_tag <- feature_sources[which(feature_sources %in% picrust2)]
		tag <- str_replace(pc_tag, "pc2", "")

		picrust_file_name <- paste0("data/picrust2/crc.", tag, ".shared")

		data <- read_shared(picrust_file_name) %>%
			inner_join(data, ., by=c("sample"="Group"))

	}

	if(any(feature_sources %in% metagenomics)){
		mg_tag <- feature_sources[which(feature_sources %in% metagenomics)]
		mg_file_name <- paste0("data/metagenome/metag.", mg_tag, ".shared")

		data <- read_shared(mg_file_name, min_abundance=5*78) %>%
			inner_join(data, ., by=c("sample"="Group"))

	}

	# Read in SCFAs spread columns
	read_tsv('data/scfa/scfa_composite.tsv', col_types=cols(study_id=col_character())) %>%
		spread(key=scfa, value=mmol_kg) %>%
		select(study_id, target_scfa) %>%
		inner_join(., data, by=c("study_id" = "sample")) %>%
		rename(regress = target_scfa) %>%
		drop_na() %>%
		select(-study_id)
}

######################## RUN PIPELINE #############################
# Get the cv and test AUCs for 100 data-splits

input <- commandArgs(trailingOnly=TRUE) # recieve input from model
# Get variables from command line
seed <- as.numeric(input[1])
path <- input[2]

if(!dir.exists(path)){
	dir.create(path, recursive=TRUE)
}

start_time <- Sys.time()

set.seed(seed)
data <- get_data(path)
get_RMSE_R2_MAE(data, seed, path) # model will be "Random_Forest"

end_time <- Sys.time()
print(end_time - start_time)
###################################################################
