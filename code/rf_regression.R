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
		filter(v != 0) %>%
		pull(feature)

	training <- training %>% select(regress, zero_variance)
	testing <- testing %>% select(regress, zero_variance)

	n_features <- ncol(training) - 1


	if(n_features > 20000) n_features <- 20000

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
		data <- fread('data/asv/crc.asv.shared', header=T, colClasses=c(Group="character")) %>%
						as_tibble()
		 select(-label, -numOtus) %>%
		 inner_join(data, ., by=c("sample"="Group"))
	}

	# Read in OTU table and remove label and numOtus columns
	if("otu" %in% feature_sources){
		data <- read_tsv('data/mothur/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.shared', col_types=cols(Group=col_character())) %>%
		  select(-label, -numOtus) %>%
			inner_join(data,., by=c("sample"="Group"))
	}

	if(any(feature_sources %in% tax_levels)){
		taxon <- feature_sources[which(feature_sources %in% tax_levels)]
		shared_taxon <- paste0('data/phylotype/crc.', taxon, ".shared")

		data <- read_tsv(shared_taxon, col_types=cols(Group=col_character())) %>%
			select(-label, -numOtus) %>%
			inner_join(data, ., by=c("sample"="Group"))
	}

	if(any("picrust1" %in% feature_sources)){

		picrust_file_name <- "data/picrust1/crc.picrust1.shared"

		data <- fread(picrust_file_name, header=T) %>%
			as_tibble() %>%
			mutate(Group=as.character(Group)) %>%
			select(-label, -numOtus) %>%
			inner_join(data, ., by=c("sample"="Group"))
	}

	if(any(feature_sources %in% picrust2)){
		pc_tag <- feature_sources[which(feature_sources %in% picrust2)]
		tag <- str_replace(pc_tag, "pc2", "")

		picrust_file_name <- paste0("data/picrust2/crc.", tag, ".shared")

		data <- fread(picrust_file_name, header=T) %>%
			as_tibble() %>%
			mutate(Group=as.character(Group)) %>%
			select(-label, -numOtus) %>%
			inner_join(data, ., by=c("sample"="Group"))
	}

	if(any(feature_sources %in% metagenomics)){
		mg_tag <- feature_sources[which(feature_sources %in% metagenomics)]
		mg_file_name <- paste0("data/metagenome/metag.", mg_tag, ".shared")

		data <- fread(mg_file_name, header=T, colClasses=c(Group="character")) %>%
			as_tibble() %>%
			select(-label, -numOtus) %>%
			inner_join(data, ., by=c("sample"="Group"))
	}

	# Read in SCFAs spread columns
	read_tsv('data/scfa/scfa_composite.tsv', col_types=cols(study_id=col_character())) %>%
		spread(key=scfa, value=mmol_kg) %>%
		mutate(pooled = 4*butyrate + 4*isobutyrate + 3*propionate + 2*acetate) %>%
		select(study_id, target_scfa) %>%
		inner_join(., data, by=c("study_id" = "sample")) %>%
		rename(regress = target_scfa) %>%
		drop_na() %>%
		select(-study_id)
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
