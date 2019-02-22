######################################################################
# Author: Begum Topcuoglu; modified by Pat Schloss
# Date: 2018-12-20
# Title: Main pipeline in R programming language
######################################################################
#
# Usage in command-line:
#   Rscript code/main_RF.R $seed data/rf/cancer_otu_fit
#
######################################################################


################### IMPORT LIBRARIES and FUNCTIONS ###################
# The dependinces for this script are consolidated in the first part
deps = c("randomForest", "pROC", "caret", "tidyverse", "data.table");
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE, repos = "http://cran.us.r-project.org");
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

######################################################################

pipeline <- function(dataset){
  # Create vectors to save cv and test AUC values for every data-split
  results_total <-  data.frame()
  test_aucs <- c()
  cv_aucs <- c()
  # Do the 80-20 data-split

  # Stratified data partitioning 80% training - 20% testing
  inTraining <- createDataPartition(dataset$classes, p = 0.80, list = FALSE)


  training <- dataset[ inTraining,]
  testing  <- dataset[-inTraining,]

	# remove columns that have no variance within the training set. These are likely to be all zero
	# and will not enter into the model
	zero_variance <- training %>%
		gather(-classes, key="feature", value="value") %>%
		group_by(feature) %>%
		summarize(v=var(value)) %>%
		filter(v == 0) %>%
		pull(feature)

	training <- training %>% select(-zero_variance)
	testing <- testing %>% select(-zero_variance)


  # Scale all features between 0-1
  preProcValues <- preProcess(training, method = "range")
  trainTransformed <- predict(preProcValues, training)
  testTransformed <- predict(preProcValues, testing)

	# Define hyper-parameter tuning grid and the training method
	n_features <- ncol(training) - 1

	cv <- trainControl(method="repeatedcv",
                     repeats = 100,
                     number=5,
                     returnResamp="final",
                     classProbs=TRUE,
                     summaryFunction=twoClassSummary,
                     indexFinal=NULL,
                     savePredictions = TRUE)


	# seems like max mtry value should be less than the nubmer of samples
	if(n_features > 20000) n_features <- 20000

	if(n_features < 19){ mtry <- 1:6
	} else { mtry <- floor(seq(1, n_features/3, length=6)) }

  grid <-  expand.grid(mtry = mtry[mtry <= n_features])

	# Train the model
  trained_model <-  train(classes ~ .,
                          data=trainTransformed,
                          method = "rf",
                          trControl = cv,
                          metric = "ROC",
                          tuneGrid = grid,
                          ntree=1000)

	# Mean AUC value over repeats of the best cost parameter during training
	cv_auc <- getTrainPerf(trained_model)$TrainROC

	# Predict on the test set and get predicted probabilities
  rpartProbs <- predict(trained_model, testTransformed, type="prob")
  test_roc <- roc(ifelse(testTransformed$classes == "case", 1, 0),
                  rpartProbs[[2]])
  test_auc <- test_roc$auc

	# Save all the test AUCs over iterations in test_aucs
  test_aucs <- c(test_aucs, test_auc)

	# Cross-validation mean AUC value
  # Save all the cv meanAUCs over iterations in cv_aucs
  cv_aucs <- c(cv_aucs, cv_auc)

	# Save all results of hyper-parameters and their corresponding meanAUCs for each iteration
  results_individual <- trained_model$results
  results_total <- rbind(results_total, results_individual)

  results <- list(cv_aucs, test_aucs, results_total)
  return(results)
}

######################################################################

get_AUCs <- function(dataset, split_number, path){

  # Save results of the modeling pipeline as a list
  results <- pipeline(dataset)

  # ------------------------------------------------------------------
  # Create a matrix with cv_aucs and test_aucs from 100 data splits
  aucs <- matrix(c(results[[1]], results[[2]]), ncol=2)
  # Convert to dataframe and add a column noting the model name
  aucs_dataframe <- data.frame(aucs) %>%
    rename(cv_aucs=X1, test_aucs=X2) %>%
		write_csv(path=paste0(path, "/optimum_mtry.", split_number, ".csv"))

  # ------------------------------------------------------------------
  # Save all tunes from 100 data splits and corresponding AUCs
  all_results <- results[3]
  # Convert to dataframe and add a column noting the model name
  dataframe <- data.frame(all_results) %>%
		write_csv(path=paste0(path, "/all_mtry.", split_number, ".csv"))

  # ------------------------------------------------------------------
}

######################## DATA PREPARATION #############################
# Features: Hemoglobin levels and 16S rRNA gene sequences in the stool
# Labels: - Colorectal lesions of 490 patients.
#         - Defined as lesion or not.(lesin here means: adenoma and cancer)
#######################################################################

get_data <- function(path) {

	parse_path <- unlist(str_split(unlist(str_split(path, '/'))[3], "_"))
	classify <- parse_path[1]
	feature_sources <- parse_path[-1]

	metagenomics <- c("opf", "kegg")
	tax_levels <- c("kingdom", "phylum", "class", "order", "family", "genus")
	picrust2 <- c("pc2ko", "pc2ec", "pc2pathways")

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
		data <- data %>% select(sample, dx, fit_result)
	} else {
		data <- data %>% select(sample, dx)
	}

	# Read in SCFAs spread columns
	if("scfa" %in% feature_sources){
		data <- read_tsv('data/scfa/scfa_composite.tsv',
										col_types=cols(study_id=col_character())) %>%
		spread(key=scfa, value=mmol_kg) %>%
		inner_join(data, ., by=c("sample"="study_id"))
	}

	if("asv" %in% feature_sources){
		data <- fread('data/asv/crc.asv.shared', header=T, colClasses=c(Group="character")) %>%
						as_tibble()
		 select(-label, -numOtus) %>%
		 inner_join(data, ., by=c("sample"="Group"))
	}

	# Read in OTU table and remove label and numOtus columns
	if("otu" %in% feature_sources){
		data <- read_tsv('data/mothur/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.shared',
									col_types=cols(Group=col_character())) %>%
		 select(-label, -numOtus) %>%
		 inner_join(data, ., by=c("sample"="Group"))
	}

	# Read in phylotype data
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

	if(classify %in% c("adenoma", "cancer")){
		data <- data %>%
			mutate(classes = case_when(
								dx == "normal" ~ "control",
								dx == classify ~ "case",#lesion
								TRUE ~ NA_character_
							)
						)
	} else {
		data <- data %>%
			mutate(classes = case_when(
								dx == "normal" ~ "control",
								dx == "adenoma" ~ "case",#lesion
								dx == "cancer" ~ "case",#lesion
								TRUE ~ NA_character_
							)
						)
	}

	# Then remove the sample ID and dx columns
	data %>%
		select(-sample, -dx) %>%
		select(classes, everything()) %>%
	  drop_na()
}

######################## RUN PIPELINE #############################

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
get_AUCs(data, seed, path)

end_time <- Sys.time()
print(end_time - start_time)

###################################################################
