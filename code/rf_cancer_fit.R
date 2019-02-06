######################################################################
# Author: Begum Topcuoglu; modified by Pat Schloss
# Date: 2018-12-20
# Title: Main pipeline in R programming language
######################################################################
#
# Usage in command-line:
#   Rscript code/main_RF.R $seed data/rf/otu_fit
#
######################################################################
# Description:

# This script will read in data from Baxter et al. 2016
#     - 0.03 subsampled OTU dataset
#     - CRC metadata: SRN information


# It will run the following machine learning pipelines:
#   Random Forest
######################################################################

######################################################################
# Dependencies and Outputs:

# Be in the project directory.

# The outputs are:
#   (1) AUC values for cross-validation and testing for each data-split
#   (2) meanAUC values for each hyper-parameter tested during each split.
######################################################################


################### IMPORT LIBRARIES and FUNCTIONS ###################
# The dependinces for this script are consolidated in the first part
deps = c("randomForest", "reshape2", "kernlab","LiblineaR", "doParallel","pROC", "caret", "gtools", "tidyverse");
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE, repos = "http://cran.us.r-project.org");
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Load in needed functions and libraries
source('code/rf_pipeline/classification_model_selection.R')
source('code/rf_pipeline/classification_pipeline.R')
source('code/rf_pipeline/classification_generate_AUCs.R')

######################################################################

######################## DATA PREPARATION #############################
# Features: Hemoglobin levels and 16S rRNA gene sequences in the stool
# Labels: - Colorectal lesions of 490 patients.
#         - Defined as lesion or not.(lesin here means: adenoma and cancer)
# Read in metadata
meta <- read_csv('data/raw/metadata/cross_section.csv', col_types=cols(sample=col_character())) %>%
  select(sample, dx, fit_result)

data <- meta %>%
  mutate(classes = case_when(
		dx == "normal" ~ "control",
		dx == "adenoma" ~ "NA",#lesion
    dx == "cancer" ~ "case",#lesion
		TRUE ~ "NA"
  )) %>%
	mutate(classes = factor(classes, levels=c("control", "case"))) %>%
	select(classes, fit_result) %>%
  drop_na()

###################################################################

######################## RUN PIPELINE #############################
start_time <- Sys.time()

input <- commandArgs(trailingOnly=TRUE) # recieve input from model
# Get variables from command line
seed <- as.numeric(input[1])
path <- input[2]
model <- "Random_Forest"

if(!dir.exists(path)){
	dir.create(path, recursive=TRUE)
}

set.seed(seed)
get_AUCs(data, model, seed, path)

end_time <- Sys.time()
print(end_time - start_time)
###################################################################
