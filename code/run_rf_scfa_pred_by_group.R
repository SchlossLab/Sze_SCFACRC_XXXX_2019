### Investigate how important OTUs based on MDA change for each group
### Using RF models do OTUs change if in non-tumor, adenoma, or carcinoma samples
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "caret", "pROC", "randomForest"))

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

log_shared <- read_tsv("data/process/final.0.03.subsample.shared") %>% 
  mutate(Group = as.character(Group)) %>% 
  select(-label, -numOtus) %>% 
  mutate_at(vars(contains("Otu")), funs(log10(. + 1)))

# Generate IDs that are initial and follow up
initial_IDs <- as.character(metaF$initial)
follow_IDs <- as.character(metaF$followUp)

# read in all the respective scfa data tables
scfa_data <- sapply(scfas, 
                    function(x) upload_scfa_data(x, "data/process/tables/", "_final_data.csv"), 
                    simplify = F)

# Combine all data into a single list of lists and remove follow samples
all_data <- pare_down_data(c(initial_IDs, follow_IDs), scfa_data, shared, picrust_meta)