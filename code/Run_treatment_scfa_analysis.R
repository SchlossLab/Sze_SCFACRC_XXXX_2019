### Test with paired wilcoxson if scfas change after treatment
### Before and after treatment
### Marc Sze



# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "gridExtra"))

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
  tempData <- read.csv(paste(path_to_file, scfa_name, ending, sep = ""), 
                       header = T, stringsAsFactors = F)
  # sends the data to the global work environment
  return(tempData)
}


# Function to merge the respective scfa data with the needed metadata for visualization
merge_data <- function(scfa_name, meta_data, dataList, not_enough){
  # scfa_name is a variable that stores the scfa of interest e.g. "acetate"
  # meta_data contains the needed meta data on disease groups, etc.
  # dataList is a list with all the scfa data tables of interests
  # not_enough is a vector with the samples that could not be run on the HPLC
  
  # Create a temp data table variable
  tempData <- dataList[[scfa_name]]
  
  # Modify and merge the data together for samples that matched
  temp_combined_data <- meta_data %>% 
    select(study_id, Dx_Bin, dx, EDRN, time_point) %>% 
    mutate(study_id = as.character(study_id)) %>% 
    inner_join(tempData, by = "study_id") %>% 
    distinct(study_id, .keep_all = TRUE)
  
  # get sample IDs that were measured but had no signal
  # and remove those that were never run on HPLC
  test <- meta_data %>% 
    select(study_id, Dx_Bin, dx) %>% 
    mutate(study_id = as.character(study_id)) %>% 
    filter(!(study_id %in% temp_combined_data$study_id) & 
             !(study_id %in% not_enough)) %>% 
    mutate(mmol_kg = rep(0, length(study_id)))
  # create a final combined table with appropriate label names
  final_temp_combined <- temp_combined_data %>% 
    bind_rows(test) %>% 
    mutate(Dx_Bin = ifelse(Dx_Bin == "High Risk Normal", invisible("Normal"), 
                           ifelse(Dx_Bin == "adv Adenoma", 
                                  invisible("adv_Adenoma"), invisible(Dx_Bin))))
  # Send the final data table to the global work environment
  return(final_temp_combined)
  
}


# Function to run the paired wilcoxson test
run_wilcox_test <- function(int_scfa, dataList, disease){
  
  tempData <- as.data.frame(filter(dataList[[int_scfa]], dx == disease, !is.na(time_point)) %>% 
                              group_by(EDRN) %>% filter(n() == 2) %>% 
                              mutate(mmol_kg = ifelse(mmol_kg < 0, invisible(0), invisible(mmol_kg))))
  
  tempTest <- wilcox.test(filter(tempData, time_point == "initial")[, "mmol_kg"], 
                          filter(tempData, time_point == "followUp")[, "mmol_kg"], 
                          paired = TRUE)$p.value
  
  return(tempTest)
  
}



##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################


# Load in needed metadata and rename sample column
metaF <- read.csv("data/raw/metadata/good_metaf_final.csv", 
                  header = T, stringsAsFactors = F) %>% 
  gather("time_point", "study_id", initial, followUp) %>% 
  select(study_id, time_point, everything())

# read in all the respective scfa data tables
scfa_data <- sapply(scfas, 
                    function(x) upload_scfa_data(x, "data/process/tables/", "_final_data.csv"), 
                    simplify = F)

# combine all the necessary data together into a single file
combined_data <- sapply(scfas, 
                        function(x) merge_data(x, metaF, scfa_data, not_measured), simplify = F)

# Run the paired wilcox test
adn_test <- t(sapply(scfas, 
                     function(x) run_wilcox_test(x, combined_data, "adenoma"), simplify = F) %>% 
                bind_rows()) %>% as.data.frame() %>% 
  mutate(scfa = rownames(.), 
         bh = p.adjust(V1, method = "BH"), 
         disease = "adenoma") %>% 
  rename(pvalue = V1) %>% 
  select(scfa, disease, pvalue, bh)


crc_test <- t(sapply(scfas, 
                     function(x) run_wilcox_test(x, combined_data, "cancer"), simplify = F) %>% 
                bind_rows()) %>% as.data.frame() %>% 
  mutate(scfa = rownames(.), 
         bh = p.adjust(V1, method = "BH"), 
         disease = "carcinoma") %>% 
  rename(pvalue = V1) %>% 
  select(scfa, disease, pvalue, bh)

combined_table <- adn_test %>% bind_rows(crc_test)


write_csv(combined_table, "data/process/tables/scfa_treatment_comparison_pvalues.csv")





