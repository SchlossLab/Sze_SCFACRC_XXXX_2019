### Prelim analysis of measured SCFAs by HPLC
### cursory look at SCFAs and Tumor groups
### Marc Sze



# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "ggplot2", "gridExtra"))

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
    select(study_id, Dx_Bin, dx) %>% 
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


#Run an ANOVA with Tukey post hoc test to look for differences within each group
get_anova_comparisons <- function(scfa_name, dataList, 
                                  set_variable = "mmol_kg", set_groups = "Dx_Bin"){
  # scfa_name is a variable that stores the scfa of interest e.g. "acetate"
  # dataList is a list of combined data needed for the comparisons
  # set_variable is the column name of the variable of interest
  # set_groups specifies what column to use for the groups variable
  
  # Create a temporary data table with the scfa of interest
  tempData <- dataList[[scfa_name]]
  
  # this runs an ANOVA with Tukey post hoc test and makes a readable table 
  test_data <- summary(aov(lm(
    as.formula(paste(set_variable, "~", set_groups)), data = tempData)))[[1]]
  # Returns the test results to the global working environment
  return(test_data)
  
}



# Function to save ANOVA results
save_anova_comparisons <- function(scfa_name, dataList, 
                                   path_to_save, ending){
  # scfa_name is a variable that stores th scfa of interest e.g. "acetate"
  # dataList is a list of all the comparison tables by scfa
  # path_to_save is the directory path where the file will go
  # ending is the unique file name ending
  
  # pulls the respective comparison of interest and transforms the column name
  tempData <- dataList[[scfa_name]] %>% 
    rename(pvalue = `Pr(>F)`)
  # writes file to computer to directory that was specified
  write.csv(tempData, 
            paste(path_to_save, scfa_name, ending, sep = ""), 
            row.names = F)
  
}



##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

# Load in needed metadata and rename sample column
metaI <- read.csv("data/raw/metadata/metaI_final.csv", 
                  header = T, stringsAsFactors = F) %>% 
  rename(study_id = sample)
# read in all the respective scfa data tables
scfa_data <- sapply(scfas, 
                    function(x) upload_scfa_data(x, "data/process/tables/", "_final_data.csv"), 
                    simplify = F)
# combine all the necessary data together into a single file
combined_data <- sapply(scfas, 
                        function(x) merge_data(x, metaI, scfa_data, not_measured), simplify = F)

# run a krsukal wallis test on every scfa of interest
test <- sapply(scfas, 
               function(x) get_anova_comparisons(x, combined_data), simplify = F)

# write out all the kruskal tests to computer
lapply(scfas, 
       function(x) 
         save_anova_comparisons(x, test, 
                                "data/process/tables/", 
                                "_kruskal_crc_groups.csv"))









