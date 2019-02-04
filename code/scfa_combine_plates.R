### Combine the different scfa plates into a single table
### Will perform this for every scfa analyzed
### Marc Sze



# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "ggplot2"))


# Set up the scfas of interest
scfas <- c("acetate", "butyrate", "isobutyrate", "propionate")

#Set up plate name
plates <- paste("transformed_plate", c(1:8), "_scfa_crc", sep = "")

# Set up raw text file path
raw_path <- "data/raw/Raw_hplc_files/"

# Create write out path
write_out_path <- "data/process/tables/"




##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################

# Function to read in data and select out specific columns
read_in_data <- function(scfa_of_int, plates, raw_path){
  
  tempData <- sapply(plates, 
                     function(x) 
                       read.csv(paste(raw_path, scfa_of_int, "/", 
                                      x, "_", scfa_of_int, ".csv", sep = ""), 
                                header = T, stringsAsFactors = F) %>% 
                       select(study_id, mmol_kg) %>% 
                       mutate(study_id = as.character(study_id)), simplify = F)
  
  return(tempData)
  
}






##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################


# combines data by type of scfa
testData <- sapply(scfas, 
                   function(x) 
                     read_in_data(x, plates = plates, raw_path = raw_path) %>% 
                     bind_rows(), simplify = F)
  
# writes out data based on scfa
sapply(scfas, 
       function(x) 
         write.csv(testData[[x]], 
                   paste(write_out_path, x, "_final_data.csv", sep = ""), 
                   row.names = F))










