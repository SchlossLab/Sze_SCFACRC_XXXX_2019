### Combine the different scfa plates into a single table
### Will perform this for every scfa analyzed
### Marc Sze

# Load needed libraries
library("tidyverse")


# Set up the scfas of interest
scfas <- c("acetate", "butyrate", "isobutyrate", "propionate")

#Set up plate name
plates <- paste("transformed_plate", c(1:8), "_scfa_crc", sep = "")

# Set up file path
path <- "data/scfa/"



##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################

# Function to read in data and select out specific columns
read_in_data <- function(scfa_of_int, plates, path){

  tempData <- sapply(plates,
                     function(x)
                       read.csv(paste(path, "/",
                                      x, "_", scfa_of_int, ".csv", sep = ""),
                                header = T, stringsAsFactors = F) %>%
                       select(study_id, mmol_kg) %>%
                       mutate(study_id = as.character(study_id),
										 					scfa = scfa_of_int), simplify = F)

  return(tempData)

}


##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

# combines data by type of scfa
sapply(scfas,
         function(x)
           read_in_data(x, plates = plates, path = path) %>%
           bind_rows(), simplify = F) %>%
			bind_rows() %>%
			group_by(study_id, scfa) %>%
			summarize(mmol_kg = mean(mmol_kg)) %>%
			ungroup() %>%
			spread(scfa, mmol_kg, fill=0) %>%
			mutate(pooled = 4*butyrate + 4*isobutyrate + 3*propionate + 2*acetate) %>%
			gather(scfa, mmol_kg, -study_id) %>%
			write_tsv("data/scfa/scfa_composite.tsv")
