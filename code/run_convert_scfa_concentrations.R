### Creation of function to make the conversions
### Need to adjust readings to get the final SCFA concentrations
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "ggplot2"))

# Load in data file
dataTable <- read.csv("exploratory/scratch/reference_for_calcs/test_data.csv", 
                      header = T, stringsAsFactors = F)

# Set up the locations of the standards
stds_location <- list(std1 = c(2:9), std2 = c(43:49), 
                      std3 = c(79:85))

# Set up the locastions of the samples
sample_location <- list(g1 = c(11:41), g2 = c(50:77), 
                        g3 = c(86:109))


##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################

# Function to generate the areas needed for the respective correlation analysis
get_areas <- function(i, all_stds_data, locationList){
  # i is for the standard of interest
  # all_stds_data is a data table with only standards
  # locationList is the location of the different groups of standards
  
  # Grabs the areas of specific group of data based on location
  temp_areas <- all_stds_data %>% slice(match(locationList[[i]], location)) %>% 
    select(Area) %>% mutate(Area = as.numeric(Area))
  # print this new table out
  return(temp_areas)
}


# Function to get the standard values used
get_standards <- function(i, all_stds_data, locationList){
  # i is the standard of interest
  # all_stds_data is a data table with only the standards
  # locationList is a list of which standards group to grab
  
  # This grabs the specific standard concentrations measured
  temp_stands <- all_stds_data %>% 
    slice(match(locationList[[i]], location)) %>% 
    select(conc)
  # Prints this out to working global environment
  return(temp_stands)
}

# Function to get specific set of samples
get_samples <- function(g, locationList, all_data){
  # g is the group name of interest
  # locationList is where the respective members of the group are 
  # all_data is the data table with the raw data
  
  # This creates a new data table with location, sample name, and area variables
  tempData <- all_data %>% 
    separate(Data.Filename, c("temp1", "temp2"), sep = "_") %>% 
    rename(location = temp2) %>% 
    separate(location, c("location", "extra"), sep = "\\.") %>% 
    mutate(location = as.numeric(location), 
           Area = as.numeric(Area)) %>% 
    slice(match(locationList[[g]], location)) %>% 
    select(location, Sample.Name, Area)
  # prints the data table to the global working environment
  return(tempData)
  
}



# Function to generate the needed coefficients for std transformation
get_corr_coefficients <- function(i, areas, stds_amount){
  # i is the std of interest
  # areas is the list of data with the standards area measurement
  # stds_amount is the actual concentrations of each respective standard
  
  # this generates the needed coefficients for the transformation of area to
  # uncorrected concentrations
  tempCoeff <- lm(areas[[i]]$Area ~ stds_amount[[i]]$conc)$coefficients
  # prints the data to the global working environment
  return(tempCoeff)
  
}


# Function to convert the standards based on provided area
run_conversion <- function(group_length, std_areas, 
                           stds_coeff = coeff_list, std_list_length = NULL, stds = T){
  # group_length represents the group number of interest
  # std_areas represents the dataList with the areas of interest
  # std_coeff is defaulted to coeff_list since it is always the same
  # std_list_length needs a number of samples rather than standards are being converted
  # stds represents whether the samples are standards or not
  
  # checks to see if conversion is for standards
  if(stds == T){
    # Convert concentrations by standard groups
    tempConc <- (std_areas[[group_length]] - 
                   filter(stds_coeff, variables == "Intercept")[, group_length]) / 
      filter(stds_coeff, variables == "m")[, group_length]
    
  } else {
    # control variable
    x = 1
    tempList <-  NULL # create empty storage variable
    # continuously iterate through every standard group for each sample 
    while (x - std_list_length <= 0){
      # saves conversions for groups of samples by standard group to a list
      tempList[[x]] <- (std_areas[[group_length]]$Area - 
                          filter(stds_coeff, variables == "Intercept")[, x]) / 
        filter(stds_coeff, variables == "m")[, x]
      # increase the control variable by 1
      x = x + 1

    }
    # convert the list to a data frame for the non-standard samples
    tempConc <- as.data.frame.list(tempList)
    colnames(tempConc) <- paste("std", c(1:std_list_length), sep = "") # set column names
    # add location and sample names to the concentrations
    tempConc <- cbind(std_areas[[group_length]], tempConc)
    
  }
  # print the tempConc to the global working environment
  return(tempConc)
  
}


# Function to create correction factor for the samples
get_correction_factor <- function(g, samplesList){
  # g is the group of interest
  # samplesList is the respective samples by grouping
  
  # set up vector to store the amount the corrective factor will be multiplied by
  tempVector <- c(1:length(rownames(samplesList[[g]])))
  # get the data table of interest based on group
  tempData <- samplesList[[g]]
  # generate the first correction factor
  first_cf <- 1/(max(tempData$location) - min(tempData$location))
  # generate the resulting correction factor for every sample
  tempVector <- tempVector * first_cf
  # print the correction factor to the global working environment
  return(tempVector)
}

# corrected_concentration = std1_value + correction_factor*(std2-std1)

# Function to get the correction concentrations
get_corr_conc <- function(list_number, conv_sampleData, 
                          correctionFactors){
  # list_number is a variable that represents the numerical group of interest
  # conv_sampleData is a list of data tables with the converted values
  # correctionFactors is a list with the unique correction facto separated by group location
  
  # sets up the initial data needed for the conversion
  value1 <- paste("std", list_number, sep = "") # first measures
  value2 <- paste("std", list_number + 1, sep = "") # second measures
  tempData <- conv_sampleData[[list_number]] # data table of interest
  # checks to see if the standards are at the end of the run or not
  if(length(conv_sampleData) != as.numeric(gsub("std", "", value1))){
    # Generate the corrected concentrations
    tempConc <- tempData[, value1] + 
      correctionFactors[[list_number]]*(tempData[, value2] - tempData[, value1])
    # merge the corrected concentrations with location and sample name
    finalTable <- tempData %>% select(location, Sample.Name) %>% 
      mutate(conc_mM = tempConc)
    
  } else {
    # sets table to NULL if standards are at the end
    finalTable <- NULL
  }
  # prints the data to the global working environment
  return(finalTable)
  
}



##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

# Make a standard only table
std_areas <- dataTable %>% 
  filter(grepl("STD", Sample.ID) == TRUE & grepl("biorad", Sample.ID) != TRUE) %>% 
  separate(Data.Filename, c("temp1", "temp2"), sep = "_") %>% 
  rename(location = temp2) %>% 
  separate(Sample.Name, c("conc", "unit"), sep = " ") %>% 
  mutate(location = as.numeric(gsub("\\.lcd", "", location)), 
         conc = as.numeric(conc))

# Generate the standard areas
area_list <- sapply(names(stds_location), 
               function(x) get_areas(x, std_areas, stds_location), simplify = F)
# Generate the standard concentrations
stand_list <- sapply(names(stds_location), 
                     function(x) get_standards(x, std_areas, stds_location), 
                     simplify = F)
# Generate the coefficients for each standard grouping
coeff_list <- sapply(names(area_list), 
                     function(x) get_corr_coefficients(x, area_list, stand_list)) %>% 
  as.data.frame() %>% mutate(variables = c("Intercept", "m"))

# generate the converted standards values based on generated coefficients
converted_stds <- sapply(c(1:length(area_list)), 
                         function(x) run_conversion(x, area_list))

# generate sample groupings based on location on run
sample_list <- sapply(names(sample_location), 
                      function(x) 
                        suppressWarnings(
                          get_samples(x, sample_location, dataTable)), simplify = F)
# generate correction factors based on location on run
correction_factor_list <- sapply(names(sample_location), 
                                 function(x) get_correction_factor(x, sample_list))
# get the converted sample values based on location of standards run and location on run
converted_samples <- sapply(c(1:length(sample_location)), 
                           function(x) 
                             run_conversion(x, sample_list, 
                                            std_list_length = length(sample_location), 
                                            stds = F), simplify = F)
# get the final corrected concentrations
corrected_conc <- sapply(c(1:length(sample_location)), 
                         function(x) 
                           get_corr_conc(x, converted_samples, correction_factor_list), 
                         simplify = F) %>% bind_rows()









# TO DO LIST
### Convert to mmol/kg based on weight measured and total volume originally suspended in







