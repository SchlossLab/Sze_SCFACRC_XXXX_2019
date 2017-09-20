### Creation of function to make the conversions
### Need to adjust readings to get the final SCFA concentrations
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "ggplot2"))


dataTable <- read.csv("exploratory/scratch/reference_for_calcs/test_data.csv", 
                      header = T, stringsAsFactors = F)


stds_location <- list(std1 = c(2:9), std2 = c(43:49), 
                      std3 = c(79:85))

sample_location <- list(g1 = c(11:41), g2 = c(50:77), 
                        g3 = c(86:109))


##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################

# Function to generate the areas needed for the respective correlation analysis
get_areas <- function(i, all_stds_data, locationList){
  
  temp_areas <- all_stds_data %>% slice(match(locationList[[i]], location)) %>% 
    select(Area) %>% mutate(Area = as.numeric(Area))
  
  return(temp_areas)
}


# Function to get the standard values used
get_standards <- function(i, all_stds_data, locationList){
  
  temp_stands <- all_stds_data %>% 
    slice(match(locationList[[i]], location)) %>% 
    select(conc)
  
  return(temp_stands)
}

# Function to get specific set of samples
get_samples <- function(g, locationList, all_data){
  
  tempData <- all_data %>% 
    separate(Data.Filename, c("temp1", "temp2"), sep = "_") %>% 
    rename(location = temp2) %>% 
    separate(location, c("location", "extra"), sep = "\\.") %>% 
    mutate(location = as.numeric(location), 
           Area = as.numeric(Area)) %>% 
    slice(match(locationList[[g]], location)) %>% 
    select(location, Sample.Name, Area)
  
  return(tempData)
  
}



# Function to generate the needed coefficients for std transformation
get_corr_coefficients <- function(i, areas, stds_amount){
  
  tempCoeff <- lm(areas[[i]]$Area ~ stds_amount[[i]]$conc)$coefficients
  
  return(tempCoeff)
  
}


# Function to convert the standards based on provided area
run_conversion <- function(i, stds_coeff, std_areas){
  
  tempConc <- (std_areas[[i]] - filter(stds_coeff, variables == "Intercept")[, i]) / 
    filter(stds_coeff, variables == "m")[, i]
  
  return(tempConc)
  
}





# Function to create correction factor for the samples
get_correction_factor <- function(){
  
  
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

# Make a sample only list


area_list <- sapply(names(stds_location), 
               function(x) get_areas(x, std_areas, stds_location), simplify = F)

stand_list <- sapply(names(stds_location), 
                     function(x) get_standards(x, std_areas, stds_location), 
                     simplify = F)

coeff_list <- sapply(names(area_list), 
                     function(x) get_corr_coefficients(x, area_list, stand_list)) %>% 
  as.data.frame() %>% mutate(variables = c("Intercept", "m"))


converted_stds <- sapply(names(area_list), 
                         function(x) run_conversion(x, coeff_list, area_list))


sample_list <- sapply(names(sample_location), 
                      function(x) 
                        suppressWarnings(
                          get_samples(x, sample_location, dataTable)), simplify = F)



















# TO DO LIST
### Create correction based on location run
### Correct concentration based on location of sample relative to three different standards used
    ### e.g. use std 1 and 2 if samples between them and std 2 and 3 if samples between them (mM)
### Convert to mmol/kg based on weight measured and total volume originally suspended in