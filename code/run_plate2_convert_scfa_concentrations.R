### Creation of function to make the conversions
### Need to adjust readings to get the final SCFA concentrations
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "ggplot2"))

# Set up the locations of the standards
stds_positions <- list(std1 = c(1:9), std2 = c(61:69))

# Set up the locastions of the samples
sample_positions <- list(g1 = c(11:60))

# Set up the scfas of interest
scfas <- c("acetate", "butyrate", "isobutyrate", "propionate")

# Set up the plate name
plate <- "plate2_scfa_crc"

# Set up mapping file nam
map_name <- "scfa_plate_metadata.csv"

# Set up raw text file path
raw_path <- "data/raw/Raw_hplc_files/"

# Set up plate mapping path
platemap_path <- "data/raw/metadata/"

# Set up final data name differentiator
final_name <- "transformed"



##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################

# Function to read in data tables for scfa or mapping file
get_scfa_data <- function(scfa_name, uniq_name, ending, path_to_file){
  
  if(!is.null(scfa_name)){
    
    tempData <- read.delim(
      paste(path_to_file, scfa_name, "/", uniq_name, "_", ending, sep = ""), 
      skip = 2, header = T, stringsAsFactors = F, row.names = 1) %>% 
      select(Data.Filename, Sample.ID, Sample.Name, Ret..Time, 
             Area, Peak.Start, Peak.End, Conc.) %>% 
      mutate(
        Area = ifelse(Area == "-----", invisible(0), invisible(Area)), 
        Ret..Time = as.numeric(Ret..Time), 
        Area = as.numeric(Area), 
        Peak.Start = as.numeric(Peak.Start), 
        Peak.End = as.numeric(Peak.End), 
        Conc. = as.numeric(Conc.))
  } else{
    
    tempData <- read.csv(paste(path_to_file, "/", uniq_name, sep = ""), 
                         header = T, stringsAsFactors = F)
  }
  
 
  
  
  
  
  return(tempData)
  
  
}


# Function to remove outliers from the samples. 
# Must be contained in min peak start and max peak end. 
# Allow for 5% variation on either end)
perform_peak_check <- function(dataTable){
  
  tempStds <- dataTable %>% filter(grepl("mM", Sample.Name) == TRUE)
  
  min_start <- min(tempStds$Peak.Start)*0.95
  max_end <- max(tempStds$Peak.End)*1.05
  
  tempData <- dataTable %>% filter(Peak.Start >= min_start, Peak.End <= max_end)
  
  return(tempData)
  
}

# Function to join plate mapping with actual samples
combine_data <- function(peak_check_data, metadata, combine_by){
  
  tempData <- peak_check_data %>% 
    inner_join(metadata, by = combine_by) %>% 
    select(-Sample.Name)
  
  return(tempData)
  
}


# Function to make a standards only table
make_standard_table <- function(peak_check_data){
  
  # Make a standard only table
  tempData <- peak_check_data %>% 
    filter(grepl("mM", Sample.Name) == TRUE) %>% 
    separate(Data.Filename, c("temp1", "temp2"), sep = "_") %>% 
    rename(location = temp2) %>% 
    separate(Sample.Name, c("conc", "unit"), sep = " ") %>% 
    mutate(location = as.numeric(gsub("\\.lcd", "", location)), 
           conc = as.numeric(conc))
  
  return(tempData)
}



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
    slice(match(locationList[[g]], location))
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
                           stds_coeff, std_list_length = NULL, stds = T){
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
  if(as.numeric(gsub("std", "", value1)) <= length(conv_sampleData)){
    # Generate the corrected concentrations
    tempConc <- tempData[, value1] + 
      correctionFactors[[list_number]]*(tempData[, value2] - tempData[, value1])
    # merge the corrected concentrations with location and sample name
    finalTable <- tempData %>% mutate(conc_mM = tempConc)
    
  } else {
    # sets table to NULL if standards are at the end
    finalTable <- NULL
  }
  # prints the data to the global working environment
  return(finalTable)
  
}

# Function to generate final mmol/kg values
# 1000uL  = 1mL = 0.001L
# g needs to be multiplied by 1000 
# mmol/kg = (corr_conc * 0.1) / (g * 1000)
get_mmol_kg <- function(corrected_data){
  
  tempData <- corrected_data %>% 
    mutate(mmol_kg = (conc_mM * (susp_volume_uL / 1E6)) / (stool_weight_g / 1000))
  
  return(tempData)
  
}




# Control function to make the finalized corrected concentration data
get_final_concentrations <- function(scfa_of_interest, rawdataList, metadata, 
                                     stds_location, sample_location){
  
  tempData <- rawdataList[[scfa_of_interest]]
  ptest <- perform_peak_check(tempData)
  
  comb_test <- combine_data(ptest, metadata, "Sample.ID")
  std_areas <- make_standard_table(ptest)
  
  
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
                           function(x) run_conversion(x, area_list, coeff_list))
  
  # generate sample groupings based on location on run
  sample_list <- sapply(names(sample_location), 
                        function(x) 
                          suppressWarnings(
                            get_samples(x, sample_location, comb_test)), simplify = F)
  # generate correction factors based on location on run
  correction_factor_list <- sapply(names(sample_location), 
                                   function(x) get_correction_factor(x, sample_list), simplify = F)
  # get the converted sample values based on location of standards run and location on run
  converted_samples <- sapply(c(1:length(sample_location)), 
                              function(x) 
                                run_conversion(x, sample_list, coeff_list, 
                                               std_list_length = length(stds_location), 
                                               stds = F), simplify = F)
  # get the final corrected concentrations
  corrected_conc <- sapply(c(1:length(sample_location)), 
                           function(x) 
                             get_corr_conc(x, converted_samples, correction_factor_list), 
                           simplify = F) %>% bind_rows()
  
  finalData <- get_mmol_kg(corrected_conc)
  
  return(finalData)
}


# Create a function that will write out all the data
write_the_data <- function(scfa_of_interest, dataList, 
                           raw_path, plate, final_name){
  
  write.csv(dataList[[scfa_of_interest]], 
            paste(raw_path, scfa_of_interest, "/", final_name, "_", 
                  plate, "_", scfa_of_interest, ".csv", sep = ""), row.names = F)
  
}




##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################


test <- sapply(scfas, 
               function(x) 
                 get_scfa_data(x, plate, paste(x, ".txt", sep = ""), raw_path), 
                               simplify = F)
  
meta_test <- get_scfa_data(NULL, map_name, NULL, platemap_path)


final_test <- sapply(scfas, 
                     function(x) 
                       get_final_concentrations(x, test, meta_test, 
                                                stds_positions, sample_positions), simplify = F)
  
sapply(scfas, 
       function(x) write_the_data(x, final_test, 
                                  raw_path = raw_path, 
                                  plate = plate, 
                                  final_name = final_name))




