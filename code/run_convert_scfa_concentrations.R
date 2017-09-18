### Creation of function to make the conversions
### Need to adjust readings to get the final SCFA concentrations
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "ggplot2"))



dataTable <- read.csv("exploratory/scratch/reference_for_calcs/test_data.csv", header = T, stringsAsFactors = F)

g1_stds <- c(20, 10, 5, 2.5, 1.0, 0.5, 0.25, 0.1)

stds1_location <- c(2:9)
stds2_location <- c(43:49)
stds3_location <- c(79:85)


std_areas <- dataTable %>% filter(grepl("STD", Sample.ID) == TRUE & grepl("biorad", Sample.ID) != TRUE) %>% 
  separate(Data.Filename, c("temp1", "temp2"), sep = "_") %>% 
  rename(location = temp2) %>% mutate(location = as.numeric(gsub("\\.lcd", "", location)))

temp_areas <- std_areas %>% slice(match(stds1_location, location)) %>% select(Area) %>% 
  mutate(Area = as.numeric(Area))


test1 <- lm(temp_areas$Area ~ g1_stds)

std1_corr_factors <- test1$coefficients



# TO DO LIST
### Create three different linear regressions
    ### One for each standard (x = concentration, y = Area)
### Convert each standard based on unique model obtained
### Create correction based on location run
### Correct concentration based on location of sample relative to three different standards used
    ### e.g. use std 1 and 2 if samples between them and std 2 and 3 if samples between them (mM)
### Convert to mmol/kg based on weight measured and total volume originally suspended in