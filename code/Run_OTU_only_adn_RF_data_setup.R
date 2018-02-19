### Build the best adenoma model possible
### Generate test and random 80/20 split data needed for OTU only adn RF analysis
## Marc Sze

#Load needed libraries
source('code/functions.R')

loadLibs(c("tidyverse", "caret", "scales"))

# Read in necessary data frames

shared <- read.delim('data/process/final.0.03.subsample.shared', 
                     header=T, sep='\t') %>% select(Group, contains("Otu0"))

metaI <- read.csv("data/raw/metadata/metaI_final.csv", 
                  stringsAsFactors = F, header = T) 

good_metaf <- read.csv('data/raw/metadata/good_metaf_final.csv', 
                       header = T, stringsAsFactors = F) %>% select(initial)

metaI <- filter(metaI, !(sample %in% good_metaf$initial)) %>% 
  mutate(Dx_Bin = gsub("adv Adenoma", "adv_adenoma", Dx_Bin))



#################################################################################
#                                                                               #
#                                                                               #
#               Data Clean up for better and faster modeling                    #
#                                                                               #
#################################################################################

# Remove follow up samples and join metadata, scfa, and 16S data
test_data <- metaI %>% 
  select(sample, dx) %>% 
  inner_join(shared, by = c("sample" = "Group")) %>% 
  filter(dx != "cancer") %>% as.data.frame() %>% 
  select(-sample) %>% 
  select(dx, everything())


#Filter out rows that are not all complete
test_data <- test_data[complete.cases(test_data), ]

stored_cases_controls <- test_data$dx

# Remove those with near zero variance 
# Similar to removal of singletons or removal based on percentage in sample
# Using near zero variance should cover the above two so could use only this instead

nzv <- nearZeroVar(test_data)
# By default, nearZeroVar will return the positions of the variables that are flagged to be problematic
# using the saveMetrics argument will output a table with more in-depth info
# Group stays in this case (because it is numberic) 
# but maybe better in future to transfer this to row names if not a numberic label

test_data <- test_data[, -nzv]


# Write data table for future use
write.csv(test_data, "data/process/tables/adn_OTU_only_test_data.csv", row.names = F)


#################################################################################
#                                                                               #
#                                                                               #
#               Data Splitting to create train set                              #
#                                                                               #
#################################################################################

# Split data evenly 
set.seed(3457)
eighty_twenty_splits <- createDataPartition(test_data$dx, 
                                            p = 0.8, list = FALSE, times = 100)

# write out the random splits
write.csv(eighty_twenty_splits, "data/process/tables/adn_OTU_only_eighty_twenty_splits.csv", row.names = F)
