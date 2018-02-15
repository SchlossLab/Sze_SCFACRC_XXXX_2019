### Analysis of Relative Abundance of previous Top 10 important variables
### from previous RF models
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "caret"))


# Previous research top 10 (PMID 29145893)

# Adenoma: Lachnospiraceae, Clostridiales, Blautia, Bacteroides, Odoribacter, Ruminococcaceae, Flavonifractor, 
# Roseburia, Escherichia/Shigella, Clostridium XIVa

# Advanced Adenoma: Clostridiales, Clostridium XIVa, Roseburia, Odoribacter, Bacteroides, Ruminococcaceae, Streptococcus, 
# Blautia, Lachnospiraceae, Ruminococcus

# Carcinoma: Porphyromonas, Parvimonas, Fusobacterium, Gemella, Prevotella, Streptococcus, Coprobacillus, Pasteurellaceae, 
# Collinsella, Bilophila

lowest_txID <- list(adn = c("Lachnospiraceae", "Clostridiales", "Blautia", "Bacteroides", 
                            "Odoribacter", "Ruminococcaceae", "Flavonifractor", 
                            "Roseburia", "Escherichia/Shigella", "Clostridium_XlVa", "Anaerostipes", "Streptococcus", 
                            "Ruminococcus", "Clostridium_IV", "Faecalibacterium", "Bifidobacterium", "Enterococcus", 
                            "Clostridium_XI", "Coriobacteriaceae", "Coprobacillus", "Actinomyces", "Alistipes", 
                            "Dorea", "Clostridium_XlVb", "Lactococcus", "Firmicutes", "Clostridium_XVIII", "Gemella", 
                            "Collinsella", "Akkermansia", "Eggerthella", "Coprococcus"), 
                    crc = c("Porphyromonas", "Parvimonas", "Fusobacterium", "Gemella", "Prevotella", "Streptococcus", 
                            "Coprobacillus", "Pasteurellaceae", "Collinsella", "Bilophila", "Bacteroides", "Parabacteroides", 
                            "Clostridium_XlVa", "Odoribacter", "Anaerostipes", "Clostridium_XlVb", "Dorea", "Ruminococcaceae", 
                            "Lachnospiraceae", "Alistipes", "Streptococcus", "Enterobacteriaceae", "Clostridium_XI", 
                            "Blautia", "Clostridium_IV", "Ruminococcus", "Clostridiales"))


combined_meta_data <- read_csv("data/raw/metadata/metaI_final.csv") %>% 
  mutate(Dx_Bin = ifelse(Dx_Bin == "High Risk Normal" | Dx_Bin == "Normal", invisible("control"), 
                         ifelse(Dx_Bin == "Adenoma", invisible("adenoma"), 
                                ifelse(Dx_Bin == "adv Adenoma", invisible("adenoma"), 
                                       ifelse(Dx_Bin == "Cancer", invisible("cancer"), invisible(Dx_Bin)))))) %>% 
  select(sample, Dx_Bin, dx) %>% 
  bind_rows(read_csv("data/raw/metadata/good_metaf_final.csv") %>% 
              mutate(sample = initial, 
                     Dx_Bin = ifelse(Dx_Bin == "adv_adenoma", invisible("adenoma"), invisible(Dx_Bin))) %>% 
              select(sample, Dx_Bin, dx))
  

##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################

# Function to get a subsamples genus file
get_genera_subsample <- function(run_number, dataList){
  
  #tempData <- as.matrix(dataList)
  tempData <- dataList
  
  lowest_seq_count <- min(rowSums(tempData))
  
  total_genera <- length(colnames(tempData))
  
  genera_names <- colnames(tempData)
  
  stored_draws_List <- NULL
  
  # Iteratres through each sample in data set
  for(j in rownames(tempData)){
    
    tempdraw <- c()
    tempVector <- unname(tempData[j, ])
    
    # creates a temp vector with each genus (as a number)
    # repeated based on the number of counts in data frame
    for(k in 1:total_genera){
      
      tempdraw <- c(tempdraw, rep(k, tempVector[k]))
      
    }
    # saves the created data in a new list
    stored_draws_List[[j]] <- tempdraw
  }
  
  # Applies a randome sampling accross the store vector of repeats
  stored_rd <- lapply(stored_draws_List, function(x) sample(x, lowest_seq_count))
  # Generates the counts from the sampling
  num_counts <- lapply(stored_rd, function(x) as.data.frame(table(x), stringsAsFactors = FALSE))
  # Runs the assign genera function
  agg_genera <- lapply(num_counts, 
                       function(x) assign_genera(x, genera_names))
  # Converts the output to a matrix of the same orientation as the data files
  # in the inputted data list
  final_table <- t(as.data.frame.list(agg_genera))
  
  print(paste("Completed", run_number, "sampling"))
  
  # Returns the data
  return(final_table)
}


# Function to get the assignments needed from the sampling
assign_genera <- function(dataTable, generaVector){
  # dataTable is part of a list where each dataTable is an individual sample
  # generaVector is a vector of genera names
  
  # Creates a temporary variable of all taxa with 0 and names the vector
  tempVector <- rep(0, length(generaVector))
  names(tempVector) <- generaVector
  
  # Iteratres through each genera sampled and changes the 0 to the
  # correct number of counts
  updatedVector <- grab_values(tempVector, dataTable)
  
  
  # returns the final vector
  return(updatedVector)  
}


# Function to pull specific value
grab_values <- function(vec_of_int, refTable){
  
  vec_of_int[as.numeric(refTable[, "x"])] <- refTable[, "Freq"]
  
  return(vec_of_int)
  
}


# Function to run the sampling x number of times and generate an average from this
get_average_counts <- function(i, repeats, dataList){
  
  total_samples <- length(rownames(dataList))
  genera_names <- colnames(dataList)
  
  full_100_runs <- lapply(1:repeats, function(x) get_genera_subsample(x, dataList))
  
  
  temp_avg_list <- lapply(c(1:total_samples), 
                          function(x) grab_row(full_100_runs, x, i, genera_names))
  
  
  final_avg_data <- t(as.data.frame.list(temp_avg_list))
  rownames(final_avg_data) <- rownames(dataList)
  
  print(paste("Completed study ", i, ": taxa subsampling.", sep = ""))
  
  return(final_avg_data)
  
}


# Function to grab rows for averaging 
grab_row <- function(list_of_int, j, study, genera_file){
  
  test <- lapply(list_of_int, function(x) x[j, ])
  
  test <- t(as.data.frame.list(test))
  
  colnames(test) <- genera_file
  rownames(test) <- c(1:length(rownames(test)))
  
  average_vector <- colMeans(test)
  
  return(average_vector)
}


get_merge_table <- function(tumor_type, metafile, genus_sharedfile, select_taxa_list){
  
  tumor_names <- list(adn = "adenoma", adv_adn = "adv_adenoma", crc = "cancer")
  
  tempData <- genus_sharedfile[[tumor_type]] %>% mutate(Group = as.numeric(Group))
  
  
  temp_df <- metafile %>% select(sample, Dx_Bin) %>% 
    left_join(tempData, by = c("sample" = "Group")) %>% 
    filter(Dx_Bin == "control" | Dx_Bin == tumor_names[[tumor_type]]) %>% 
    select(sample, Dx_Bin, one_of(select_taxa_list[[tumor_type]]))
  
  return(temp_df)
}



# Function to run the wilcoxson rank sum test on each provided OTU
get_wilcox_testing <- function(tumor_type, test_data){
  
  tempData <- as.data.frame(test_data[[tumor_type]])

  tumor_names <- list(adn = "adenoma", adv_adn = "adv_adenoma", crc = "cancer")
  
  vars_to_test <- colnames(select(test_data[[tumor_type]], -sample, -Dx_Bin))
  
  test <- t(sapply(vars_to_test, 
         function(x) wilcox.test(
           filter(tempData, Dx_Bin == "control")[, x], 
           filter(tempData, Dx_Bin == tumor_names[[tumor_type]])[, x], 
           alternative = "two.sided")$p.value, simplify = F) %>% bind_rows()) %>% 
    as.data.frame() %>% 
    mutate(taxa = rownames(.), 
           bh = p.adjust(V1, method = "BH")) %>% 
    rename(pvalue = V1) %>% 
    arrange(pvalue) %>% 
    select(taxa, pvalue, bh)

  
  
  return(test)
}



##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

# Create genera table
genera_data <- get_tax_level_shared("baxter", read_tsv("data/process/final.shared"), 
                                    read_tsv("data/process/final.taxonomy"), 6)

rownames(genera_data) <- genera_data$Group
genera_data <- genera_data[, -1]


avg_subsample_table <- get_average_counts("baxter", 100, genera_data)

avg_subsample_table <- avg_subsample_table %>% 
  as.data.frame() %>% 
  mutate(Group = rownames(genera_data))


# generate vector to make sapply run
tumors <- c("adn", "crc")
# Grabs the specific taxa in the top 10 of the previous models
specific_taxa <- sapply(tumors, 
                        function(x) 
                          avg_subsample_table %>% select(Group, one_of(lowest_txID[[x]])), simplify = F)

genera_test_tables_list <- sapply(tumors, 
                                  function(x) get_merge_table(x, combined_meta_data, 
                                                              specific_taxa, lowest_txID), simplify = F) 


# Runs a standard wilcoxson rank sum test and a BH correction 
test_result_tables <- sapply(tumors, 
                             function(x) get_wilcox_testing(x, genera_test_tables_list), simplify = F)

# Write out the data
sapply(tumors, 
       function(x) write_csv(test_result_tables[[x]], 
                             paste("data/process/tables/", x, "_16S_top10percent_RF_taxa_testing.csv", sep = "")))





