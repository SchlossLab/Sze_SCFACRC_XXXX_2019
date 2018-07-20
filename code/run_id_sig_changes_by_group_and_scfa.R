### ID largest changers based on SCFA group or concentrations
### Group SCFA analysis
# Marc Sze



# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "caret"))


#setup variables that will be used
scfas <- c("acetate", "butyrate", "propionate")

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
  tempData <- read_csv(paste(path_to_file, scfa_name, ending, sep = "")) %>% 
    rename(Group = study_id)
  # sends the data to the global work environment
  return(tempData)
}

# Function to remove specific data from sets of interest
pare_down_data <- function(ids_to_remove, scfa_list, otu_data, meta_data){
  
  # Remove data from tables first
  temp_otuData <- otu_data %>% filter(!(Group %in% ids_to_remove))
  
  nzv <- nearZeroVar(temp_otuData)
  
  nzv_rm_tempOTU <- temp_otuData[, -nzv]
  
  temp_meta_data <- meta_data %>% filter(!(Group %in% ids_to_remove))
  
  # Remove data from lists
  tempList <- lapply(scfa_list, 
                     function(x) filter(x, !(Group %in% ids_to_remove)) %>% 
                       inner_join(temp_meta_data, by = "Group") %>% 
                       inner_join(nzv_rm_tempOTU, by = "Group")) 
  
  return(tempList)
}

# Function to generate the high low groupings for each respective SCFA
get_high_low <- function(all_data_list, scfa_median_list, scfas){
  
  tempList <- sapply(scfas, 
                     function(x) all_data_list[[x]] %>% 
                       mutate(high_low = ifelse(mmol_kg <= scfa_median_list[[x]], 
                                                invisible("low"), invisible("high"))) %>% 
                       select(Group, mmol_kg, high_low, Dx_Bin, dx, fit_result, everything()), 
                     simplify = F)
  
  return(tempList)
  
}

# Function to create relevent RF data frames
split_dataframe <- function(i, dataList){
  
  tempData <- dataList[[i]]
  
  temp_regression <- tempData %>% select(mmol_kg, Group, dx, contains("Otu"))
  
  temp_high_low <- tempData %>% select(high_low, Group, dx, contains("Otu"))
  
  finalList <- list(rf_regression = temp_regression, 
                    rf_groups = temp_high_low)
  return(finalList)
}



##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

taxonomy <- read_tsv("data/process/final.taxonomy") %>% 
  select(-Size) %>% 
  mutate(Taxonomy = str_replace_all(Taxonomy, "\\(\\d*\\)", "")) %>% 
  separate(Taxonomy, c("kingdom", "phyla", "class", "order", "family", "genus", "species"), sep = ";")

# Load in needed metadata and rename sample column
picrust_meta <- read_tsv("data/process/picrust_metadata") %>% 
  mutate(Group = as.character(Group))

# Load in more specific follow up data
metaF <- read_csv("data/raw/metadata/good_metaf_final.csv")

# Load in shared data
shared <- read_tsv("data/process/final.0.03.subsample.shared") %>% 
  mutate(Group = as.character(Group)) %>% 
  select(-label, -numOtus)


# Generate IDs that are initial and follow up
initial_IDs <- as.character(metaF$initial)
follow_IDs <- as.character(metaF$followUp)

# read in all the respective scfa data tables
scfa_data <- sapply(scfas, 
                    function(x) upload_scfa_data(x, "data/process/tables/", "_final_data.csv"), 
                    simplify = F)

# Combine all data into a single list of lists and remove follow samples
all_data <- pare_down_data(c(initial_IDs, follow_IDs), scfa_data, shared, picrust_meta)


# Get medians for each scfa measured
scfa_medians <- lapply(all_data, function(x) median(x$mmol_kg))

# get final data tables to be used in the RF analysis
final_data <- get_high_low(all_data, scfa_medians, scfas)

#split the data for the downstream nzv that needs to occur
final_sp_data <- sapply(scfas, 
                        function(x) split_dataframe(x, final_data), simplify = F)

# Set up the separate analysis lists
class_data <- list(acetate = c(), butyrate = c(), propionate = c())
reg_data <- list(acetate = c(), butyrate = c(), propionate = c())

norm_class <- list(acetate = c(), butyrate = c(), propionate = c())
adn_class <- list(acetate = c(), butyrate = c(), propionate = c())
crc_class <- list(acetate = c(), butyrate = c(), propionate = c())

norm_reg <- list(acetate = c(), butyrate = c(), propionate = c())
adn_reg <- list(acetate = c(), butyrate = c(), propionate = c())
crc_reg <- list(acetate = c(), butyrate = c(), propionate = c())

# Set up and run the analysis
for(i in scfas){
  
  # Generate the overall comparisons and pvalues for each scfa by group
  class_data[[i]] <- final_sp_data[[i]]$rf_groups %>%
    gather("otu", "rel_abund", contains("Otu")) %>%
    group_by(otu, dx) %>%
    nest() %>%
    mutate(wil_test = map(data, ~wilcox.test(rel_abund ~ high_low, data = .x)),
           summary_data = map(wil_test, broom::tidy)) %>%
    select(otu, dx, summary_data) %>%
    unnest(summary_data) %>% 
    left_join(select(taxonomy, OTU, genus), by = c("otu" = "OTU")) %>% 
    mutate(scfa = i)
  # Separate by normal
  norm_class[[i]] <- class_data[[i]] %>%
    filter(dx == "normal") %>%
    mutate(bh = p.adjust(p.value, method = "BH")) %>%
    arrange(bh)
  # Separate by adenoma
  adn_class[[i]] <- class_data[[i]] %>%
    filter(dx == "adenoma") %>%
    mutate(bh = p.adjust(p.value, method = "BH")) %>%
    arrange(bh)
  # Separate by carcinoma
  crc_class[[i]] <- class_data[[i]] %>%
    filter(dx == "cancer") %>%
    mutate(bh = p.adjust(p.value, method = "BH")) %>%
    arrange(bh)
  # Status message
  print(paste("Completed classification high/low group tests of ", i, ".", sep = ""))
  
  # Generate the overall correlations (spearman) and pvalues for all scfas by group
  reg_data[[i]] <- final_sp_data[[i]]$rf_regression %>% 
    gather("otu", "rel_abund", contains("Otu")) %>% 
    group_by(otu, dx) %>% 
    nest() %>% 
    mutate(spear_test = map(data, ~cor.test(.x$rel_abund, .x$mmol_kg, method = "spearman")), 
           summary_data = map(spear_test, broom::tidy)) %>% 
    select(otu, dx, summary_data) %>% 
    unnest(summary_data) %>% 
    left_join(select(taxonomy, OTU, genus), by = c("otu" = "OTU")) %>% 
    mutate(scfa = i)
  # Separate by normal
  norm_reg[[i]] <- reg_data[[i]] %>%
    filter(dx == "normal") %>%
    mutate(bh = p.adjust(p.value, method = "BH")) %>%
    arrange(bh)
  # Separate by adenoma
  adn_reg[[i]] <- reg_data[[i]] %>%
    filter(dx == "adenoma") %>%
    mutate(bh = p.adjust(p.value, method = "BH")) %>%
    arrange(bh)
  #Separate by carcinoma
  crc_reg[[i]] <- reg_data[[i]] %>%
    filter(dx == "cancer") %>%
    mutate(bh = p.adjust(p.value, method = "BH")) %>%
    arrange(bh)
  # Status message
  print(paste("Completed correlation of mmol/kg tests of ", i, ".", sep = ""))
  
}


sig_class_summary <- norm_class %>% bind_rows() %>% 
  filter(bh < 0.05) %>% 
  mutate(model = "norm_class") %>% 
  select(otu, genus, dx, statistic, p.value, bh, scfa, model) %>% 
  bind_rows(adn_class %>% bind_rows() %>% 
              filter(bh < 0.05) %>% 
              mutate(model = "adn_class") %>% 
              select(otu, genus, dx, statistic, p.value, bh, scfa, model), 
            crc_class %>% bind_rows() %>% 
              filter(bh < 0.05) %>% 
              mutate(model = "crc_class") %>% 
              select(otu, genus, dx, statistic, p.value, bh, scfa, model))

sig_reg_summary <- norm_reg %>% bind_rows() %>% 
  filter(bh < 0.05) %>% 
  mutate(model = "norm_reg") %>% 
  select(otu, genus, dx, statistic, estimate, p.value, bh, scfa, model) %>% 
  bind_rows(adn_reg %>% bind_rows() %>% 
              filter(bh < 0.05) %>% 
              mutate(model = "adn_reg") %>% 
              select(otu, genus, dx, statistic, estimate, p.value, bh, scfa, model), 
            crc_reg %>% bind_rows() %>% 
              filter(bh < 0.05) %>% 
              mutate(model = "crc_reg") %>% 
              select(otu, genus, dx, statistic, estimate, p.value, bh, scfa, model))



