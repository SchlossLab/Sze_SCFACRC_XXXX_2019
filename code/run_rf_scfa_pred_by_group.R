
##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

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


# Split into normal, adenoma, and cancer groups
split_final_data <- list(
  final_data_normal = sapply(scfas, 
                             function(x) final_data[[x]] %>% 
                               filter(dx == "normal"), simplify = F), 
  final_data_adenoma = sapply(scfas, 
                              function(x) final_data[[x]] %>% 
                                filter(dx == "adenoma"), simplify = F), 
  final_data_carcinoma = sapply(scfas, 
                                function(x) final_data[[x]] %>% 
                                  filter(dx == "cancer"), simplify = F))

# Set up initial store for summary data
class_summary_model_data <- list(
  final_data_normal = list(acetate = data_frame(), butyrate = data_frame(), 
                           isobutyrate = data_frame(), propionate = data_frame()), 
  final_data_adenoma = list(acetate = data_frame(), butyrate = data_frame(), 
                            isobutyrate = data_frame(), propionate = data_frame()), 
    final_data_carcinoma = list(acetate = data_frame(), butyrate = data_frame(), 
                              isobutyrate = data_frame(), propionate = data_frame()))
  
  
class_important_vars <- list(
  final_data_normal = list(acetate = data_frame(), butyrate = data_frame(), 
                           isobutyrate = data_frame(), propionate = data_frame()), 
  final_data_adenoma = list(acetate = data_frame(), butyrate = data_frame(), 
                            isobutyrate = data_frame(), propionate = data_frame()), 
  final_data_carcinoma = list(acetate = data_frame(), butyrate = data_frame(), 
                              isobutyrate = data_frame(), propionate = data_frame()))


# Run a 100 different splits for each group and SCFA on classification
for(i in names(split_final_data)){
  
  print(paste("Getting model data from ", i, " only...."))
  
  for(j in 1:1){
    
    # Generate an 80/20 data split
    rf_data <- sapply(scfas, function(x) eighty_twenty_split(x, "rf_groups", split_final_data[[i]]), simplify = F)
    
    class_test <- sapply(scfas, function(x)
      make_rf_model(x, j, "training_data", "high_low", "rf", "ROC", rf_data),
      simplify = F)

    class_pred <- sapply(scfas, function(x)
      run_prediction(x, class_test, "test_data", "raw", "high_low", rf_data),
      simplify = F)

    ROC_data <- sapply(scfas, function(x)
      get_test_roc(x, "high_low", "tempPredictions", class_pred, classif = T), simplify = F)

    class_summary_model_data <- sapply(scfas, function(x)
      add_model_summary_data(x, class_test, ROC_data, class_summary_model_data[[i]], classif = T), simplify = F)

    class_important_vars <- sapply(scfas, function(x)
      grab_importance(x, j, "classification", class_test, class_important_vars[[i]]), simplify = F)
  }
  
  print(paste("Completed getting all model data for ", i, " only."))
}






