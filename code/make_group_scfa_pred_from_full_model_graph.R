### Graph the aggregated data
### Visualize how well do the predictions do by group (normal, adenoma, cancer)?
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse"))

#setup variables that will be used
scfas <- c("acetate", "butyrate", "propionate")


# Read in the data
# The data is obtained from the model tracking of the run_rf_scfa_predictions.R script
# Put in a temporary run number until the analysis is re-run with it in the data frame.

class_data <- map(scfas, function(x) 
  read_csv(paste("data/process/tables/", x, "_classification_RF_train_probs_summary.csv", sep = "")) %>% 
    mutate(scfa = x, 
           model = "class_train")) %>% bind_rows() %>% 
  bind_rows(map(scfas, function(x) 
    read_csv(paste("data/process/tables/", x, "_classification_RF_test_probs_summary.csv", sep = "")) %>% 
      mutate(run = rep(1:100, each = 85), 
             scfa = x, 
             model = "class_test")) %>% 
      bind_rows() %>% 
      rename(pred = tempPredictions, obs = high_low, sample_id = Group)) %>% 
  mutate(correct_class = case_when(
    pred == obs ~ "yes", 
    TRUE ~ "no")) %>% 
  group_by(model, scfa, run, dx) %>% 
  summarise(yes = table(correct_class)[2]/(table(correct_class)[2] + table(correct_class)[1]), 
            no = table(correct_class)[1]/(table(correct_class)[2] + table(correct_class)[1])) %>% 
  ungroup() %>% 
  group_by(model, scfa, dx) %>% 
  summarise(median_yes = median(yes, na.rm = T), 
            median_no = median(no, na.rm = T), 
            min_yes = min(yes, na.rm = T), 
            max_yes = max(yes, na.rm = T))

  
  
reg_data <- map(scfas, function(x) 
  read_csv(paste("data/process/tables/", x, "_regression_RF_train_conc_summary.csv", sep = "")) %>% 
  mutate(scfa = x, 
         model = "reg_train")) %>% 
  bind_rows() %>% 
  bind_rows(map(scfas, function(x) 
    read_csv(paste("data/process/tables/", x, "_regression_RF_test_conc_summary.csv", sep = "")) %>% 
      mutate(run = rep(1:100, each = 85), 
             scfa = x, 
             model = "reg_test")) %>% 
      bind_rows() %>% 
      rename(obs = mmol_kg, pred = tempPredictions, sample_id = Group)) %>% 
  mutate(difference = obs - pred) %>% 
  group_by(model, scfa, run, dx) %>% 
  summarise(median_diff = median(difference), 
            min_dff = min(difference), 
            max_diff = max(difference)) %>% 
  ungroup()

log_reg_data <- map(scfas, function(x) 
  read_csv(paste("data/process/tables/log_", x, "_regression_RF_train_conc_summary.csv", sep = "")) %>% 
    mutate(scfa = x, 
           model = "reg_train")) %>% 
  bind_rows() %>% 
  bind_rows(map(scfas, function(x) 
    read_csv(paste("data/process/tables/log_", x, "_regression_RF_test_conc_summary.csv", sep = "")) %>% 
      mutate(run = rep(1:100, each = 85), 
             scfa = x, 
             model = "reg_test")) %>% 
      bind_rows() %>% 
      rename(obs = mmol_kg, pred = tempPredictions, sample_id = Group)) %>% 
  mutate(difference = obs - pred) %>% 
  group_by(model, scfa, run, dx) %>% 
  summarise(median_diff = median(difference), 
            min_dff = min(difference), 
            max_diff = max(difference)) %>% 
  ungroup()

### Create the classification graph







