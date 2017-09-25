### Prelim analys of butyrate
### cursory look at butyrate and CRC
### Marc Sze



# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "ggplot2", "gridExtra"))

#setup variables that will be used
scfas <- c("acetate", "butyrate", "isobutyrate", "propionate")

# Samples that had to little volume to actually be measured
not_measured <- c("3367653", "2041650", "2043650", "3027650")


##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################

# Function to load in the respective scfa data
upload_scfa_data <- function(scfa_name, path_to_file, ending){
  
  tempData <- read.csv(paste(path_to_file, scfa_name, ending, sep = ""), 
                       header = T, stringsAsFactors = F)
  
  return(tempData)
  
}

# Function to merge the respective scfa data with the needed metadata for visualization
merge_data <- function(scfa_name, meta_data, dataList, not_enough){
  
  tempData <- dataList[[scfa_name]]
  
  # Modify and merge the data together
  
  temp_combined_data <- meta_data %>% 
    select(study_id, Dx_Bin, dx) %>% 
    mutate(study_id = as.character(study_id)) %>% 
    inner_join(tempData, by = "study_id") %>% 
    distinct(study_id, .keep_all = TRUE)
  
  # get sample IDs that were measured but had no signal
  test <- meta_data %>% 
    select(study_id, Dx_Bin, dx) %>% 
    mutate(study_id = as.character(study_id)) %>% 
    filter(!(study_id %in% temp_combined_data$study_id) & 
             !(study_id %in% not_enough)) %>% 
    mutate(mmol_kg = rep(0, length(study_id)))
  
  final_temp_combined <- temp_combined_data %>% 
    bind_rows(test) %>% 
    mutate(Dx_Bin = ifelse(Dx_Bin == "High Risk Normal", invisible("Normal"), 
                           ifelse(Dx_Bin == "adv Adenoma", 
                                  invisible("adv_Adenoma"), invisible(Dx_Bin))))
  
  
  return(final_temp_combined)
  
}

# Function to generate first pass visualizations
create_graphs <- function(scfa_name, dataList){
  
  tempData <- dataList[[scfa_name]]
  
  # Graph the results to visualize if patterns exits 
  tempPlot <- ggplot(tempData, aes(factor(Dx_Bin, 
                                   levels = c("Normal", "Adenoma", "adv_Adenoma", "Cancer"), 
                                   labels = c("Control", "Adenoma", "Advanced\nAdenoma", "Carcinoma")), 
                            mmol_kg)) + 
    geom_jitter(aes(color = Dx_Bin), width = 0.3, show.legend = F) + 
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
                 colour = "black", geom = "crossbar", size = 0.5, width = 0.5)  + 
    xlab("") + ylab("mmol per Kg") + theme_bw()
  
  
  tempPlot2 <- tempData %>% 
    mutate(
      Dx_Bin = factor(Dx_Bin, 
                      levels = c("Normal", "Adenoma", "adv_Adenoma", "Cancer"), 
                      labels = c("Control", "Adenoma", "Advanced\nAdenoma", "Carcinoma"))) %>% 
    ggplot(aes(mmol_kg, group = Dx_Bin, fill = Dx_Bin)) + 
    geom_density(alpha = 0.3, show.legend = F) + 
    xlab("mmol per Kg") + ylab("Density") + theme_bw() + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
    ggtitle(paste(scfa_name))
  
  
  tempFinalGraph <- grid.arrange(tempPlot2, tempPlot, nrow = 2)
  
  return(tempFinalGraph)
  
}


#Run an ANOVA with Tukey post hoc test to look for differences within each group
get_anova_comparisons <- function(scfa_name, dataList, 
                                  set_variable = "mmol_kg", set_groups = "Dx_Bin"){
  # alpha_metric represents the alpha variables of interest
  # set_groups specifies what column to use for the groups variable
  # data_set is the combined data set (used in all testing functions)
  
  tempData <- dataList[[scfa_name]]
  
  # this runs an ANOVA with Tukey post hoc test and makes a readable table 
  test_data <- summary(aov(lm(
    as.formula(paste(set_variable, "~", set_groups)), data = tempData)))[[1]]
  
  return(test_data)
  
}


# Function to save plots as pdfs to be viewed as needed
save_gg_plots <- function(scfa_name, graphLists, path_to_save, ending){
  
  ggsave(paste(path_to_save, scfa_name, ending), graphLists[[scfa_name]], 
         device = "pdf", width = 8, height = 6)
  
  print(paste("completed saving graph for ", scfa_name, sep = ""))
  
}

# Function to save ANOVA results
save_anova_comparisons <- function(scfa_name, dataList, 
                                   path_to_save, ending){
  
  tempData <- dataList[[scfa_name]] %>% 
    rename(pvalue = `Pr(>F)`)
  
  write.csv(tempData, 
            paste(path_to_save, scfa_name, ending, sep = ""), 
            row.names = F)
  
}



##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

# Load in needed data

metaI <- read.csv("data/raw/metadata/metaI_final.csv", 
                  header = T, stringsAsFactors = F) %>% 
  rename(study_id = sample)

scfa_data <- sapply(scfas, 
                    function(x) upload_scfa_data(x, "data/process/tables/", "_final_data.csv"), 
                    simplify = F)

combined_data <- sapply(scfas, 
                        function(x) merge_data(x, metaI, scfa_data, not_measured), simplify = F)

graphs <- sapply(scfas, 
                 function(x) create_graphs(x, combined_data), simplify = F)

test <- sapply(scfas, 
               function(x) get_anova_comparisons(x, combined_data), simplify = F)

lapply(scfas, 
       function(x) 
         save_gg_plots(x, graphs, 
                       "exploratory/notebook/exploratory_graphs/", 
                       "_general_distributions.pdf"))

lapply(scfas, 
       function(x) 
         save_anova_comparisons(x, test, 
                                "data/process/tables/", 
                                "_anova_crc_groups.csv"))









