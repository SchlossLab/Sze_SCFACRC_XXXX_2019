### Graphs of measured SCFAs by HPLC
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
  # scfa_name is a variable that stores the scfa of interest e.g. "acetate"
  # path_to_file stores the directory path where the files are
  # ending is the common ending of the file of interest
  
  # the command that reads in the files
  tempData <- read.csv(paste(path_to_file, scfa_name, ending, sep = ""), 
                       header = T, stringsAsFactors = F)
  # sends the data to the global work environment
  return(tempData)
}


# Function to merge the respective scfa data with the needed metadata for visualization
merge_data <- function(scfa_name, meta_data, dataList, not_enough){
  # scfa_name is a variable that stores the scfa of interest e.g. "acetate"
  # meta_data contains the needed meta data on disease groups, etc.
  # dataList is a list with all the scfa data tables of interests
  # not_enough is a vector with the samples that could not be run on the HPLC
  
  # Create a temp data table variable
  tempData <- dataList[[scfa_name]]
  
  # Modify and merge the data together for samples that matched
  temp_combined_data <- meta_data %>% 
    select(study_id, Dx_Bin, dx) %>% 
    mutate(study_id = as.character(study_id)) %>% 
    inner_join(tempData, by = "study_id") %>% 
    distinct(study_id, .keep_all = TRUE)
  
  # get sample IDs that were measured but had no signal
  # and remove those that were never run on HPLC
  test <- meta_data %>% 
    select(study_id, Dx_Bin, dx) %>% 
    mutate(study_id = as.character(study_id)) %>% 
    filter(!(study_id %in% temp_combined_data$study_id) & 
             !(study_id %in% not_enough)) %>% 
    mutate(mmol_kg = rep(0, length(study_id)))
  # create a final combined table with appropriate label names
  final_temp_combined <- temp_combined_data %>% 
    bind_rows(test) %>% 
    mutate(Dx_Bin = ifelse(Dx_Bin == "High Risk Normal", invisible("Normal"), 
                           ifelse(Dx_Bin == "adv Adenoma", 
                                  invisible("adv_Adenoma"), invisible(Dx_Bin))))
  # Send the final data table to the global work environment
  return(final_temp_combined)
  
}

##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################


# Load in needed metadata and rename sample column
metaI <- read.csv("data/raw/metadata/metaI_final.csv", 
                  header = T, stringsAsFactors = F) %>% 
  rename(study_id = sample)
# read in all the respective scfa data tables
scfa_data <- sapply(scfas, 
                    function(x) upload_scfa_data(x, "data/process/tables/", "_final_data.csv"), 
                    simplify = F)
# combine all the necessary data together into a single file
combined_data <- sapply(scfas, 
                        function(x) merge_data(x, metaI, scfa_data, not_measured), simplify = F)


# Create graph for each SCFA

### Acetate overall values by group
acetate <- combined_data[["acetate"]] %>% 
  ggplot(aes(
    factor(dx, 
           levels = c("normal", "adenoma", "cancer"), 
           labels = c("Control", "Adenoma", "Carcinoma")), mmol_kg)) + 
  geom_jitter(aes(color = dx), width = 0.3, size = 2.5, alpha = 0.75, show.legend = F) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.75, width = 0.5)  + 
  scale_color_manual(values = c('#FFD700', '#DC143C', '#228B22')) + 
  annotate("text", label = paste("Acetate"), x = 3, y = 300, size = 4) + 
  
  xlab("") + ylab("mmol per Kg") + theme_bw() + ggtitle("A") + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))


### Butyrate overall values by group
butyrate <- combined_data[["butyrate"]] %>% 
  ggplot(aes(
    factor(dx, 
           levels = c("normal", "adenoma", "cancer"), 
           labels = c("Control", "Adenoma", "Carcinoma")), mmol_kg)) + 
  geom_jitter(aes(color = dx), width = 0.3, size = 2.5, alpha = 0.75, show.legend = F) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.75, width = 0.5)  + 
  scale_color_manual(values = c('#FFD700', '#DC143C', '#228B22')) + 
  annotate("text", label = paste("Butyrate"), x = 3, y = 80, size = 4) + 
  xlab("") + ylab("mmol per Kg") + theme_bw() + ggtitle("B") + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))


### Isobutyrate overall values by group
isobutyrate <- combined_data[["isobutyrate"]] %>% 
  ggplot(aes(
    factor(dx, 
           levels = c("normal", "adenoma", "cancer"), 
           labels = c("Control", "Adenoma", "Carcinoma")), mmol_kg)) + 
  geom_jitter(aes(color = dx), width = 0.3, size = 2.5, alpha = 0.75, show.legend = F) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.75, width = 0.5)  + 
  scale_color_manual(values = c('#FFD700', '#DC143C', '#228B22')) + 
  annotate("text", label = paste("Isobutyrate"), x = 3, y = 30, size = 4) + 
  xlab("") + ylab("mmol per Kg") + theme_bw() + ggtitle("C") + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))


### Propionate overall values by group
propionate <- combined_data[["propionate"]] %>% 
  ggplot(aes(
    factor(dx, 
           levels = c("normal", "adenoma", "cancer"), 
           labels = c("Control", "Adenoma", "Carcinoma")), mmol_kg)) + 
  geom_jitter(aes(color = dx), width = 0.3, size = 2.5, alpha = 0.75, show.legend = F) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.75, width = 0.5)  + 
  scale_color_manual(values =  c('#FFD700', '#DC143C', '#228B22')) + 
  annotate("text", label = paste("Propionate"), x = 3, y = 110, size = 4) + 
  xlab("") + ylab("mmol per Kg") + theme_bw() + ggtitle("D") + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))

# Combine the graphs together

scfa_plot <- grid.arrange(acetate, butyrate, isobutyrate, propionate)

# Write out to specific directory
ggsave("results/figures/hplc_scfa_graph.tiff", scfa_plot, width = 8, height = 8)







