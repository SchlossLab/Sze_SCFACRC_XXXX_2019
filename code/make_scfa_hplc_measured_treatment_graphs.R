### Graphs of measured SCFAs by HPLC
### Before and after treatment
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
    select(study_id, Dx_Bin, dx, EDRN, time_point) %>% 
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
metaF <- read.csv("data/raw/metadata/good_metaf_final.csv", 
                  header = T, stringsAsFactors = F) %>% 
  gather("time_point", "study_id", initial, followUp) %>% 
  select(study_id, time_point, everything())

# read in all the respective scfa data tables
scfa_data <- sapply(scfas, 
                    function(x) upload_scfa_data(x, "data/process/tables/", "_final_data.csv"), 
                    simplify = F)

# combine all the necessary data together into a single file
combined_data <- sapply(scfas, 
                        function(x) merge_data(x, metaF, scfa_data, not_measured), simplify = F)


acetate_graph <- combined_data[["acetate"]] %>% 
  filter(!is.na(time_point)) %>% 
  mutate(time_point =  factor(time_point, 
                      levels = c("initial", "followUp"), 
                      labels = c("Pre-", "Post-")), 
         dx = factor(dx, levels = c("adenoma", "cancer"), labels = c("Adenoma", "Carcinoma"))) %>% 
  ggplot(aes(time_point, log10(mmol_kg+1), color = dx, group = EDRN)) + 
  geom_point(show.legend = F) + geom_line(show.legend = F) + 
  stat_summary(mapping = aes(group = dx), 
               fun.y = median, colour = "black", geom = "point", size = 3) + 
  stat_summary(mapping = aes(group = dx), 
               fun.y = median, colour = "black", geom = "line", size = 1) + 
  facet_grid(~dx) + coord_cartesian(ylim = c(0, 2.5)) + 
  theme_bw() + labs(x = "", y = expression(Log["10"]~mmol~per~Kg)) + 
  scale_color_manual(values = c('#FFD700', '#DC143C')) + ggtitle("A") + 
  annotate("text", label = paste("Acetate"), x = 1.5, y = 2.5, size = 3) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))
  

butyrate_graph <- combined_data[["butyrate"]] %>% 
  filter(!is.na(time_point)) %>% 
  mutate(time_point =  factor(time_point, 
                              levels = c("initial", "followUp"), 
                              labels = c("Pre-", "Post-")), 
         dx = factor(dx, levels = c("adenoma", "cancer"), labels = c("Adenoma", "Carcinoma"))) %>% 
  ggplot(aes(time_point, log10(mmol_kg+1), color = dx, group = EDRN)) + 
  geom_point(show.legend = F) + geom_line(show.legend = F) + 
  stat_summary(mapping = aes(group = dx), 
               fun.y = median, colour = "black", geom = "point", size = 3) + 
  stat_summary(mapping = aes(group = dx), 
               fun.y = median, colour = "black", geom = "line", size = 1) + 
  facet_grid(~dx) + coord_cartesian(ylim = c(0, 2)) + 
  theme_bw() + labs(x = "", y = expression(Log["10"]~mmol~per~Kg)) + 
  scale_color_manual(values = c('#FFD700', '#DC143C')) + ggtitle("B") + 
  annotate("text", label = paste("Butyrate"), x = 1.5, y = 2, size = 3) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))


isobutyrate_graph <- combined_data[["isobutyrate"]] %>% 
  filter(!is.na(time_point)) %>% 
  mutate(time_point =  factor(time_point, 
                              levels = c("initial", "followUp"), 
                              labels = c("Pre-", "Post-")), 
         dx = factor(dx, levels = c("adenoma", "cancer"), labels = c("Adenoma", "Carcinoma")), 
         mmol_kg = ifelse(mmol_kg < 0, invisible(0), invisible(mmol_kg))) %>% 
  ggplot(aes(time_point, log10(mmol_kg+1), color = dx, group = EDRN)) + 
  geom_point(show.legend = F) + geom_line(show.legend = F) + 
  stat_summary(mapping = aes(group = dx), 
               fun.y = median, colour = "black", geom = "point", size = 3) + 
  stat_summary(mapping = aes(group = dx), 
               fun.y = median, colour = "black", geom = "line", size = 1) + 
  facet_grid(~dx) + coord_cartesian(ylim = c(0, 2)) + 
  theme_bw() + labs(x = "", y = expression(Log["10"]~mmol~per~Kg)) + 
  scale_color_manual(values = c('#FFD700', '#DC143C')) +  ggtitle("C") + 
  annotate("text", label = paste("Isobutyrate"), x = 1.5, y = 2, size = 3) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))


propionate_graph <- combined_data[["propionate"]] %>% 
  filter(!is.na(time_point)) %>% 
  mutate(time_point =  factor(time_point, 
                              levels = c("initial", "followUp"), 
                              labels = c("Pre-", "Post-")), 
         dx = factor(dx, levels = c("adenoma", "cancer"), labels = c("Adenoma", "Carcinoma"))) %>% 
  ggplot(aes(time_point, log10(mmol_kg + 1), color = dx, group = EDRN)) + 
  geom_point(show.legend = F) + geom_line(show.legend = F) + 
  stat_summary(mapping = aes(group = dx), 
               fun.y = median, colour = "black", geom = "point", size = 3) + 
  stat_summary(mapping = aes(group = dx), 
               fun.y = median, colour = "black", geom = "line", size = 1) + 
  facet_grid(~dx) + coord_cartesian(ylim = c(0, 2.5)) + 
  theme_bw() + labs(x = "", y = expression(Log["10"]~mmol~per~Kg)) + 
  scale_color_manual(values = c('#FFD700', '#DC143C')) +  ggtitle("D") + 
  annotate("text", label = paste("Propionate"), x = 1.5, y = 2.5, size = 3) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10, face = "bold"))


# Combine the graphs together

treatment_scfa_plot <- grid.arrange(acetate_graph, butyrate_graph, 
                                    isobutyrate_graph, propionate_graph)

# Write out to specific directory
ggsave("results/figures/hplc_treatment_scfa_graph.tiff", treatment_scfa_plot, width = 8, height = 8)


