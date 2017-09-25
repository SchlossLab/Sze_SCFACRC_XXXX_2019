### Prelim analys of butyrate
### cursory look at butyrate and CRC
### Marc Sze



# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "ggplot2"))


# Load in needed data

metaI <- read.csv("data/raw/metadata/metaI_final.csv", 
                  header = T, stringsAsFactors = F) %>% 
  rename(study_id = sample)

butyrate_data <- read.csv("data/process/tables/butyrate_final_data.csv", 
                          header = T, stringsAsFactors = F)

# Modify and merge the data together

combined_data <- metaI %>% 
  select(study_id, Dx_Bin, dx) %>% 
  mutate(Dx_Bin = ifelse(Dx_Bin == "High Risk Normal", invisible("Normal"), 
                         ifelse(Dx_Bin == "adv Adenoma", 
                                invisible("adv_Adenoma"), invisible(Dx_Bin))), 
         study_id = as.character(study_id)) %>% 
  inner_join(butyrate_data, by = "study_id") %>% 
  distinct(study_id, .keep_all = TRUE)



ggplot(combined_data, aes(factor(Dx_Bin, 
                                 levels = c("Normal", "Adenoma", "adv_Adenoma", "Cancer"), 
                                 labels = c("Control", "Adenoma", "Advanced\nAdenoma", "Carcinoma")), 
                                 mmol_kg)) + 
  geom_jitter(aes(color = Dx_Bin), width = 0.3, show.legend = F) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5)  + 
  xlab("") + ylab("mmol per Kg") + theme_bw()













