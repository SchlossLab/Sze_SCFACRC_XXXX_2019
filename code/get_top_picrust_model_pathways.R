### Identify Important Pathways in picrust model
### Pull top pathways 
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

source("https://bioconductor.org/biocLite.R")

library(KEGGREST)

# Load needed libraries
loadLibs(c("tidyverse", "stringr"))

# Load in needed kegg pathway data
kegg_pathway <- read_csv("data/process/tables/kegg_id_key.csv") 
# Load in model information
adn_model_imp <- read_csv("data/process/tables/adenoma_imp_otus_classification_RF_summary.csv")
crc_model_imp <- read_csv("data/process/tables/cancer_imp_otus_classification_RF_summary.csv")


##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################

# Function to summarize all the Keggs by median and quantiles and order highest to lowest
group_stats <- function(dataTable){
  
  tempData <- get(dataTable)
  
  tempData <- tempData %>% 
    group_by(kegg_id) %>% 
    summarise(median_MDA = median(Overall), 
              quantile_25 = quantile(Overall, probs = 0.25), 
              quantile_75 = quantile(Overall, probs = 0.75)) %>% 
    arrange(desc(median_MDA))
  
  return(tempData)
}


# Function to ID the top X pathways
top_model_pathways <- function(dataTable, mapFile, topNumber){
  
  tempData <- dataTable %>% 
    slice(1:topNumber) %>% 
    left_join(mapFile, by = "kegg_id")
  
  return(tempData)
}


# Function to pull specific gene names, pathway, and kegg ortholog
get_gene_ids <- function(dataTable, column_w_keggs){
  
  tempIDs <- as.data.frame(dataTable)[, column_w_keggs]
  
  tempData <- as.data.frame(t(sapply(tempIDs, 
                     function(x) pull_respective_kegg_data(x))))
  
  names(tempData) <- c("kegg_id", "gene_name", "pathways", "kegg_orthologs")
  
  temp_combined_data <- dataTable %>% 
    left_join(mutate(tempData, kegg_id = as.character(kegg_id)), by = "kegg_id")
  
  return(temp_combined_data)
}

# Function to do the actual data pulling
pull_respective_kegg_data <- function(ID_of_int){
  
  tempKEGG <- try(keggGet(ID_of_int))
  
  tempEntry <- try(c(unname(tempKEGG[[1]]$ENTRY), tempKEGG[[1]]$NAME, 
                     paste(unname(tempKEGG[[1]]$PATHWAY), "_", 
                           collapse = "", sep = ""), 
                     paste(names(tempKEGG[[1]]$PATHWAY), "_", collapse = "", sep = "")))
  
  if(length(tempEntry) <= 1){
    
    tempEntry <- c(ID_of_int, NA, NA, NA)
  }
  
  print(paste("Completed processing KEGG ID: ", ID_of_int, sep = ""))
  
  return(tempEntry)
}


##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

# Create vectors to be used in aggregating the adenoma and crc models
models_used <- c("adn_model_imp", "crc_model_imp")

# Get the summary imp data
summary_imp_lists <- sapply(models_used, 
                            function(x) group_stats(x), simplify = F)

# Get top 25 most predictive genes according to picrust for adenoma
adn_imp_w_gene_data <- get_gene_ids(summary_imp_lists[["adn_model_imp"]], "kegg_id")
adn_imp_w_gene_data <- adn_imp_w_gene_data %>% slice(1:25)

# Get top 25 most predictive genes according to picrust for carcinoma
crc_imp_w_gene_data <- summary_imp_lists[["crc_model_imp"]] %>% 
  left_join(
    select(adn_imp_w_gene_data, kegg_id, gene_name, pathways, kegg_orthologs), by = "kegg_id")

crc_imp_w_gene_data <- crc_imp_w_gene_data %>% slice(1:25)














# CRC gene functions
# gingipains R - amino acid metabolism
# MFS Transporter - drug resistance
# taurine-pyruvate aminotransferase - bile acid
# 3-oxo-5-alpha-steroid 4-dehydrogenase 2  - steroid metabolism (human)
# aldehyde dehydrogenase - acetate (aerobic bacteria)
# PTS system (phosphotransferase system) - carbohydrate translocation and phosphorylation
# AraC family transcriptional regulator - Transcription 
# streptolysin S associated protein - toxin
# phosphotransacylase - acetate
# glycerol-1-phosphate dehydrogenase - membrane phospholipids (Archaea)
# cystine reductase - amino acid metabolism
# betaine reductase - amino acid metabolism
# 1-propanol dehydrogenase - alcohol 
# glucosyltransferase - glucose transfer
# galactitol-1-phosphate 5-dehydrogenase - galactose metabolism
# xanthine dehydrogenase - purine, pyrimidines, pterins, aldehydes
# LysR family transcriptional regulator - transcription

# Adn gene functions
# phosphotransacylase - acetate
# 1-propanol dehydrogenase - alcohol
# adenosylcobinamide hydrolase - cofactor biosynthesis
# fructose 1,6-bisphosphate aldolase/phosphatase - gluconeogenesis (thermophilic) (propionate)
# 3-dehydro-L-gulonate-6-phosphate decarboxylase - ascorbate and aldarate, pentose and glucuronate interconversion
# L-fuculokinase - fructose and mannose metabolism
# Bax protein - Cell death (human)
# gingipain R - amino acid metabolism
# crotonobetainyl-CoA:carnitine CoA-transferase  - CoA transferase (butyrate)
# DeoR family transcriptional regulator - Transcription
# tRNA-intron endonuclease - tRNA maturation
# anaerobic dimethyl sulfoxide reductase subunit C - anaerobic respiration
# trimethylamine corrinoid protein - methanogenesis
# cinnamoyl-CoA:phenyllactate CoA-transferase - CoA transferase
# phosphoribosylaminoimidazolecarboxamide formyltransferase - purine metabolism


