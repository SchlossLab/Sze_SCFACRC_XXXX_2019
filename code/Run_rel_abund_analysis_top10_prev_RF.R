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
                            "Roseburia", "Escherichia/Shigella", "Clostridium_XIVa"), 
                    adv_adn = c("Clostridiales", "Clostridium_XIVa", "Roseburia", "Odoribacter", "Bacteroides", 
                                "Ruminococcaceae", "Streptococcus", "Blautia", "Lachnospiraceae", "Ruminococcus"), 
                    crc = c("Porphyromonas", "Parvimonas", "Fusobacterium", "Gemella", "Prevotella", "Streptococcus", 
                            "Coprobacillus", "Pasteurellaceae", "Collinsella", "Bilophila"))

shared_data <- read_tsv("data/process/final.0.03.subsample.shared")


otus_present <- colnames(select(shared_data, -label, -Group, -numOtus))


taxonomic_data <- read_tsv("data/process/final.taxonomy") %>% 
  mutate(Taxonomy = str_replace_all(Taxonomy, "\\((\\d{2,3})\\)", ""), 
         Taxonomy = str_replace_all(Taxonomy, "_unclassified", ""), 
         Taxonomy = str_replace_all(Taxonomy, "2", "")) %>% 
  separate(Taxonomy, c("kingdom", "phyla", "class", "order", "family", "genus", "species"), sep = ";") %>% 
  filter(OTU %in% otus_present)

combined_meta_data <- read_csv("data/raw/metadata/metaI_final.csv") %>% 
  mutate(Dx_Bin = ifelse(Dx_Bin == "High Risk Normal" | Dx_Bin == "Normal", invisible("control"), 
                         ifelse(Dx_Bin == "Adenoma", invisible("adenoma"), 
                                ifelse(Dx_Bin == "adv Adenoma", invisible("adv_adenoma"), 
                                       ifelse(Dx_Bin == "Cancer", invisible("cancer"), invisible(Dx_Bin)))))) %>% 
  select(sample, Dx_Bin, dx) %>% 
  bind_rows(read_csv("data/raw/metadata/good_metaf_final.csv") %>% 
              mutate(sample = initial) %>% 
              select(sample, Dx_Bin, dx))
  


##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################

# Function to generate tax lists for select IDs by tumor type
get_relevent_taxa <- function(tumor_type, taxa_list, tax_file){
  
  temp_df <- tax_file %>% filter(genus %in% taxa_list[[tumor_type]])
  
  return(temp_df)
}


# Function to prune shared file based on taxIDs and by meta_data then remove nzv data
get_specific_otu <- function(tumor_type, metafile, sharedfile, taxaList){
  
  tumor_names <- list(adn = "adenoma", adv_adn = "adv_adenoma", crc = "cancer")
  temp_taxa_file <- as.data.frame(taxaList[[tumor_type]])
  
  temp_df <- metafile %>% select(sample, Dx_Bin) %>% 
    left_join(sharedfile, by = c("sample" = "Group")) %>% 
    filter(Dx_Bin == "control" | Dx_Bin == tumor_names[[tumor_type]]) %>% 
    select(sample, Dx_Bin, one_of(temp_taxa_file$OTU))
  
  nzv <- nearZeroVar(temp_df)
  
  temp_df <- temp_df[, -nzv]
  
  return(temp_df)
}


# Function to run the wilcoxson rank sum test on each provided OTU
get_wilcox_testing <- function(tumor_type, test_data, taxaList){
  
  tempData <- as.data.frame(test_data[[tumor_type]])
  tempTaxa <- taxaList[[tumor_type]] %>% 
    select(OTU, genus)
  
  tumor_names <- list(adn = "adenoma", adv_adn = "adv_adenoma", crc = "cancer")
  
  vars_to_test <- colnames(select(test_data[[tumor_type]], -sample, -Dx_Bin))
  
  test <- t(sapply(vars_to_test, 
         function(x) wilcox.test(
           filter(tempData, Dx_Bin == "control")[, x], 
           filter(tempData, Dx_Bin == tumor_names[[tumor_type]])[, x], 
           alternative = "two.sided")$p.value, simplify = F) %>% bind_rows()) %>% 
    as.data.frame() %>% 
    mutate(otus = rownames(.), 
           bh = p.adjust(V1, method = "BH")) %>% 
    left_join(tempTaxa, by = c("otus" = "OTU")) %>% 
    rename(pvalue = V1) %>% 
    select(otus, genus, pvalue, bh) %>% 
    arrange(pvalue)

  
  
  return(test)
}



##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

# generate vector to make sapply run
tumors <- c("adn", "adv_adn", "crc")
# Grabs the specific taxa in the top 10 of the previous models
specific_taxa <- sapply(tumors, 
                        function(x) get_relevent_taxa(x, lowest_txID, taxonomic_data), simplify = F)
# Pares down the candidate taxa OTUs to a manageable number
taxa_for_testing <- sapply(tumors, 
                           function(x) get_specific_otu(x, combined_meta_data, shared_data, specific_taxa), simplify = F)
# Runs a standard wilcoxson rank sum test and a BH correction 
test_result_tables <- sapply(tumors, 
                             function(x) get_wilcox_testing(x, taxa_for_testing, specific_taxa), simplify = F)

# Write out the data
sapply(tumors, 
       function(x) write_csv(test_result_tables[[x]], 
                             paste("data/process/tables/", x, "_16S_top10_RF_taxa_testing.csv", sep = "")))





