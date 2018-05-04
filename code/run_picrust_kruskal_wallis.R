### Investigate how well the picrust predicted metagenomes work
### compare most differentially expressed genes by disease (kruskal - wallis)
### Marc Sze



# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "biomformat", "caret"))


# Load in needed metadata and rename sample column
picrust_meta <- read_tsv("data/process/picrust_metadata") %>% 
  mutate(Group = as.character(Group))

# Load in more specific follow up data
metaF <- read_csv("data/raw/metadata/good_metaf_final.csv")

# Load in predicted metagenome data
all_biom_data <- read_biom("data/process/predicted_metagenomes.biom")

# Converts the needed biom data into a matrix
full_pi_data <- t(as(biom_data(all_biom_data), "matrix")) %>% 
  as.data.frame() %>% 
  mutate(sample_id = rownames(.)) %>% 
  select(sample_id, everything())

# Get the kegg mappings for each coded sample
kegg_table <- observation_metadata(all_biom_data) %>% 
  mutate(kegg_id = rownames(.)) %>% 
  select(kegg_id, everything())

# Generate the near zero variance
nzv <- nearZeroVar(full_pi_data)

# Remove these from the data set
nzv_full_pi_data <- full_pi_data[, -nzv]


##############################################################################################
############### List of functions to get things to run nice ##################################
##############################################################################################





















##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################






















