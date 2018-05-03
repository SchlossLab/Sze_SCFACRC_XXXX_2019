### Getting the ID of the key OPF clusters
### Meant to be run after the diamond_kegg_search.py program
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse"))

library(KEGGREST)

# Diamond output file has the following default columns (left to right) --- from uniprot search
# qseqid = Query Seq - id, sseqid = Subject Seq - id, 
# pident = Percentage of identical matches,  length = Alignment length, 
# mismatch = Number of mismatches, gapopen = Number of gap openings, 
# qstart = Start of alignment in query, qend = End of alignment in query, 
# sstart = Start of alignment in subject, send = End of alignment in subject, 
# evalue = Expect value, bitscore = Bit score

# Load in needed data
comparison_data <- read_csv("data/process/opf_dunn_select_results.csv")

# Create function to pull bacterial species  
get_bacterium <- function(kegg_species_id){
  
  test <- keggInfo(kegg_species_id)
  
  tempdf <- as.data.frame(test) %>% 
    rename(info = test) %>% 
    separate(info, c("X1", "X2", "X3"), sep = "           ") %>% 
    separate(X2, c("species", "extra"), sep = " KEGG") %>% 
    pull(species)
  
  return(tempdf)
}

# Load in Kegg data
kegg_data <- read_tsv("data/process/sig_dunn_kegg_protein_matches.tsv", col_names = F) %>% 
  rename(query_seq = X1, subject_seq = X2, pident = X3, align_length = X4, mismatch = X5, gapopen = X6, 
         query_start = X7, query_end = X8, subject_start = X9, subject_end = X10, evalue = X11, bitscore = X12) %>% 
  separate(subject_seq, c("X1", "X2", "X3", "X4"), sep = "\\|") %>% 
  mutate(X4 = case_when(str_detect(X3, "K0") == TRUE ~ X3, 
                        str_detect(X3, "none") == TRUE ~ X3, 
                        TRUE ~ X4), 
         X2 = str_replace_all(X2, "_", " ")) %>% 
  rename(species = X1, gene = X2, pathway_or_complex = X3, kegg_ortholog = X4) %>% 
  separate(species, c("species_id", "extra_info"), sep = ":") %>% 
  mutate(specific_species = map(species_id, function(x) get_bacterium(x))) %>% 
  unnest(specific_species) %>% 
  select(query_seq, species_id, specific_species, gene, pathway_or_complex, kegg_ortholog, everything())





