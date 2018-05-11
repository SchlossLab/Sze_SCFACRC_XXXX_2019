### Selective search and test of specific OPFs
### Generating the relative abundances
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "optparse"))

# Command line input 
option_list = list(
  make_option(c("-mf", "--match_file"), type="character", default="data/process/select_scfa_opf_matches.tsv", 
              help="KEGG to contig match file name", metavar="character"),
  make_option(c("-of", "--opf_file"), type="character", default="data/process/opf_shared.tsv", 
              help="opf shared file [default= %default]", metavar="character"), 
  make_option(c("-md", "--metadata"), type="character", default="data/process/sra_meta_conversion.txt", 
              help="metadata of the samples [default= %default]", metavar="character"), 
  make_option(c("-o", "--output_file"), type="character", default="data/process/select_scfa_opf_kruskal_summary.csv", 
              help="Output kruskal summary name [default= %default]", metavar="character"))

opt_parser = OptionParser(option_list=option_list)
command_line_input = parse_args(opt_parser)

# Unzip needed files
system(paste("gunzip ", command_line_input[["opf_file"]], ".gz", sep = ""), intern = TRUE, ignore.stderr = TRUE)














