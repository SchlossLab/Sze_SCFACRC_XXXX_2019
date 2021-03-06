library(tidyverse)

options(warn = 1)

cv_test_compare_files <- function(file) {
	print(file)

	read_tsv(file, col_types=cols(condition=col_character(), .default=col_double())) %>%
		mutate(condition=str_replace(condition, "_", ".")) %>%
		separate(condition, into=c("class", "input"), sep="\\.")

}

pool_files <- function(input, output){


	input %>%
		map_dfr(., cv_test_compare_files) %>%
		write_tsv(., output)

}


input <- commandArgs(trailingOnly=TRUE) # recieve input from model
# Get variables from command line
output_file_name <- input[1]
input_file_names <- input[-1]

pool_files(input_file_names, output_file_name)
