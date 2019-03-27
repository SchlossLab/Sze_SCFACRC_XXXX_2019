library(tidyverse)

read_optimum_mtry_file <- function(file) {

	read_csv(file, col_types=cols(.default=col_double())) %>%
		mutate(file=str_replace(file, "//", "/"),
					iteration=as.numeric(str_replace(file, ".*optimum.mtry.(\\d*).csv", "\\1")),
					condition=str_replace(file, "data/rf/(.*)/optimum.*", "\\1")
				) %>%
		select(condition, iteration, everything(), -file)

}

pool_files <- function(path){

	output_file <- paste0(path, "/cv_test_compare.tsv")

	list.files(path, "optimum_mtry", full.names=TRUE) %>%
		map_dfr(., read_optimum_mtry_file) %>%
		arrange(iteration) %>%
		write_tsv(., output_file)

}


input <- commandArgs(trailingOnly=TRUE) # recieve input from model
# Get variables from command line
path <- input[1]

pool_files(path)
