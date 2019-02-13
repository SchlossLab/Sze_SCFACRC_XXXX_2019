# Create a OPF abundance table
# Need files clu.tsv, and orf_abund.tsv
# Within workflow pathing for the files are
# 	data/metagenome/mmseq2/clu.tsv
# 	data/metagenome/diamond/orf_abund.tsv
#
#	Previously, MS outputted a table with a column that was
#		cor_total_counts = sum_counts/total_opf_members
#
# This doesn't make sense since the output from orf_abund.tsv only allows each read to map to one
# orf sequence. I'm going to go with the total number of reads mapped to each ORF as the count

library(tidyverse)

# Command line input
input <- commandArgs(TRUE)
opf_filename <- input[1]		#opf_filename <- "data/metagenome/mmseq2/clu.tsv"
orf_filename <- input[2]		#orf_filename <- "data/metagenome/diamond/orf_abund.tsv"

# Load in needed data
opf_clusters <- read_tsv(opf_filename, col_names=c("opf_cluster", "seq_name"))
orf_count_data <- read_tsv(orf_filename, col_names=c("count", "seq_name", "sample_id"))

# Create a tidy OPF abundance table
combined_data <- inner_join(orf_count_data, opf_clusters, by = "seq_name") %>%
  group_by(sample_id, opf_cluster) %>%
  summarise(total_counts = sum(count)) %>%
	ungroup()

#	combined_data %>% group_by(sample_id) %>% summarize(N=sum(total_counts)) %>% arrange(N)
#	combined_data %>% group_by(sample_id) %>% summarize(N=sum(total_counts)) %>% arrange(desc(N))
#	These commands show ~10-fold difference in the number of OPFs across the 85 samples. Will
# remove samples with fewer than 300,000 hits to OPFs and will subsample to that many hits/reads

subsample <- function(n, df){
	unroll <- rep(df$opf_cluster, df$total_counts)
	sample(unroll, n) %>% table() %>% as_tibble(.name_repair=function(x){c("opf", "count")})
}

min_hits <- 300000

combined_data %>%
	group_by(sample_id) %>%
	summarize(hits = sum(total_counts)) %>%
	filter(hits > min_hits) %>%
	select(sample_id) %>%
	inner_join(., combined_data, by="sample_id") %>%
	group_by(sample_id) %>%
	nest() %>%
	mutate(sub = map(data, ~subsample(min_hits, .x))) %>%
	unnest(sub) %>%
	write_tsv("data/metagenome/opf.tidy_shared.tsv")
