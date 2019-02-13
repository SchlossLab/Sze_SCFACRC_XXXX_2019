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
library(data.table)

opf_filename <- "data/metagenome/mmseq2/clu.tsv"
orf_filename <- "data/metagenome/diamond/orf_abund.tsv"
kegg_filename <- "data/metagenome/mmseq2/resultDB.m8"

# Load in needed data
opf_clusters <- fread(opf_filename, header=FALSE) %>%
		as_tibble(.name_repair=function(x){c("opf_cluster", "seq_name")})

orf_count_data <- fread(orf_filename, header=FALSE) %>%
		as_tibble(.name_repair=function(x){c("count", "seq_name", "sample_id")})

kegg_data <- fread(kegg_filename, header=FALSE) %>%
		as_tibble(.name_repair=function(x){c("seq_name", "target", "identity", "aligment_length",
		 																		"n_mismatches", "n_gap_openings", "q_start", "q_end",
																				"t_start", "t_end", "e_value", "bit_score")
																			}
							) %>%
		select(seq_name, target) %>%
		mutate(k=str_replace(target, ".*\\|", "")) %>%
		group_by(seq_name) %>%
		count(k, sort=TRUE) %>%
		summarize(k = k[1]) %>%
		ungroup() %>%
		select(seq_name, k)

#	need to do a left join between kegg_data and the other data to map annotations to OPFs and get
#	consensus annotation for each OPF
get_consensus <- function(df){
	tibble(k=rep(df$k, df$count)) %>% count(k, sort=TRUE) %>% first()
}

annotated_reads <- left_join(opf_clusters, kegg_data, by="seq_name") %>%
	mutate(k=ifelse(is.na(k), "none", k)) %>%
	inner_join(orf_count_data, ., by = "seq_name") %>%
	group_by(opf_cluster) %>%
	nest() %>%
	mutate(annotation = map(data, ~get_consensus(.x))) %>%
	unnest(annotation) %>%
	select(opf_cluster, k) %>%
	write_tsv("data/metagenome/opf_k.annotations.tsv")


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

subsampled_data <- combined_data %>%
	group_by(sample_id) %>%
	summarize(hits = sum(total_counts)) %>%
	filter(hits > min_hits) %>%
	select(sample_id) %>%
	inner_join(., combined_data, by="sample_id") %>%
	group_by(sample_id) %>%
	nest() %>%
	mutate(sub = map(data, ~subsample(min_hits, .x))) %>%
	unnest(sub) %>%
	inner_join(., annotated_reads, by=c("opf"="opf_cluster")) %>%
	select(sample_id, opf, k, count) %>%
	write_tsv("data/metagenome/opf_k.tidy_shared.tsv")
