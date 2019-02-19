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


# Need to match sample ids between the SRA and EDRN
hannigan_metadata_filename <- "data/metadata/hannigan_metadata.tsv"
sra_metadata_filename <- "data/raw/metagenome/SRP108915_info.csv"

#	generate table containing Hannigan's MG file names and Zackular sample ID
hannigan <- read_tsv(hannigan_metadata_filename, col_type=cols(LibraryDate=col_character())) %>%
	select(MetaID, SampleID) %>%
	mutate(SampleID = str_replace(SampleID, "^(.)0(.)$", "\\1\\2"))

# generate table containing the SRA sample ID and the Hannigan MG file name
sra <- read_csv(sra_metadata_filename) %>% select(Run, SampleName)

# merge the three tables so we can convert between the SRA sample ID and the EDRN
metadata_lookup <- inner_join(sra, hannigan, by=c("SampleName"="MetaID")) %>%
										rename(sample_id=Run) %>%
										select(-SampleName)




#Need to read in various data files for merging and generating output
opf_filename <- "data/metagenome/mmseq2/clu.tsv"
orf_filename <- "data/metagenome/diamond/orf_abund.tsv"
kegg_filename <- "data/metagenome/mmseq2/resultDB.m8"


# Load in needed data. The three following files that I read in need column names added.
opf_clusters <- fread(opf_filename, header=FALSE) %>%
		as_tibble(.name_repair=function(x){c("opf_cluster", "seq_name")})

orf_count_data <- fread(orf_filename, header=FALSE) %>%
		as_tibble(.name_repair=function(x){c("count", "seq_name", "sample_id")}) %>%
		inner_join(., metadata_lookup, by="sample_id") %>%
		select(-sample_id) %>%
		rename(sample_id=SampleID)

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
	write_tsv("data/metagenome/metag.ko_lookup.tsv")


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
	select(sample_id, opf, k, count)


#output opf shared file
opf_shared <- subsampled_data %>%
	select(-k) %>%
	spread(opf, count, fill=0) %>%
	mutate(label="opf", numOtus=ncol(.)-1) %>%
	select(label, sample_id, numOtus, everything()) %>%
	rename(Group = sample_id) %>%
	write_tsv("data/metagenome/metag.opf.shared")


#output kegg shared file
kegg_shared <- subsampled_data %>%
	select(-opf) %>%
	group_by(sample_id, k) %>%
	summarize(count = sum(count)) %>%
	ungroup() %>%
	spread(k, count, fill=0) %>%
	mutate(label="kegg", numOtus=ncol(.)-1) %>%
	select(label, sample_id, numOtus, everything()) %>%
	rename(Group = sample_id) %>%
	write_tsv("data/metagenome/metag.kegg.shared")
