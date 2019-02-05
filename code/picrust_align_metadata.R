### Align Metadata with shared file
### Preps a metadata file for creation of a biom file for picrust
### Marc Sze; modified by Pat Schloss


# Load needed libraries
library("tidyverse")

## Read in the conserved taxonomy file
read_tsv("data/mothur/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons_gg.taxonomy") %>%
	mutate(Taxonomy=str_replace_all(Taxonomy, ".*unknown.*", "k__Bacteria(100);p__(100);c__(100);o__(100);f__(100);g__(100);s__(100);")) %>%
	write_tsv("data/picrust/crc.cons.taxonomy")

## Read in needed data tables
metaInitial <- read_csv("data/raw/metadata/cross_section.csv",
												col_types=list(sample=col_character())) %>%
  select(sample, dx, fit_result) %>%
	mutate(time_point = "initial") %>%
  rename(Group = sample)

metaFollowUp <- read_csv("data/raw/metadata/follow_up.csv",
												col_types=list(initial=col_character(), followUp=col_character())) %>%
  gather(key = "time_point", value = "Group", initial, followUp) %>%
  select(Group, dx, fit_result, time_point)


# Merge the metadata together into a single table
metadata <- full_join(metaFollowUp, metaInitial)


# Get the shared file and join the metadata
meta_shared <- read_tsv("data/mothur/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.shared",
									col_types=list(Group=col_character())) %>%
					inner_join(metadata, shared, by="Group")




# reorder combined meta file to match the shared and write to data/picrust
meta_shared %>%
	select(Group, dx, fit_result, time_point) %>%
	write_tsv("data/picrust/crc.metadata")

meta_shared %>%
	select(-dx, -fit_result, -time_point) %>%
	write_tsv("data/picrust/crc.shared")
