
library("tidyverse")


clean_tax <- function(){

# Correct the taxonomy file for cases where greengenes doesn't know what the taxonomy is
## Read in the conserved taxonomy file
	read_tsv("data/mothur/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons_gg.taxonomy") %>%
		mutate(Taxonomy=str_replace_all(Taxonomy, ".*unknown.*", "k__Bacteria(100);p__(100);c__(100);o__(100);f__(100);g__(100);s__(100);")) %>%
		write_tsv("data/picrust1/crc.cons.taxonomy")

}



get_shared_annotation <- function(){

	pc_data <- read_tsv("data/picrust1/crc.metagenomes.tsv", skip=1) %>%
		rename(kegg = `#OTU ID`,
					description = KEGG_Description)

	keggs <- pc_data %>% select(kegg, description)

	freqs <- pc_data %>%
		select(-description) %>%
		gather(Group, count, -kegg)

	good_otus <- freqs %>%
		group_by(kegg) %>%
		summarize(n=sum(count)) %>%
		filter(n!=0) %>%
		select(kegg)

	inner_join(freqs, good_otus, by="kegg") %>%
		spread(kegg, count) %>%
		mutate(label="picrust1", numOtus=ncol(.)-1) %>%
		select(label, Group, numOtus, everything()) %>%
		write_tsv('data/picrust1/crc.picrust1.shared')

	inner_join(keggs, good_otus, by="kegg") %>%
		write_tsv('data/picrust1/crc.picrust1.keggs')

}
