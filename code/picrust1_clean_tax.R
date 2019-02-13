# Correct the taxonomy file for cases where greengenes doesn't know what the taxonomy is

library("tidyverse")

## Read in the conserved taxonomy file
read_tsv("data/mothur/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons_gg.taxonomy") %>%
	mutate(Taxonomy=str_replace_all(Taxonomy, ".*unknown.*", "k__Bacteria(100);p__(100);c__(100);o__(100);f__(100);g__(100);s__(100);")) %>%
	write_tsv("data/picrust1/crc.cons.taxonomy")
