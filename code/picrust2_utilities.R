library(tidyverse)
library(data.table)

sort_and_unalign_sequences <- function(){

	fasta <- scan("data/asv/crc.asv.fasta", quiet=T, what=character(), sep="\n")

	names <- fasta[c(T,F)] %>% str_replace(">(\\S*)\\s*", "\\1")
	sequences <- fasta[c(F,T)] %>% str_replace_all("[-.]", "")

	tibble(names, sequences) %>%
		arrange(names) %>%
		mutate(fasta=paste0(">", names, "\n", sequences)) %>%
		pull(fasta) %>%
		write("data/picrust2/crc.asv.fasta")

}


convert_tsv_to_shared <- function(tsv_file_name){

	label_tag <- str_replace(tsv_file_name, ".*/(.*)_out.*", "\\1") %>% str_to_lower() %>% str_replace("_metagenome", "")

	fread(tsv_file_name, header=TRUE) %>%
		as_tibble() %>%
		rename(annotation=1) %>%
		gather(Group, count, -annotation) %>%
		spread(annotation, count) %>%
		mutate(label=label_tag, numOtus=ncol(.)-1) %>%
		select(label, Group, numOtus, everything()) %>%
		write_tsv(paste0('data/picrust2/crc.', label_tag, ".shared"))

}
