library(tidyverse)
library(data.table)

scfa_keggs <- c(
					"K00929", #butyrate kinase [EC:2.7.2.7]
					"K01034", #acetate CoA/acetoacetate CoA-transferase alpha subunit [EC:2.8.3.8 2.8.3.9]
					"K01035", #acetate CoA/acetoacetate CoA-transferase beta subunit [EC:2.8.3.8 2.8.3.9]
					"K01908", #propionyl-CoA synthetase [EC:6.2.1.17]
					"K01895", #acetyl-CoA synthetase [EC:6.2.1.1]
					"K00932", #propionate kinase [EC:2.7.2.15]
					"K19697", #propionate kinase [EC:2.7.2.15]
					"K18372", #methyl acetate hydrolase [EC:3.1.1.-]
					"K00467", #lactate 2-monooxygenase [EC:1.13.12.4]
					"K00156", #pyruvate dehydrogenase (quinone) [EC:1.2.5.1]
					"K00925", #acetate kinase [EC:2.7.2.1]
					"K01512", #acylphosphatase [EC:3.6.1.7]
					"K01067", #acetyl-CoA hydrolase [EC:3.1.2.1]
					"K18118", #succinyl-CoA:acetate CoA-transferase [EC:2.8.3.18]
					"K01026", #propionate CoA-transferase [EC:2.8.3.1]
					"K01905", #acetate---CoA ligase (ADP-forming) subunit alpha [EC:6.2.1.13]
					"K22224", #acetate---CoA ligase (ADP-forming) subunit beta [EC:6.2.1.13]
					"K01895", #acetyl-CoA synthetase [EC:6.2.1.1]
					"K00128", #aldehyde dehydrogenase (NAD+) [EC:1.2.1.3]
					"K14085", #aldehyde dehydrogenase family 7 member A1 [EC:1.2.1.31 1.2.1.8 1.2.1.3]
					"K00149", #aldehyde dehydrogenase family 9 member A1 [EC:1.2.1.47 1.2.1.3]
					"K00138"	#aldehyde dehydrogenase [EC:1.2.1.-]
				) %>% sort()

get_scfa_keggs <- function(shared_file_name, lookup_file_name=NULL){

	full_shared <- fread(shared_file_name, header=T) %>% as_tibble()

	if(!is.null(lookup_file_name)){
		scfa_opfs <- read_tsv(lookup_file_name) %>% filter(k %in% scfa_keggs) %>% pull(opf_cluster)
		overlapping <- scfa_opfs[scfa_opfs %in% colnames(full_shared)]
	} else {
		overlapping <- scfa_keggs[scfa_keggs %in% colnames(full_shared)]
	}

	output_shared <- str_replace(shared_file_name, ".shared", "scfa.shared")

	full_shared %>%
		select(label, Group, numOtus, overlapping) %>%
		mutate(numOtus=ncol(.)-3) %>%
		write_tsv(output_shared)

}
