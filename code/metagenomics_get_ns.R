library(tidyverse)
library(data.table)

fread("data/metagenome/metag.kegg.shared")$Group %>%
	tibble(group=.) %>%
	mutate(dx = str_replace_all(group, "\\d", ""),
				dx = case_when(
					dx == "C" ~ "cancer",
					dx == "A" ~ "adenoma",
					dx == "H" ~ "normal"
				)
		) %>%
	count(dx) %>%
	write_tsv("data/metagenome/metag.sample.counts")
