library(tidyverse)
library(RColorBrewer)
library(cowplot)

### Helper function
capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}


data_file_name <- "data/rf/classification_data_pool.tsv"

model_data <- read_tsv(data_file_name) %>%
	mutate(scfa = str_detect(input, "scfa"),
				microbiome=str_replace(input, "scfa_", ""),
				microbiome=str_replace(microbiome, "_scfa", "")
		) %>%
	select(class, microbiome, scfa, cv_aucs, test_aucs) %>%
	filter(class != "lesion")


model_data %>%
	filter(microbiome %in% c("scfa", "otu", "genus")) %>%
	mutate(microbiome = factor(microbiome, levels=c("scfa", "genus", "otu"))) %>%
	ggplot(aes(x=microbiome, fill=scfa, y=test_aucs)) +
	 geom_boxplot(position=position_dodge2(preserve = "single"))+
	#	geom_violin(draw_quantiles=0.5)+#, position=position_dodge(preserve = "single")) +
		facet_wrap(~class) +
		coord_cartesian(ylim=c(0.35,1.0)) +
		ggsave("test.pdf")
