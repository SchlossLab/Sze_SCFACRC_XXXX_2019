library(tidyverse)
library(RColorBrewer)

### Helper function
capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}


classification_file_name <- "data/rf/classification_data_pool.tsv"

classification_data <- read_tsv(classification_file_name,
												col_types=cols(class=col_character(),
																				input=col_character(),
																				.default=col_double()
																			)
											) %>%
	mutate(scfa = str_detect(input, "scfa"),
				microbiome=str_replace(input, "scfa_", ""),
				microbiome=str_replace(microbiome, "_scfa", "")
		) %>%
	select(class, microbiome, scfa, cv_aucs, test_aucs) %>%
	filter(class != "lesion") %>%
	filter(microbiome %in% c("scfa", "otu", "genus")) %>%
	mutate(microbiome = factor(microbiome, levels=c("scfa", "genus", "otu"))) %>%
	mutate(scfa = ifelse(scfa, "With SCFAs", "Without SCFAs")) %>%
	gather(test, auc, cv_aucs, test_aucs)





classification_data %>%
	ggplot(aes(x=microbiome, fill=test, y=auc)) +
		geom_hline(yintercept=0.5, col="#888888") +
		geom_boxplot(position=position_dodge2(preserve = "single"))+
		facet_grid(scfa~class, labeller=labeller(.default=capitalize)) +
		coord_cartesian(ylim=c(0.35,1.0)) +
		scale_x_discrete(
			breaks = c("scfa", "genus", "otu"),
			labels = c("None", "Genus", "OTU")
		) +
		scale_fill_manual(
			name = NULL,
			breaks = c("cv_aucs", "test_aucs"),
			labels = c("Cross validation", "Held out data"),
			values = brewer.pal("RdBu", n=3)[-2]
		) +
		labs(x="Microbiome data included", y="AUCROC for diagnosis classification") +
		theme_classic() +
		theme(
			strip.background = element_blank(),
			panel.spacing.x=unit(1, "lines"),
			axis.ticks.x = element_blank(),
			strip.text = element_text(size=10)
		) + ggsave("results/figures/classification_testing.pdf", width=6.875, height=5)
