library(tidyverse)
library(RColorBrewer)

### Helper function
capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}




regression_file_name <- "data/rf/regression_data_pool.tsv"

regression_features <- c("genus", "otu", "kegg", "opf", "pc2pathways", "pc2ko")
scfa_types <- c("acetate", "propionate", "butyrate")

regression_data <- read_tsv(regression_file_name,
					col_types=cols(class=col_character(),
													input=col_character(),
													.default=col_double()
												)
				) %>%
				mutate(microbiome=str_replace(input, "fit_", ""),
							microbiome=str_replace(microbiome, "_fit", "")
		) %>%
	select(class, microbiome, train_Rsquared, test_Rsquared) %>%
	filter(class %in% scfa_types) %>%
	filter(microbiome %in% regression_features) %>%
	mutate(class = factor(class, levels=scfa_types),
					microbiome = factor(microbiome, levels=regression_features),
					data_type = case_when(
			       microbiome %in% c("genus", "otu")  ~ "16S rRNA gene",
						 microbiome %in% c("kegg", "opf")  ~ "Metagenomics",
						 microbiome %in% c("pc2ko", "pc2pathways")  ~ "PICRUSt"
					 )
				) %>%
	gather(test, auc, train_Rsquared, test_Rsquared) %>%
	mutate(test=factor(test, levels=c("train_Rsquared", "test_Rsquared")))

regression_data %>%
	ggplot(aes(x=microbiome, y=auc, fill=test)) +
	geom_boxplot() +
	facet_grid(class~data_type, labeller=labeller(.default=capitalize), scale="free_x") +
	scale_x_discrete(
		breaks = c("genus", "otu", "kegg", "opf", "pc2ko", "pc2pathways"),
		labels = c("Genus", "OTU", "KEGG", "OPF", "KEGG", "Paths")
	) +
	scale_fill_manual(
		name = NULL,
		breaks = c("train_Rsquared", "test_Rsquared"),
		labels = c("Cross validation", "Held out data"),
		values = brewer.pal("RdBu", n=3)[-2]
	) +
	labs(x="Microbiome data included", y=bquote(R^{2}~'for regression of SCFA concentrations')) +
	theme_classic() +
	theme(
		strip.background = element_blank(),
		panel.spacing.x=unit(1, "lines"),
		axis.ticks.x = element_blank(),
		strip.text = element_text(size=10)
	) + ggsave("results/figures/regression_testing.pdf", width=6.875, height=5)
