library(tidyverse)
library(RColorBrewer)
library(cowplot)

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
	filter(class != "lesion")


classification_plot <- classification_data %>%
	filter(microbiome %in% c("scfa", "otu", "genus")) %>%
	mutate(microbiome = factor(microbiome, levels=c("scfa", "genus", "otu"))) %>%
	ggplot(aes(x=microbiome, fill=scfa, y=test_aucs)) +
		geom_hline(yintercept=0.5, col="#888888") +
		geom_boxplot(position=position_dodge2(preserve = "single"))+
		facet_wrap(~class, labeller=labeller(.default=capitalize)) +
		coord_cartesian(ylim=c(0.35,1.0)) +
		scale_x_discrete(
			breaks = c("scfa", "genus", "otu"),
			labels = c("None", "Genus", "OTU")
		) +
		scale_fill_manual(
			name = NULL,
			breaks = c(FALSE, TRUE),
			labels = c("Without SCFAs", "With SCFAs"),
			values = brewer.pal("PRGn", n=3)[-2]
		) +
		labs(x="Microbiome data included", y="AUCROC for\ndiagnosis classification") +
		theme_classic() +
		theme(
			strip.background = element_blank(),
			panel.spacing.x=unit(1, "lines"),
			axis.ticks.x = element_blank(),
			strip.text = element_text(size=10)
		)


regression_file_name <- "data/rf/regression_data_pool.tsv"

regression_features <- c("genus", "otu", "kegg", "opf", "pc2pathways", "pc2ec")
scfa_types <- c("acetate", "propionate", "isobutyrate", "butyrate")

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
			       microbiome %in% c("genus", "otu")  ~ "16S",
						 microbiome %in% c("kegg", "opf")  ~ "metagenome",
						 microbiome %in% c("pc2ec", "pc2pathways")  ~ "picrust"
					 )
				)

regression_features_plot <- regression_data %>%
	ggplot(aes(x=microbiome, y=test_Rsquared, fill=data_type)) +
	geom_boxplot() +
	facet_wrap(~class, nrow=2, labeller=labeller(.default=capitalize)) +
	scale_x_discrete(
		breaks = c("genus", "otu", "kegg", "opf", "pc2ec", "pc2pathways"),
		labels = c("Genus", "OTU", "KEGG", "OPF", "KEGG", "Paths")
	) +
	scale_fill_manual(
		name = NULL,
		breaks = c("16S", "metagenome", "picrust"),
		labels = c("16S rRNA gene", "Metagenomic", "PICRUSt"),
		values = brewer.pal("RdYlBu", n=3)
	) +
	labs(x="Microbiome data included", y=bquote(R^{2}~'for regression of SCFA concentrations')) +
	theme_classic() +
	theme(
		strip.background = element_blank(),
		panel.spacing.x=unit(1, "lines"),
		axis.ticks.x = element_blank(),
		strip.text = element_text(size=10)
	)

plot_grid(classification_plot, regression_features_plot, nrow=2, labels=c("A", "B"), rel_heights=c(0.5, 1))

ggsave("results/figures/scfa_modeling.pdf", width=6.875, height=6, units="in")
