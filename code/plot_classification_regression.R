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


classification_plot <- model_data %>%
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
			labels = c("Without SCFA data", "With SCFA data"),
			values = brewer.pal("PRGn", n=3)[-2]
		) +
		labs(x="Microbiome data included", y="AUCROC for held out data") +
		theme_classic() +
		theme(
			strip.background = element_blank(),
			panel.spacing.x=unit(1, "lines"),
			axis.ticks.x = element_blank(),
			strip.text = element_text(size=10),
			legend.position=c(0.15, 0.8)
		)

#acetate, propionate, isobutyrate, butyrate
#Genus, OTU, KEGG, OPF, PC2Pathways, PC2KO
