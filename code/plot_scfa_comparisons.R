library(tidyverse)
library(RColorBrewer)
library(cowplot)

### Helper function
capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}


### Read in scfa data and metadata
scfa_data <- read_tsv('data/scfa/scfa_composite.tsv') %>%
	mutate(scfa=factor(scfa, levels=c("acetate", "propionate", "isobutyrate", "butyrate"))) %>%
	filter(scfa != "pooled")

cross_section_metadata <- read_csv("data/metadata/cross_section.csv",
																col_type=cols(sample=col_character())) %>%
														mutate(dx=factor(dx, levels=c("normal", "adenoma", "cancer")))
follow_up_metadata <- read_csv("data/metadata/follow_up.csv",
														col_type=cols(initial=col_character(), followUp=col_character())) %>%
												mutate(dx=factor(dx, levels=c("normal", "adenoma", "cancer")))


### Plot SCFAs for each SCFA and diagnosis
cross_section <- inner_join(scfa_data, cross_section_metadata, by=c("study_id"="sample")) %>%
	select(scfa, mmol_kg, dx) %>%
	ggplot(aes(x=dx, y=mmol_kg, fill=dx)) +
		# geom_violin(show.legend=FALSE, draw_quantiles=0.5) +
		geom_boxplot(show.legend=FALSE) +
		facet_wrap(~scfa, scales="free_y", nrow=1, labeller=labeller(.default=capitalize)) +
		labs(x=NULL, y="mmol of SCFA/kg of feces") +
		scale_x_discrete(labels = NULL) +
		scale_fill_manual(
			breaks=c("normal", "adenoma", "cancer"),
			labels=c("Normal", "Adenoma", "Cancer"),
			values=brewer.pal(n=3, "YlOrRd"),
			name=NULL
		) +
		theme_classic() +
		theme(
			strip.background = element_blank(),
			axis.ticks.x = element_blank(),
			strip.text = element_text(size=10)
		)


### Plot SCFAs for each SCFA and diagnosis
longitudinal <- inner_join(scfa_data, follow_up_metadata, by=c("study_id"="initial")) %>%
	rename(initial = mmol_kg) %>%
	inner_join(scfa_data, ., by=c("study_id"="followUp", "scfa")) %>%
	rename(final = mmol_kg) %>%
	select(EDRN, dx, scfa, initial, final) %>%
	mutate(diff = final-initial) %>%
	filter(diff < 120) %>%
	group_by(scfa) %>%
	mutate(ymax = max(abs(diff)), ymin = -1 * ymax) %>%
	ungroup() %>%
	ggplot(aes(x=dx, y=diff, fill=dx)) +
		# geom_violin(draw_quantiles=0.5) +
		geom_boxplot() +
		geom_blank(aes(y = ymin)) + geom_blank(aes(y = ymax)) +
		facet_wrap(~scfa, scales="free_y", nrow=1, labeller=labeller(.default=capitalize)) +
		labs(x=NULL, y="Difference in SCFA concentration\nbefore and after treatment") +
		scale_x_discrete(labels = NULL) +
		scale_fill_manual(
			breaks=c("normal", "adenoma", "cancer"),
			labels=c("Normal", "Adenoma", "Cancer"),
			values=brewer.pal(n=3, "YlOrRd"),
			name=NULL,
			drop=FALSE
		) +
		theme_classic() +
		theme(
			strip.background = element_blank(),
			axis.ticks.x = element_blank(),
			axis.line.x = element_blank(),
			strip.text = element_text(size=10)
		)


### Merge the figures together
plot_grid(cross_section, longitudinal, nrow=2, labels=c("A", "B"))
ggsave("results/figures/scfa_comparisons.pdf", width=6.875, height=6, units="in")
