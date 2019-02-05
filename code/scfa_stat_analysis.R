library(tidyverse)
library(broom)

scfa <- read_tsv("data/scfa/scfa_composite.tsv") %>%
				group_by(study_id, scfa) %>% # there are duplicate readings, average to get a single read
				summarize(mmol_kg = mean(mmol_kg)) %>%
				ungroup()


#cross section analysis
get_mean_dx <- function(x){
														x %>%
															group_by(dx) %>%
															summarize(mean=mean(mmol_kg)) %>%
															spread(key=dx, value=mean)
													}

read_csv("data/raw/metadata/cross_section.csv", col_type=cols(sample=col_character())) %>%
	inner_join(scfa, ., by=c("study_id"="sample")) %>%
	select(scfa, dx, mmol_kg) %>%
	group_by(scfa) %>%
	nest() %>%
	mutate(model=map(data, ~kruskal.test(x=.x$mmol_kg, g=as.factor(.x$dx)) %>% tidy())) %>%
	mutate(mean = map(data, get_mean_dx)) %>%
	unnest(model, mean) %>%
	mutate(model=map(data,
								~pairwise.wilcox.test(x=.x$mmol_kg, g=as.factor(.x$dx),  p.adjust.method="BH") %>%
									tidy() %>%
									mutate(compare=paste(group1, group2, sep="-")) %>%
									select(-group1, -group2) %>%
									spread(key=compare, value=p.value)
								)
				) %>%
	unnest(model) %>%
	select(-data, -parameter, -statistic) %>%
	write_tsv("data/process/scfa_cross_section_stats.tsv")


#pre-post analysis
read_csv("data/raw/metadata/follow_up.csv",
				col_type=cols(initial=col_character(), followUp=col_character(), EDRN=col_character())) %>%
	gather(sample, study_id, initial, followUp) %>%
	select(sample, study_id, EDRN, dx) %>%
	inner_join(scfa, ., by=c("study_id"="study_id")) %>%
	select(-study_id) %>%
	spread(key="sample", value="mmol_kg") %>%
	drop_na() %>%
	group_by(scfa, dx) %>%
	nest() %>%
	mutate(model=map(data, ~wilcox.test(x=.x$initial, y=.x$followUp, paired=TRUE, conf.int=TRUE) %>% tidy())) %>%
	unnest(model) %>%
	select(scfa, dx, estimate, p.value, conf.low, conf.high, method) %>%
	write_tsv("data/process/scfa_pre_post_stats.tsv")
