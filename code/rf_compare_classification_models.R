library(tidyverse)
library(broom)

data_file_name <- "data/rf/classification_data_pool.tsv"
model_data <- read_tsv(data_file_name)

# compare cv and test models' AUCs
model_data %>%
	group_by(class, input) %>%
	nest() %>%
	mutate(test = map(data, ~tidy(wilcox.test(x=.x$cv_aucs, y=.x$test_aucs, exact=FALSE))),
				cv_lci = map(data, ~quantile(.x$cv_aucs, 0.025)),
				cv_q1 = map(data, ~quantile(.x$cv_aucs, 0.25)),
				cv_median = map(data, ~median(.x$cv_aucs)),
				cv_q3 = map(data, ~quantile(.x$cv_aucs, 0.75)),
				cv_uci = map(data, ~quantile(.x$cv_aucs, 0.975)),
				test_lci = map(data, ~quantile(.x$test_aucs, 0.025)),
				test_q1 = map(data, ~quantile(.x$test_aucs, 0.25)),
				test_median = map(data, ~median(.x$test_aucs)),
				test_q3 = map(data, ~quantile(.x$test_aucs, 0.75)),
				test_uci = map(data, ~quantile(.x$test_aucs, 0.975))
			) %>%
	select(class, input, test,
				cv_lci, cv_q1, cv_median, cv_q3, cv_uci,
				test_lci, test_q1, test_median, test_q3, test_uci
			) %>%
	unnest() %>%
	write_tsv("data/rf/classification_cv_test_compare.tsv")


# compare models constructed w/ and w/o SCFA data
model_data %>%
	mutate(scfa = str_detect(input, "scfa"), microbiome=str_replace(input, "scfa_", ""), microbiome=str_replace(microbiome, "_scfa", "")) %>%
	filter(microbiome != "scfa") %>%
	select(class, microbiome, scfa, cv_aucs, test_aucs) %>%
	group_by(class, microbiome) %>%
	nest() %>%
	mutate(test = map(data, ~tidy(wilcox.test(.x$test_aucs~as.factor(.x$scfa), exact=FALSE, alternative="less"))),
				wo_lci = map(data, ~quantile(.x$test_aucs[!.x$scfa], probs=0.025)),
				wo_q1 = map(data, ~quantile(.x$test_aucs[!.x$scfa], probs=0.25)),
				wo_median = map(data, ~median(.x$test_aucs[!.x$scfa])),
				wo_q3 = map(data, ~quantile(.x$test_aucs[!.x$scfa], probs=0.75)),
				wo_uci = map(data, ~quantile(.x$test_aucs[!.x$scfa], probs=0.975)),
				w_lci = map(data, ~quantile(.x$test_aucs[.x$scfa], probs=0.025)),
				w_q1 = map(data, ~quantile(.x$test_aucs[.x$scfa], probs=0.25)),
				w_median = map(data, ~median(.x$test_aucs[.x$scfa])),
				w_q3 = map(data, ~quantile(.x$test_aucs[.x$scfa], probs=0.75)),
				w_uci = map(data, ~quantile(.x$test_aucs[.x$scfa], probs=0.975))
			) %>%
	select(class, microbiome, test,
				wo_lci, wo_q1, wo_median, wo_q3, wo_uci,
				w_lci, w_q1, w_median, w_q3, w_uci
			) %>%
	unnest() %>%
	write_tsv("data/rf/classification_w_wo_SCFA.tsv")



# compare SCFA model performance to 0.5
model_data %>%
	filter(input == "scfa") %>%
	group_by(class) %>%
	summarize(test=list(tidy(wilcox.test(test_aucs, mu=0.5, alternative="great"))),
						summary=list(tidy(summary(test_aucs)))) %>%
	unnest() %>%
	write_tsv("data/rf/classification_SCFA_to_random.tsv")

# they're significantly above 0.50, but not by much
