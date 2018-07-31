## Revisiting the relationship between short-chain fatty acids, the microbiota, and colorectal tumors

**Background.** Colorectal cancer (CRC) is a continuing health concern with the majority of the
risk for developing disease being due to environmental factors. The microbiota is one of these
environmental factors with certain bacterial community members being associated with CRC and
other taxa being associated to colons without tumors. Some of the bacterial species in taxa
associated to colons without tumors can use fiber to produce short-chain fatty acids (SCFAs) that can inhibit tumor growth in model systems. However, the data supporting the importance of SCFAs in human CRC is less certain. Here, we test the hypothesis that SCFA concentrations are different in individuals with colorectal tumors.


**Methods.** We analyzed a cross-sectional (n=490) and longitudinal pre- and post-treatment (n=67) group for their concentrations of acetate, butyrate, and propionate. Analysis also included imputed gene relative abundance with PICRUSt, metagenomic sequencing on a subset (n=85) of the total cross-sectional group, and tumor classification and SCFA prediction models using Random Forest.



**Results.** No difference in SCFA concentrations were found between individuals without tumors and patients with adenomas or carcinomas (P-value > 0.15). Using metagenomic sequencing, there was also no difference in genes involved with SCFA synthesis between individuals without tumors and patients with adenomas or carcinomas (P-value > 0.70). Finally, there was no difference between classification models with or without SCFAs in their ability to predict patients with adenomas or carcinomas versus individuals without tumors (P-value > 0.05).



**Conclusions.** Although our data does not support the hypothesis that SCFAs are different in
individuals that have colorectal tumors, there may be context specific scenarios where SCFAs may still be beneficial for treatment of CRC. Alternatively, our observations also support the hypothesis that there may be other mechanisms that have not been thoroughly investigated that are more involved with the development of human CRC.




### Overview

	project
	|- README          # the top level description of content (this doc)
	|- CONTRIBUTING    # instructions for how to contribute to your project
	|- LICENSE         # the license for this project
	|
	|- submission/
	| |- study.Rmd    # executable Rmarkdown for this study, if applicable
	| |- study.md     # Markdown (GitHub) version of the *.Rmd file
	| |- study.tex    # TeX version of *.Rmd file
	| |- study.pdf    # PDF version of *.Rmd file
	| |- header.tex   # LaTeX header file to format pdf version of manuscript
	| |- references.bib # BibTeX formatted references
	| |- XXXX.csl     # csl file to format references for journal XXX
	|
	|- data           # raw and primary data, are not changed once created
	| |- references/  # reference files to be used in analysis
	| |- raw/         # raw data, will not be altered
	| |- mothur/      # mothur processed data
	| +- process/     # cleaned data, will not be altered once created;
	|                 # will be committed to repo
	|
	|- code/          # any programmatic code
	|
	|- results        # all output from workflows and analyses
	| |- tables/      # text version of tables to be rendered with kable in R
	| |- figures/     # graphs, likely designated for manuscript figures
	| +- pictures/    # diagrams, images, and other non-graph graphics
	|
	|- exploratory/   # exploratory data analysis for study
	| |- notebook/    # preliminary analyses
	| +- scratch/     # temporary files that can be safely deleted or lost
	|
	+- Makefile       # executable Makefile for this study, if applicable


### How to regenerate this repository

#### Dependencies and locations
* Gnu Make should be located in the user's PATH
* mothur (v1.39.5) should be located in the user's PATH
* R (v. 3.4.4) should be located in the user's PATH
* PICRUSt (v. 1.1.2) should be located on the user's computer.


#### Running analysis

```
git clone https://github.com/SchlossLab/Sze_scfa_crc_XXXX_2017.git
make write.paper
```
