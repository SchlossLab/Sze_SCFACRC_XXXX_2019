## Revisiting the relationship between short-chain fatty acids, the microbiota, and colorectal tumors

**Background.** Colorectal cancer (CRC) is increasing in prevalence in individuals under 50 and
because of this will be a continuing health concern for the foreseeable future. The majority of the risk for developing CRC is attributable to environmental factors. One of these environmental factors is the microbiota, with certain bacterial community members being more prevalent in CRC while other taxa are higher in individuals without tumors. Some of the bacterial species in taxa that are more abundant in individuals without tumors can use fiber to produce short-chain fatty acids (SCFAs) that inhibit tumor growth in model systems. However, the data supporting the importance of SCFAs in human CRC is less certain. Here, we test the hypothesis that SCFA concentrations and taxa associated with their production are different in individuals with colorectal tumors.


**Methods.** We analyzed a cross-sectional (n=490) and longitudinal pre- and post-treatment (n=67) group for their fecal concentrations of acetate, butyrate, and propionate. Analysis also included imputed gene relative abundance with PICRUSt, metagenomic sequencing on a subset (n=85) of the total cross-sectional group, and tumor classification and SCFA prediction models using Random Forest.


**Results.** No difference in SCFA concentrations were found between individuals without tumors
and patients with adenomas or carcinomas (P-value > 0.15). Using metagenomic sequencing, there was also no difference in genes involved with SCFA synthesis between individuals without tumors and patients with adenomas or carcinomas (P-value > 0.70). Finally, there was no difference between the ability of Random Forest models to predict patients with adenomas or carcinomas versus individuals without tumors (P-value > 0.05).


**Conclusions.** Although our data does not support the hypothesis that fecal SCFA concentrations are different in the general CRC population, there still may be specific types of colorectal tumors where SCFAs may be beneficial for treatment of CRC. Alternatively, our observations also support the hypothesis that there may be other metabolites or mechanisms (e.g. bacterial niche exclusion) that may be more protective against tumorigenesis and have not been thoroughly investigated in the context of human CRC.





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
* mothur (v1.42) should be located in the user's PATH
* R (v. 3.X.X) should be located in the user's PATH
* PICRUSt (v. 1.1.3) should be located in the user's PATH.
* mmseqs2
* various modules that are on flux

Dependencies:
* SRA toolkit v.2.8.2 (in path)
* mothur v.1.42.0 (in path)
* picrust v.1.53.0 (in path)
* R packages
	- tidyverse v.1.2.1
* Other
	- Need KEGG's genes.pep.format.fasta in data/references
	- Analysis assumes the use of 12 processors


#### Running analysis

Don't do this, but you'll get the idea...

```
git clone https://github.com/SchlossLab/Sze_SCFACRC_XXXX_2019.git
make write.paper
```
