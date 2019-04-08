## Fecal short-chain fatty acids are not predictive of colonic tumor status and cannot be predicted based on bacterial community structure


Colonic bacterial populations are thought to have a role in the development of colorectal cancer with some protecting against inflammation and others exacerbating inflammation. Short-chain fatty acids (SCFAs), including butyrate, have been shown to have anti-inflammatory properties and are produced in large quantities by colonic bacteria which produce SCFAs by fermenting fiber. We assessed whether there was an association between fecal SCFA concentrations and the presence of colonic adenomas or carcinomas in a cohort of individuals using 16S rRNA gene and metagenomic shotgun sequence data. We measured the fecal concentrations of acetate, propionate, and butyrate within the cohort and found that there were no significant associations between SCFA concentration and tumor status. When we incorporated these concentrations into random forest classification models trained to differentiate between people with normal colons and those with adenomas or carcinomas, we found that they did not significantly improve the ability of 16S rRNA gene or metagenomic gene sequence-based models to classify individuals. Finally, we generated random forest regression models trained to predict the concentration of each SCFA based on 16S rRNA gene or metagenomic gene sequence data from the same samples. These models performed poorly and were able to explain at most 14% of the observed variation in the SCFA concentrations. These results support the broader epidemiological data that questions the value of fiber consumption for reducing the risks of colorectal cancer. Although other bacterial metabolites may serve as biomarkers to detect adenomas or carcinomas, fecal SCFA concentrations have limited predictive power.





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
* picrust v.1.53.0 (in picrust1 conda environment)
* picrust v.2.1.0-b (in picrust2 conda environment)
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
