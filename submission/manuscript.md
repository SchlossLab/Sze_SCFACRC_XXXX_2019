---
title: "Revisiting the Relationship between Short-Chain Fatty Acids, the Microbiota, and Colorectal Tumors"
bibliography: references.bib
output:
  pdf_document:
    keep_tex: true
    includes:
      in_header: header.tex
csl: mbio.csl
fontsize: 11pt
geometry: margin=1.0in
---




\vspace{35mm}

Fecal short chain fatty acids are not predictive of colorectal cancer status and cannot be predicted based on bacterial community structure

\vspace{35mm}


Marc A. Sze${^1}$, Begüm Topçuoğlu${^1}$, Nicholas A. Lesniak${^1}$, Mack T. Ruffin IV${^2}$, Patrick D. Schloss${^1}$${^\dagger}$

\vspace{40mm}

$\dagger$ To whom correspondence should be addressed: pschloss@umich.edu

$1$ Department of Microbiology and Immunology, University of Michigan, Ann Arbor, MI 48109

$2$ Department of Family Medicine and Community Medicine, Penn State Hershey Medical Center, Hershey, PA


\newpage
\linenumbers


## Abstract (250 words)

Something clever


## Importance (150 words)

Something clever

\newpage

Colorectal cancer is the third leading cancer-related cause of death within the United States [@Haggar2009; @Siegel2016]. Less than 10% of cases can be attributed to genetic risk factors [@Fearon1990]. This leaves a significant role for environmental,  behavioral, and other factors such as smoking and diet [@FlissIsakov2017; @Lee2015]. Colorectal cancer is thought to be initiated by a series of mutations that accumulate as mutated cells proliferate leading to adenomatous lesions, which are succeeded by carcinomas [@Fearon1990]. Throughout this progression, there are ample opportunities for bacterial populations to create mutations, induce inflammation, and accelerate tumorigenesis [@Flynn2016]. Numerous studies in murine models have supported this model [@Zackular2013; @Baxter2014; @Zackular2015; @DeStefanoShields2016; @Tomkovich2017]. Additional cross sectional studies in humans have identified microbiome-based biomarkers of disease [@Kostic2011; @Sze2018; @Baxter2016; @Zackular2014; @Hannigan2017; @Zeller2014; @Shah2017; @Dejea2018; @Arthur2012]. These studies suggest that in some cases, it is the loss of bacterial populations that produce short-chain fatty acids (SCFAs) that results in increased inflammation and tumorigenesis.

SCFAs have have anti-inflammatory and anti-proliferative activities [@Encarnao2018; @Verma2018; @Zheng2017; @OKeefe2016]. Furthermore, manipulation of SCFAs in mouse models of colorectal cancer by direct supplementation or feeding of fiber caused an overall reduction in tumor burden [@Bishehsari2018; @Tian2018]. These results suggest that supplementation with substrates that bacteria can ferment to produce SCFAs may confer beneficial effects against colorectal cancer. Regardless, there is a lack of evidence that increasing SCFA concentrations can protect against colorectal cancer in humans. Based on similar observations, many microbiome studies use the concentrations of SCFAs and the presence of 16S rRNA gene sequences from organisms and the genes involved in producing them as a biomarker of a healthy microbiota [@Vital2014; @Sanna2019; @Liu2019; @Meisel2016]. Case-control studies that have investigated SCFA concentrations in colorectal cancer found that patients with carcinomas had lower concentrations of SCFAs versus patients with adenomas or individuals without colon tumors [@Ohigashi2013]. Although this would argue that increasing SCFA concentrations could be protective against tumorigenesis, in randomized controlled trials fiber supplementation has been inconsistently associated with protection against tumor formation and recurrence [@Schatzkin2000; @Yao2017; @Kunzmann2015; @Murphy2012; @Gianfredi2018]. These findings temper enthusiasm for treatments that aim to use SCFAs as biomarkers or protection against tumorigenesis.




**SCFA concentrations do not meaningfully vary with diagnosis or treatment.** To quantify the associations between colorectal cancer, the microbiome, and SCFAs, we quantified the concentration of acetate, propionate, isobutyrate, and butyrate in feces of previously characterized individuals with normal colons (N=172) and those with colonic adenomas (N=198) or carcinomas (N=120) [@Baxter2016]. The only SCFA that had a significantly different concentration across the diagnoses was isobutyrate (P=0.0091; Figure 1A). The median concentration of isobutyrate was 3.30 mmol/kg in people with normal colons and it was 3.00 and 3.84 mmol/kg in people with colonic adenomas or carcinomas. The difference in isobutyrate concentration between people with adenomas and carcinomas was significantly different (P=0.0065); however, the differences between people with normal colons and those with adenomas or carcinomas was not significant (P=0.19 and P=0.11). The median concentration of isobutyrate in people with normal colons (3.30 mmol/kg) was between those with adenomas (3.00 mmol/kg) and carcinomas (3.84 mmol/kg). Among the subjects with adenomas and carcinomas, a subset (N~adenoma~=41, N~carcinoma~=26) were treated and sampled a year later [@Sze2017]. The only SCFA that changed following treatment was isobutyrate, which decreased by 0.99 mmol/kg (P=0.002; Figure 1B). For both the pre-treatment cross-sectional data and the pre/post treatment data, we pooled the SCFA concentrations on a per molecule of carbon basis and again failed to see any significant differences (P>0.15). The low concentration of isobutyrate relative to the other SCFAs, inconsistent concentrations, and unexpected decrease in concentration with treatment makes it difficult to ascribe much biological relevance to this observation.



**Combining SCFA and microbiome data does not improve the ability to diagnose individual as having adenomas or carcinomas.** We previously found that binning 16S rRNA gene sequence data into operational taxonomic units based on 97% similarity or into genera enabled us to classify individuals as having adenomas or carcinomas using Random Forest machine learning models. We repeated that analysis but added the concentration of the individual or pooled SCFAs as possible features to train a model. A model trained on the individual or pool SCFAs alone had a median area under the receiver operator characteristic curve (AUROC) of **0.XX** (IQR=**0.XX**). When we trained the models with the SCFAs concentrations and OTU or genus-level relative abundances (AUROC~OTU+SCFA~=**0.XX**, AUROC~genus+SCFA~=**0.XX**), the AUROC values were no different from the models trained without the SCFA concentrations (AUROC~OTU~=**0.XX**, AUROC~genus~=**0.XX**; Figure 2A). We also trained models using a smaller dataset that generated shotgun metagenomic sequencing data from a subset of our cohort (N~normal~=**XX**, N~adenoma~=**XX**, and N~cancer~=**XX**). We binned genes extracted from the assembled metagenomes into operational protein families (OPFs) or KEGG categories [**XXXXXX**]. Again, the performance of the models trained with the meteagenomic data did not improve when the SCFA concentrations were added as possible features when training the model (Figure 2A). These data demonstrate that knowledge of the SCFA profile from a patient's fecal sample does not improve the ability to diagnose a colonic lesion.



**Knowledge of microbial community structure does not predict SCFA concentrations.** Regardless of a person's diagnosis, we next asked whether the fecal community structure was predictive of fecal SCFA concentrations. We trained Random Forest regression models using 16S rRNA gene sequence data binned into OTUs and genera to predict the concentration of the individual and pooled SCFAs. Regardless of the binning method or SCFA, the largest amount of variation that the model could explain was **0.XXXX** (Figure 2B). Next, we trained Random Forest regression models using metagenomic sequence data binned into OPFs and KEGG categories to predict the concentration of the individual and pooled SCFAs. Similar to the analysis using 16S rRNA gene sequence data, the metagenomic data was not predictive of SCFA concentration; the largest amount of variation that the models could explain was **0.XXXX** (Figure 2B). Because of the limited number of samples that we were able to generate metagenomic sequence data from, we used our 16S rRNA gene sequence data to impute metagenomes that were binned into metabolic pathways, enzyme commission numbers, or KEGG categories using picrust2. Again, SCFA concentrations could not be predicted based on the imputed metagenomic data; the largest amount of variation that the models could explain was **0.XXXX** (Figure 2B). The inability to model SCFA concentrations from microbiome data indicates that the knowledge of the abundance of organisms and genes is insufficient to predict SCFA concentrations.


**Conclusion.**

* Unlikely that SCFAs are the primary mechanism that limits tumorigenesis
* Limitations
	- Observing the community once the tumor is developed
	- Not looking at mucosal communities
* Need to identify other possible mechanisms that drive and prevent tumorigenesis
* Supports evidence that fiber consumption has a null effect on tumorigenesis


\newpage

## Materials and Methods

**Study design and sampling.** The overall protocol has been described in detail previously [@Sze2017; @Baxter2016]. In brief, this study used fecal samples obtained at either a single cross-sectional time point (n=**r length(metai$sample)**) or from before (pre-) and after (post-) treatment of a patient's tumor (adenoma n =**r length(filter(metaf, dx == "adenoma") %>% select(initial) %>% pull())** and carcinoma n = ** length(filter(metaf, dx == "cancer") %>% select(initial) %>% pull())**). For patients undergoing treatment for their tumor the length of time between their initial and follow up sample ranged from **r min(metaf$time)** - **r max(metaf$time)** days. Our use of treatment has been previously defined as encompassing removal of a tumor (surgery or colonoscopy) with or without chemotherapy and radiation [@Sze2017]. Diagnosis of tumor was made by colonoscopic examination and histopathological review of biopsies obtained [@Sze2017; @Baxter2016]. The University of Michigan Institutional Review Board approved the study and informed consent was obtained from all participants in accordance to the guidelines set out by the Helsinki Declaration.

**Measuring specific SCFAs.** Our protocol for the measurement of acetate, butyrate, and propionate followed a previously published protocol that used a High-Performance Liquid Chromatography (HPLC) machine [@Venkataraman2016]. The following changes to this protocol included the use of frozen fecal samples suspended in 1ml of PBS instead of fecal suspensions in DNA Genotek OmniGut tubes, and the use of the actual weight of fecal samples instead of the average weight for SCFA concentration normalizations. These methodological changes did not affect the overall median concentrations of these SCFAs between the two studies (see Table 1 [@Venkataraman2016] and Figure 1 here).    

**16S rRNA gene sequencing.** The workflow and processing have been previously described [@Schloss2009; @Kozich2013; @Sze2017]. In brief, sequences were quality filtered and contigs created from the paired end reads. Any sequences with ambiguous base calls were discarded. Contigs were then checked for matches to the V4 region of the 16S rRNA gene using the SILVA database [@Quast2012]. Chimeras were identified and removed using UCHIME and OTUs clustered at 97% similarity [@Edgar2011]. The major differences from these previous reports include: the use of version 1.39.5 of the mothur software package and clustering Operational Taxonomic Units (OTUs) at 97% similarity using the OptiClust algorithm [@Westcott2017].  

**Generating imputed metagenomes.** The use of PICRUSt version 1.1.2 with the recommended standard operating protocol [@Langille2013] was used. Briefly, the mothur shared file and metadata was converted into a biom formatted table using the biom convert function, the subsequent biom file was processed with the 'normalize_by_copy_number.py' function, and subsequent imputed metagenomes created using the 'predict_metagenomes.py' function.

**Obtaining Operational Protein Families from metagenomes.** A subset of the cross-sectional group (n=**r length(metai$sample)**) containing a total of **r length(geof_samples$X1)** individuals (normal n=**r filter(geof_samples, X30 == "Healthy") %>% pull(X30) %>% length()**, adenoma n=**r filter(geof_samples, X30 == "Adenoma") %>% pull(X30) %>% length()**, and carcinoma n=**r filter(geof_samples, X30 == "Cancer") %>% pull(X30) %>% length()**) was shotgun sequenced on an Illumina HiSeq using 125 bp paired end reads and a previously described method [@Hannigan2017]. Briefly, the sequences were quality filtered and sequences aligning to the human genome were removed prior to contig assembly with MEGAHIT [@Li2015]. Open Reading Frames (ORFs) were identified using Prodigal [@Hyatt2010], counts generated using Diamond [@Buchfink2014], subsequent clustering into Operational Protein Families (OPFs) used mmseq2 [@Steinegger2017], and OPF alignment used the KEGG database [@Kanehisa2015].

**Random Forest models.** The model was first trained on 80% of the data and then tested on the held out 20% (80/20 split) using the Random Forest algorithm for classification and regression models via the caret package [@Liaw2002; @Kuhn2017]. This was repeated on 100 different 80/20 splits of the data to generate a reasonable range for the AUC of the model. The reported AUCs, unless otherwise specified, are for the test sets. The classification models were built to group normal versus adenoma, normal versus carcinoma, and high versus low SCFA concentrations. The regression models were built to classify the SCFA concentrations of acetate, butyrate, and propionate regardless of disease status.

**Statistical analysis workflow.** All analysis was performed using the statistical language R [@r_citation_2017]. Generally, a Kruskal-Wallis rank sum test with a Dunn's post-hoc test was used to assess differences between individuals without colon tumors, patients with adenomas, and patients with carcinomas. Where appropriate Benjamini-Hochberg was used to correct for multiple comparisons [**XXXXXX**]. First, we assessed differences in SCFA concentrations measured by HPLC between individuals with normal colons and patients with tumors (adenoma or carcinoma). We then analyzed whether SCFA concentrations changed in patients with an adenoma or carcinoma pre- versus post-treatment. Next, the imputed gene counts of genes encoding enzymes involved in SCFA synthesis was tested. Additionally, metagenomic sequencing counts for important genes involved with SCFA production were analyzed. From here we analyzed the number of significant positive and negative correlations between OTU relative abundance and SCFA concentrations in individuals without colon tumors and patients with adenomas or carcinomas using Spearman's rho. Next, we assessed whether OTUs alone or OTUs and SCFAs were better able to classify individuals with and without tumors using Random Forest models. Finally, models to classify the actual SCFA concentration or high/low SCFA concentration based on the median of each SCFA using 16S rRNA gene sequencing data was created using the Random Forest algorithm. For all Random Forest models, the assessment of the most important variables was based on the top 10 features (OTUs or SCFAs) using the mean decrease in accuracy.

\newpage

## Acknowledgements

The authors thank the Great Lakes-New England Early Detection Research Network for providing the fecal samples that were used in this study. We would thank the University of Michigan Center for Microbial Systems for enabling our short-chain fatty acid analysis. Support for MAS came from the Canadian Institute of Health Research and the National Institutes of Health (UL1TR002240). Support for PDS came from the National Institutes of Health (P30DK034933 and R01CA215574).   

\newpage

## References

<div id="refs"></div>

\newpage

\includegraphics{../results/figures/scfa_abundance.pdf}

**Figure 1. SCFA concentrations did not vary meaningfully with diagnosis of colonic lesions or with treatment for adenomas or carcinomas.** (A) We measured the concentration of fecal SCFAs from individuals with normal colons (N=172) or those with adenoma (N=198) or carcinomas (N=120) was for isobutyrate. (B) A subset of individuals diagnosed with adenomas (N=41) or carcinomas (N=26) who underwent treatment were resampled a year after the initial sampling; one extreme propionate value (124.4 mmol/kg) was included in the adenoma analysis but censored from the visualization for clarity.

\newpage

**Figure 2. SCFA concentrations do not improve models for diagnosing the presence of adenomas, carcinomas, or all lesions. 16S rRNA gene and metagenomic sequence data do not predict SCFAs concentrations.**


\newpage


**Figure S1. Comparison of training and testing results.**
