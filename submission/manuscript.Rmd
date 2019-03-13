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


## Abstract



\newpage

## Introduction (400 words)

Colorectal cancer is the third leading cancer-related cause of death within the United States [@Haggar2009; @Siegel2016]. Less than 10% of cases can be attributed to genetic risk factors [@Fearon1990]. This leaves a significant role for environmental,  behavioral, and other factors such as smoking and diet [@FlissIsakov2017; @Lee2015]. Colorectal cancer is thought to be initiated by a series of mutations that accumulate as mutated cells proliferate leading to adenomatous lesions, which are succeeded by carcinomas [@Fearon1990]. Throughout this progression, there are ample opportunities for bacterial populations to create mutations, induce inflammation, and accelerate tumorigenesis [@Flynn2016]. Numerous studies in murine models have supported this model [@Zackular2013; @Baxter2014; @Zackular2015; @DeStefanoShields2016; @Tomkovich2017]. Additional cross sectional studies in humans have identified microbiome-based biomarkers of disease [@Kostic2011; @Sze2018; @Baxter2016; @Zackular2014; @Hannigan2017; @Zeller2014; @Shah2017; @Dejea2018; @Arthur2012]. These studies suggest that in some cases, it is the loss of bacterial populations that produce short-chain fatty acids (SCFAs) that results in increased inflammation and tumorigenesis.

SCFAs have have anti-inflammatory and anti-proliferative activities [@Encarnao2018; @Verma2018; @Zheng2017; @OKeefe2016]. Furthermore, manipulation of SCFAs in mouse models of colorectal cancer by direct supplementation or feeding of fiber caused an overall reduction in tumor burden [@Bishehsari2018; @Tian2018]. These results suggest that supplementation with substrates that bacteria can ferment to produce SCFAs may confer beneficial effects against colorectal cancer. Regardless, there is a lack of evidence that increasing SCFA concentrations can protect against colorectal cancer in humans. Based on similar observations, many microbiome studies use the concentrations of SCFAs and the presence of 16S rRNA gene sequences from organisms and the genes involved in producing them as a biomarker of a healthy microbiota [@Vital2014; @Sanna2019; @Liu2019; @Meisel2016]. Case-control studies that have investigated SCFA concentrations in colorectal cancer found that patients with carcinomas had lower concentrations of SCFAs versus patients with adenomas or individuals without colon tumors [@Ohigashi2013]. Although this would argue that increasing SCFA concentrations could be protective against tumorigenesis, in randomized controlled trials fiber supplementation has been inconsistently associated with protection against tumor formation and recurrence [@Schatzkin2000; @Yao2017; @Kunzmann2015; @Murphy2012; @Gianfredi2018]. These findings temper enthusiasm for treatments that aim to use SCFAs as biomarkers or protection against tumorigenesis.

To better understand the connection between colorectal cancer, the microbiome, and SCFAs, we quantified the concentration of SCFAs in feces of previously characterized individuals with normal colons, adenomas, and carcinomas and subset of those individuals after they underwent treatment for their lesions [@Baxter2016; @Hannigan2017]. Our goals included: (1) testing whether there was an association between SCFA concentration and tumor status, (2) determining whether SCFA concentrations could be used as biomarkers to improve the detection of colonic lesions, and (3) predicting SCFA concentrations based on the relative abundance of fecal bacteria and their genes.

\newpage


## Results (700 words)

* There is no association between presence of tumors and SCFA concentration
	- Cross sectional
	- Aggregating SCFAs does not improve association with CRC
	- Lack of change pre and post treatment
	- Aggregating SCFAs does not improve association with CRC

* Adding SCFA data to diagnostic models does not improve performance
	- Combining 16S and SCFA does not improve prediction of CRC
	- OTU, genus, OPF, kegg

* Community structure is not predictive of SCFA production
	- Random forest regression
	- 16S
	- picrust
	- shotgun metagenomic data


\newpage

## Discussion (200 words)

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

**Pulling genes involved with SCFA synthesis.** Specific genes located near the end of the pathways involved in the synthesis of acetate, butyrate, and propionate were analyzed for any differences between individuals with normal colons and those with tumors. These genes were based on pathways from KEGG as well as previous research [@Baxter2014; @Kanehisa2015] and a list can be found in the supplemental material [Table S1].

**Random Forest models.** The model was first trained on 80% of the data and then tested on the held out 20% (80/20 split) using the Random Forest algorithm for classification and regression models via the caret package [@Liaw2002; @Kuhn2017]. This was repeated on 100 different 80/20 splits of the data to generate a reasonable range for the AUC of the model. The reported AUCs, unless otherwise specified, are for the test sets. The classification models were built to group normal versus adenoma, normal versus carcinoma, and high versus low SCFA concentrations. The regression models were built to classify the SCFA concentrations of acetate, butyrate, and propionate regardless of disease status.

**Statistical analysis workflow.** All analysis was performed using the statistical language R [@r_citation_2017]. Generally, a Kruskal-Wallis rank sum test with a Dunn's post-hoc test was used to assess differences between individuals without colon tumors, patients with adenomas, and patients with carcinomas. Where appropriate Benjamini-Hochberg was used to correct for multiple comparisons [**XXXXXX**]. First, we assessed differences in SCFA concentrations measured by HPLC between individuals with normal colons and patients with tumors (adenoma or carcinoma). We then analyzed whether SCFA concentrations changed in patients with an adenoma or carcinoma pre- versus post-treatment. Next, the imputed gene counts of genes encoding enzymes involved in SCFA synthesis was tested. Additionally, metagenomic sequencing counts for important genes involved with SCFA production were analyzed. From here we analyzed the number of significant positive and negative correlations between OTU relative abundance and SCFA concentrations in individuals without colon tumors and patients with adenomas or carcinomas using Spearman's rho. Next, we assessed whether OTUs alone or OTUs and SCFAs were better able to classify individuals with and without tumors using Random Forest models. Finally, models to classify the actual SCFA concentration or high/low SCFA concentration based on the median of each SCFA using 16S rRNA gene sequencing data was created using the Random Forest algorithm. For all Random Forest models, the assessment of the most important variables was based on the top 10 features (OTUs or SCFAs) using the mean decrease in accuracy.

\newpage

## Acknowledgements

The authors thank the Great Lakes-New England Early Detection Research Network for providing the fecal samples that were used in this study. We would also like to thank Kwi Kim and Thomas M Schmidt for their help in running the short-chain fatty acid analysis on the High-Performance Liquid Chromatography machine at the University of Michigan. We would also like to thank Ada Hagan for providing valuable feedback on earlier drafts of the manuscript. Salary support for Marc A. Sze came from the Canadian Institute of Health Research and NIH grant UL1TR002240. Salary support for Patrick D. Schloss came from NIH grants P30DK034933 and 1R01CA215574.   

\newpage

## References

<div id="refs"></div>

\newpage

**Figure 1. No differences in SCFA measurements were observed between individuals with normal colons or those with adenoma or carcinomas before treatment (A) or due to treatment (B).** Concentrations of SCFAs did not differ among individuals in the different diagnosis categories (normal=**XXX**, adenoma=**XXX**, carcinoma=**XXX**) except for isobutyrate, which was significantly higher among individuals with cancer than those with adenomas (P=**0.XXXX**), but no different from those with normal colons (all P>**0.XXX**) (A). The change in SCFA concentrations among individuals undergoing treatment who initially had adenoma (N=**XX**) and carcinoma (N=**XX**) was only significant for the isobutyrate concentration among individuals with carcinoma (P=**0.XXX**, all others P>**0.XXX**).

**Figure 2. SCFA concentrations do not improve models for diagnosing the presence of adenomas, carcinomas, or all lesions.**

**Figure 3. 16S rRNA gene and metagenomic sequence data do not predict SCFAs concentrations.**

**Figure S1. Comparison of training and testing results.**
