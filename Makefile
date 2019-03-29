REFS = data/references
FIGS = results/figures
TABLES = data/process/tables
PROC = data/process
FINAL = submission/
RAW = data/raw
MOTHUR = data/mothur

# utility function to print various variables. For example, running the
# following at the command line:
#
#	make print-BAM
#
# will generate:
#	BAM=data/raw_june/V1V3_0001.bam data/raw_june/V1V3_0002.bam ...
print-%:
	@echo '$*=$($*)'



################################################################################
#
# Part 1: Get the references
#
# We will need several reference files to complete the analyses including the
# SILVA reference alignment and RDP reference taxonomy.
#
################################################################################

# We want the latest greatest reference alignment and the SILVA reference
# alignment is the best reference alignment on the market. This version is from
# v123 and described at http://blog.mothur.org/2015/12/03/SILVA-v123-reference-files/
# We will use the SEED v. 123, which contain 12,083 bacterial sequences. This
# also contains the reference taxonomy. We will limit the databases to only
# include bacterial sequences.
$(REFS)/silva.seed.align :
	wget -N https://mothur.org/w/images/7/71/Silva.seed_v132.tgz
	tar xvzf Silva.seed_v132.tgz silva.seed_v132.align silva.seed_v132.tax
	mothur "#get.lineage(fasta=silva.seed_v132.align, taxonomy=silva.seed_v132.tax, taxon=Bacteria);degap.seqs(fasta=silva.seed_v132.pick.align, processors=8)"
	mv silva.seed_v132.pick.align $(REFS)/silva.seed.align
	rm Silva.seed_v132.tgz silva.seed_v132.*

$(REFS)/silva.v4.align : $(REFS)/silva.seed.align
	mothur "#pcr.seqs(fasta=$(REFS)/silva.seed.align, start=11894, end=25319, keepdots=F, processors=8)"
	mv $(REFS)/silva.seed.pcr.align $(REFS)/silva.v4.align


# Next, we want the RDP reference taxonomy. The current version is v10 and we
# use a "special" pds version of the database files, which are described at
# http://blog.mothur.org/2014/10/28/RDP-v10-reference-files/
$(REFS)/trainset16_022016.% :
	wget -N https://mothur.org/w/images/c/c3/Trainset16_022016.pds.tgz
	tar xvzf Trainset16_022016.pds.tgz trainset16_022016.pds/trainset16_022016.pds.*
	mv trainset16_022016.pds/* $(REFS)/
	rmdir trainset16_022016.pds
	rm Trainset16_022016.pds.tgz


# Grab greengenes reference files for picrust1 - used old version because mappings are available
$(REFS)/97_otu_map.txt $(REFS)/gg_13_5_99.gg.fasta $(REFS)/gg_13_5_99.gg.tax  :
	wget -N http://www.mothur.org/w/images/b/be/GG_13_5_otuMapTable.zip
	unzip GG_13_5_otuMapTable.zip
	mv GG_13_5_otuMapTable/*.txt $(REFS)/
	rm -rf GG_13_5_otuMapTable*
	wget -N http://www.mothur.org/w/images/9/9d/Gg_13_5_99.taxonomy.tgz
	tar xvzf Gg_13_5_99.taxonomy.tgz
	mv gg_13_5_99* $(REFS)/
	rm -rf Gg_13_5_99.taxonomy.tgz ._gg_13_5_99.pds.tax ._pds.notes __MACOSX/ pds.notes
	mv $(REFS)/gg_13_5_99.fasta $(REFS)/gg_13_5_99.gg.fasta


################################################################################
#
# Part 2: Get raw data and curate
#
#	Process data through the generation of files that will be used in the overall
# analyses.
#
################################################################################

# Targets build correctly
# Process the SCFA data
PLATES = plate1 plate2 plate3 plate4 plate5 plate6 plate7 plate8
SCFA = acetate butyrate isobutyrate propionate
PLATE_FILES = $(foreach A, $(SCFA), $(foreach S, $(PLATES),data/raw/scfa/$(A)/$(S)_scfa_crc_$(A).txt))

data/scfa/scfa_composite.tsv : $(PLATE_FILES) data/metadata/scfa_plate_metadata.csv\
									code/scfa_process_plates.R code/scfa_combine_plates.R
	Rscript code/scfa_process_plates.R
	Rscript code/scfa_combine_plates.R


# Targets build correctly
# Download files from SRA, run mothur pipeline to generate files for ASV, phylotype, and OTU-based
# analyses as well as for picrust
$(MOTHUR)/crc.otu.shared $(MOTHUR)/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons_gg.taxonomy $(MOTHUR)/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons_rdp.taxonomy $(MOTHUR)/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta $(MOTHUR)/crc.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.pick.count_table $(MOTHUR)/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.taxonomy :	code/mothur.sh\
				data/references/silva.v4.align\
				data/references/trainset16_022016.pds.fasta\
				data/references/trainset16_022016.pds.tax\
				data/references/gg_13_5_99.gg.fasta\
				data/references/gg_13_5_99.gg.tax
	bash code/mothur.sh


# Targets build correctly
# Pool OTUs into phylotypes
data/phylotype/crc.%.shared data/phylotype/crc.%.taxonomy :\
		$(MOTHUR)/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons_rdp.taxonomy\
		$(MOTHUR)/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.shared
	Rscript code/phylotype.R $*


# Targets build correctly
# Output crc.asv.shared and crc.asv.taxonomy, and crc.asv.fasta file based on screened preclustered
#	sequences in mothur pipeline
data/asv/crc.asv.% : $(MOTHUR)/crc.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.pick.count_table\
			$(MOTHUR)/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.taxonomy\
			code/asv.sh
	bash code/asv.sh


# Targets build correctly
# Run the picrust1 prediction algorithm on the gg generated file
data/picrust1/crc.picrust1.% : code/picrust1.sh\
		code/picrust1_utilities.R\
		$(MOTHUR)/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.shared\
		$(MOTHUR)/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons_gg.taxonomy
	bash code/picrust1.sh


#	Run the picrust2 prediction algorithm on the ASV data to generate crc.ec.shared, crc.ko.shared,
#	and crc.pathways.shared
data/picrust2/crc.%.shared : data/asv/crc.asv.shared\
														data/asv/crc.asv.fasta\
														code/picrust2_utilities.R\
														code/picrust2.sh
	bash code/picrust2.sh


# Targets build correctly
# Create the metag.opf.shared and metag.kegg.shared shared file along wiht the ko/opf look up file
MGSHARED = data/metagenome/metag.opf.shared data/metagenome/metag.kegg.shared
$(MGSHARED) data/metagenome/metag.ko_lookup.tsv : $(REFS)/genes.pep.format.fasta\
																			data/metadata/hannigan_metadata.tsv\
																			code/metagenomics.sh\
																			code/metagenomics_get_shared.R
	bash code/metagenomics.sh

data/metagenome/metag.sample.counts : data/metagenome/metag.kegg.shared code/metagenomics_get_ns.R
	Rscript code/metagenomics_get_ns.R




# Targets build correctly
# Get specific SCFA KEGG matches
SCFA_SHARED=data/metagenome/metag.keggscfa.shared data/picrust1/crc.picrust1scfa.shared data/picrust2/crc.picrust2scfa.shared

.SECONDEXPANSION:
$(SCFA_SHARED) : $$(subst scfa,,$$@) code/scfa_genes.R
	R -e "source('code/scfa_genes.R'); get_scfa_keggs('$<')"

data/metagenome/metag.opfscfa.shared : data/metagenome/metag.opf.shared\
																				data/metagenome/metag.ko_lookup.tsv\
																				code/scfa_genes.R
	R -e "source('code/scfa_genes.R'); get_scfa_keggs('data/metagenome/metag.opf.shared', 'data/metagenome/metag.ko_lookup.tsv')"




################################################################################
#
# Part 3: Compare SCFA concentrations across diagnosis groups
#
################################################################################

# Run the SCFA analysis for cross-sectional and pre/post-treatment
$(PROC)/scfa_cross_section_stats.tsv $(PROC)/scfa_pre_post_stats.tsv : \
										data/scfa/scfa_composite.tsv\
										data/metadata/cross_section.csv\
										data/metadata/follow_up.csv\
										code/scfa_stat_analysis.R
	Rscript code/scfa_stat_analysis.R



################################################################################
#
# Part 5: Build Random Forest algorithms to predict diagnosis groups and SCFA
# concentrations
#
################################################################################

DX = adenoma cancer lesion
SCFA = butyrate isobutyrate propionate acetate pooled
MICROBIOME = asv otu genus family phylum picrust1 pc2ko pc2ec pc2pathways opf kegg
SEED = $(shell seq 0 99)

CLASS_O_RF = $(foreach T,$(1),$(foreach D,$(DX),$(foreach S,$(SEED),data/rf/$D_$T/optimum_mtry.$S.csv)))


ANALYTE_OC=$(call CLASS_O_RF,fit scfa fit_scfa)
ANALYTE_MC=$(subst optimum,all,$(ANALYTE_OC))

$(ANALYTE_OC) $(ANALYTE_MC) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv code/rf_classification.R
	Rscript code/rf_classification.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


ASV_OC=$(call CLASS_O_RF,asv fit_asv scfa_asv fit_scfa_asv)
ASV_MC=$(subst optimum,all,$(ASV_OC))

$(ASV_OC) $(ASV_MC) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv data/asv/crc.asv.shared code/rf_classification.R
	Rscript code/rf_classification.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


OTU_OC=$(call CLASS_O_RF,otu fit_otu scfa_otu fit_scfa_otu)
OTU_MC=$(subst optimum,all,$(OTU_OC))

$(OTU_OC) $(OTU_MC) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv data/mothur/crc.otu.shared code/rf_classification.R
	Rscript code/rf_classification.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


GENUS_OC=$(call CLASS_O_RF,genus fit_genus scfa_genus fit_scfa_genus)
GENUS_MC=$(subst optimum,all,$(GENUS_OC))

$(GENUS_OC) $(GENUS_MC) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv data/phylotype/crc.genus.shared code/rf_classification.R
	Rscript code/rf_classification.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


FAMILY_OC=$(call CLASS_O_RF,family fit_family scfa_family fit_scfa_family)
FAMILY_MC=$(subst optimum,all,$(FAMILY_OC))

$(FAMILY_OC) $(FAMILY_MC) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv data/phylotype/crc.family.shared code/rf_classification.R
	Rscript code/rf_classification.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


PHYLUM_OC=$(call CLASS_O_RF,phylum fit_phylum scfa_phylum fit_scfa_phylum)
PHYLUM_MC=$(subst optimum,all,$(PHYLUM_OC))

$(PHYLUM_OC) $(PHYLUM_MC) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv data/phylotype/crc.phylum.shared code/rf_classification.R
	Rscript code/rf_classification.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


PICRUST1_OC=$(call CLASS_O_RF,picrust1 fit_picrust1 scfa_picrust1 fit_scfa_picrust1)
PICRUST1_MC=$(subst optimum,all,$(PICRUST1_OC))

$(PICRUST1_OC) $(PICRUST1_MC) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv data/picrust1/crc.picrust1.shared code/rf_classification.R
	Rscript code/rf_classification.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


PC2KO_OC=$(call CLASS_O_RF,pc2ko fit_pc2ko scfa_pc2ko fit_scfa_pc2ko)
PC2KO_MC=$(subst optimum,all,$(PC2KO_OC))

$(PC2KO_OC) $(PC2KO_MC) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv data/picrust2/crc.ko.shared code/rf_classification.R
	Rscript code/rf_classification.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


PC2EC_OC=$(call CLASS_O_RF,pc2ec fit_pc2ec scfa_pc2ec fit_scfa_pc2ec)
PC2EC_MC=$(subst optimum,all,$(PC2EC_OC))

$(PC2EC_OC) $(PC2EC_MC) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv data/picrust2/crc.ec.shared code/rf_classification.R
	Rscript code/rf_classification.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


PC2PATHWAYS_OC=$(call CLASS_O_RF,pc2pathways fit_pc2pathways scfa_pc2pathways fit_scfa_pc2pathways)
PC2PATHWAYS_MC=$(subst optimum,all,$(PC2PATHWAYS_OC))

$(PC2PATHWAYS_OC) $(PC2PATHWAYS_MC) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv data/picrust2/crc.pathways.shared code/rf_classification.R
	Rscript code/rf_classification.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


OPF_OC=$(call CLASS_O_RF,opf fit_opf scfa_opf fit_scfa_opf)
OPF_MC=$(subst optimum,all,$(OPF_OC))

$(OPF_OC) $(OPF_MC) : data/scfa/scfa_composite.tsv data/metadata/zackular_metadata.tsv data/metagenome/metag.opf.shared code/rf_classification.R
	Rscript code/rf_classification.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


KEGG_OC=$(call CLASS_O_RF,kegg fit_kegg scfa_kegg fit_scfa_kegg)
KEGG_MC=$(subst optimum,all,$(KEGG_OC))

$(KEGG_OC) $(KEGG_MC) : data/scfa/scfa_composite.tsv data/metadata/zackular_metadata.tsv data/metagenome/metag.kegg.shared code/rf_classification.R
	Rscript code/rf_classification.R $(subst .,,$(suffix $(basename $@))) $(dir $@)



REG_O_RF = $(foreach T,$(1),$(foreach D,$(SCFA),$(foreach S,$(SEED),data/rf/$D_$T/optimum_mtry.$S.csv)))

ANALYTE_OR=$(call REG_O_RF,fit)
ANALYTE_MR=$(subst optimum,all,$(ANALYTE_OR))

$(ANALYTE_OR) $(ANALYTE_MR) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv code/rf_regression.R
	Rscript code/rf_regression.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


ASV_OR=$(call REG_O_RF,asv fit_asv)
ASV_MR=$(subst optimum,all,$(ASV_OR))

$(ASV_OR) $(ASV_MR) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv data/asv/crc.asv.shared code/rf_regression.R
	Rscript code/rf_regression.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


OTU_OR=$(call REG_O_RF,otu fit_otu)
OTU_MR=$(subst optimum,all,$(OTU_OR))

$(OTU_OR) $(OTU_MR) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv data/mothur/crc.otu.shared code/rf_regression.R
	Rscript code/rf_regression.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


GENUS_OR=$(call REG_O_RF,genus fit_genus)
GENUS_MR=$(subst optimum,all,$(GENUS_OR))

$(GENUS_OR) $(GENUS_MR) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv data/phylotype/crc.genus.shared code/rf_regression.R
	Rscript code/rf_regression.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


FAMILY_OR=$(call REG_O_RF,family fit_family)
FAMILY_MR=$(subst optimum,all,$(FAMILY_OR))

$(FAMILY_OR) $(FAMILY_MR) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv data/phylotype/crc.family.shared code/rf_regression.R
	Rscript code/rf_regression.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


PHYLUM_OR=$(call REG_O_RF,phylum fit_phylum)
PHYLUM_MR=$(subst optimum,all,$(PHYLUM_OR))

$(PHYLUM_OR) $(PHYLUM_MR) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv data/phylotype/crc.phylum.shared code/rf_regression.R
	Rscript code/rf_regression.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


PICRUST1_OR=$(call REG_O_RF,picrust1 fit_picrust1)
PICRUST1_MR=$(subst optimum,all,$(PICRUST1_OR))

$(PICRUST1_OR) $(PICRUST1_MR) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv data/picrust1/crc.picrust1.shared code/rf_regression.R
	Rscript code/rf_regression.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


PC2KO_OR=$(call REG_O_RF,pc2ko fit_pc2ko)
PC2KO_MR=$(subst optimum,all,$(PC2KO_OR))

$(PC2KO_OR) $(PC2KO_MR) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv data/picrust2/crc.ko.shared code/rf_regression.R
	Rscript code/rf_regression.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


PC2EC_OR=$(call REG_O_RF,pc2ec fit_pc2ec)
PC2EC_MR=$(subst optimum,all,$(PC2EC_OR))

$(PC2EC_OR) $(PC2EC_MR) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv data/picrust2/crc.ec.shared code/rf_regression.R
	Rscript code/rf_regression.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


PC2PATHWAYS_OR=$(call REG_O_RF,pc2pathways fit_pc2pathways)
PC2PATHWAYS_MR=$(subst optimum,all,$(PC2PATHWAYS_OR))

$(PC2PATHWAYS_OR) $(PC2PATHWAYS_MR) : data/scfa/scfa_composite.tsv data/metadata/cross_section.csv data/picrust2/crc.pathways.shared code/rf_regression.R
	Rscript code/rf_regression.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


OPF_OR=$(call REG_O_RF,opf fit_opf)
OPF_MR=$(subst optimum,all,$(OPF_OR))

$(OPF_OR) $(OPF_MR) : data/scfa/scfa_composite.tsv data/metadata/zackular_metadata.tsv data/metagenome/metag.opf.shared code/rf_regression.R
	Rscript code/rf_regression.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


KEGG_OR=$(call REG_O_RF,kegg fit_kegg)
KEGG_MR=$(subst optimum,all,$(KEGG_OR))

$(KEGG_OR) $(KEGG_MR) : data/scfa/scfa_composite.tsv data/metadata/zackular_metadata.tsv data/metagenome/metag.kegg.shared code/rf_regression.R
	Rscript code/rf_regression.R $(subst .,,$(suffix $(basename $@))) $(dir $@)


# Pool relevant files to each other
.SECONDEXPANSION:
%/cv_test_compare.tsv : $$(foreach S,$$(SEED),$$*/optimum_mtry.$$(S).csv) code/rf_pool_optimum_mtry.R
	Rscript code/rf_pool_optimum_mtry.R $*


# Make monster pool files for classification and RF analysis
CLASSIFICATION_TAGS = fit scfa fit_scfa $(MICROBIOME) $(foreach M,$(MICROBIOME),fit_$M) $(foreach M,$(MICROBIOME),scfa_$M) $(foreach M,$(MICROBIOME),fit_scfa_$M)
CLASSIFICATION_POOLS = $(foreach C,$(CLASSIFICATION_TAGS),$(foreach D,$(DX),data/rf/$D_$C/cv_test_compare.tsv))

data/rf/classification_data_pool.tsv : $(CLASSIFICATION_POOLS) code/rf_pool_pools.R
	Rscript code/rf_pool_pools.R $@ $(CLASSIFICATION_POOLS)

REGRESSION_TAGS = fit $(MICROBIOME) $(foreach M,$(MICROBIOME),fit_$M)
REGRESSION_POOLS = $(foreach C,$(REGRESSION_TAGS),$(foreach S,$(SCFA),data/rf/$S_$C/cv_test_compare.tsv))

data/rf/regression_data_pool.tsv : $(REGRESSION_POOLS) code/rf_pool_pools.R
	Rscript code/rf_pool_pools.R $@ $(REGRESSION_POOLS)


################################################################################
#
# Part 5: Figure and table generation
#
#	Run scripts to generate figures and tables
#
################################################################################

data/rf/classification_cv_test_compare.tsv data/rf/classification_w_wo_SCFA.tsv data/rf/classification_SCFA_to_random.tsv : data/rf/classification_data_pool.tsv code/rf_compare_classification_models.R
	Rscript code/rf_compare_classification_models.R


results/figures/scfa_abundance.pdf : code/plot_scfa_comparisons.R\
																		data/scfa/scfa_composite.tsv\
																		data/metadata/cross_section.csv\
																		data/metadata/follow_up.csv
	Rscript code/plot_scfa_comparisons.R


results/figures/scfa_modeling.pdf : code/plot_classification_regression.R\
																		data/rf/classification_data_pool.tsv\
																		data/rf/regression_data_pool.tsv
	Rscript code/plot_classification_regression.R


results/figures/classification_testing.pdf : code/plot_classification_fit.R\
																						data/rf/classification_data_pool.tsv
	Rscript code/plot_classification_fit.R


results/figures/regression_testing.pdf : code/plot_regression_fit.R\
																						data/rf/regression_data_pool.tsv
	Rscript code/plot_regression_fit.R


################################################################################
#
# Part 6: Pull it all together
#
# Render the manuscript
#
################################################################################

submission/figure_1.ps : results/figures/scfa_comparisons.pdf
	pdf2ps $^ $@

submission/figure_2.ps : results/figures/scfa_modeling.pdf
	pdf2ps $^ $@

submission/figure_s1.ps : results/figures/classification_testing.pdf
	pdf2ps $^ $@

submission/figure_s2.ps : results/figures/regression_testing.pdf
	pdf2ps $^ $@

%.png : %.ps
	convert -density 300 $^ $@

figures : submission/figure_1.ps submission/figure_2.ps\
					submission/figure_s1.ps submission/figure_s2.ps\
					submission/figure_1.png submission/figure_2.png\
					submission/figure_s1.png submission/figure_s2.png




write.paper : \
		figures\
		data/metadata/cross_section.csv\
		data/metadata/follow_up.csv\
		data/process/scfa_cross_section_stats.tsv\
		data/process/scfa_pre_post_stats.tsv\
		data/rf/classification_SCFA_to_random.tsv\
		data/rf/classification_w_wo_SCFA.tsv\
		data/rf/regression_data_pool.tsv\
		data/metagenome/metag.sample.counts\
		$(FINAL)/manuscript.Rmd\
		$(FINAL)/references.bib $(FINAL)/mbio.csl $(FINAL)/header.tex
	R -e 'library(rmarkdown); render("submission/manuscript.Rmd", clean=FALSE)'
	mv submission/manuscript.knit.md submission/manuscript.md
	rm submission/manuscript.utf8.md
