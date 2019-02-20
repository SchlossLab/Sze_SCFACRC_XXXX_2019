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
$(MOTHUR)/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.shared $(MOTHUR)/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons_gg.taxonomy $(MOTHUR)/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons_rdp.taxonomy $(MOTHUR)/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta $(MOTHUR)/crc.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.pick.count_table $(MOTHUR)/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.taxonomy :	code/mothur.sh\
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









# Run specific SCFA OPF analysis
# Nothing seems to depend on these targets and they are redundant with the previous rule
# $(PROC)/select_scfa_opf_matches.tsv\
# $(PROC)/opf_shared.tsv\
# $(PROC)/sra_meta_conversion.txt : $(TABLES)/select_scfa_opf_kruskal_summary.csv\
# 		$(TABLES)/select_scfa_opf_data.csv code/run_opf_select_scfa_analysis.R
# 	Rscript code/run_opf_select_scfa_analysis.R






# Run the SCFA analysis for cross-sectional and pre/post-treatment
$(PROC)/scfa_cross_section_stats.tsv $(PROC)/scfa_pre_post_stats.tsv : \
										data/scfa/scfa_composite.tsv\
										data/metadata/cross_section.csv\
										data/metadata/follow_up.csv\
										code/scfa_stat_analysis.R
	Rscript code/scfa_stat_analysis.R






################################################################################
#
# Part 3: Figure and table generation
#
#	Run scripts to generate figures and tables
#
################################################################################

$(FIGS)/Figure1.pdf : $(RAW)/metadata/good_metaf_final.csv\
		$(RAW)/metadata/metaI_final.csv $(FINAL_SCFA_DATA)
	Rscript code/make_scfa_measures_hplc_figure.R


$(FIGS)/Figure2.pdf : $(TABLES)/selected_scfa_gene_data.csv\
		$(TABLES)/select_scfa_opf_data.csv code/make_pi_opf_combined_graph.R
	Rscript code/make_pi_opf_combined_graph.R


$(FIGS)/Figure3.pdf : $(PROC)/final.taxonomy $(TABLES)/significant_reg_otu_comp_summary.csv\
		code/make_reg_group_scfa_sig_graph.R
	Rscript code/make_reg_group_scfa_sig_graph.R


$(FIGS)/Figure4.pdf : $(TABLES)/adn_full_AUC_model_summary.csv\
		$(TABLES)/adn_otu_only_AUC_model_summary.csv\
		$(TABLES)/crc_full_AUC_model_summary.csv\
		$(TABLES)/crc_otu_only_AUC_model_summary.csv\
		code/make_rf_auc_graphs.R
	Rscript code/make_rf_auc_graphs.R


$(FIGS)/Figure5.pdf : $(SCFA_RF_R_TRAIN) $(SCFA_RF_R_TEST) $(SCFA_RF_R_DATA)\
		$(SCFA_RF_R_IMP) code/figure5.R
	Rscript code/figure5.R


$(FIGS)/FigureS1.pdf : $(PROC)/final.taxonomy $(TABLES)/significant_class_otu_comp_summary.csv\
		code/make_class_group_scfa_sig_graph.R
	Rscript code/make_class_group_scfa_sig_graph.R

$(FIGS)/FigureS2.pdf : $(SCFA_RF_C_TRAIN) $(SCFA_RF_C_TEST) $(SCFA_RF_C_DATA)\
		$(SCFA_RF_C_IMP) code/figureS2.R
	Rscript code/figureS2.R




################################################################################
#
# Part 4: Pull it all together
#
# Render the manuscript
#
################################################################################

write.paper : $(FINAL)/manuscript.Rmd $(FINAL)/supplement.Rmd\
		$(FINAL)/references.bib $(FINAL)/mbio.csl $(FINAL)/header.tex\
		$(FIGS)/Figure1.pdf $(FIGS)/Figure2.pdf\
		$(FIGS)/Figure3.pdf $(FIGS)/Figure4.pdf\
		$(FIGS)/Figure5.pdf $(FIGS)/FigureS1.pdf\
		$(FIGS)/FigureS2.pdf code/Run_render_paper.R
	R -e "source('code/Run_render_paper.R')"
