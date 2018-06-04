REFS = data/references
FIGS = results/figures
TABLES = data/process/tables
PROC = data/process
FINAL = submission/
RAW = data/raw
PLATES = plate1 plate2 plate3 plate4 plate5 plate6 plate7 plate8
SCFA = acetate butyrate propionate

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
	wget -N http://mothur.org/w/images/1/15/Silva.seed_v123.tgz
	tar xvzf Silva.seed_v123.tgz silva.seed_v123.align silva.seed_v123.tax
	mothur "#get.lineage(fasta=silva.seed_v123.align, taxonomy=silva.seed_v123.tax, taxon=Bacteria);degap.seqs(fasta=silva.seed_v123.pick.align, processors=8)"
	mv silva.seed_v123.pick.align $(REFS)/silva.seed.align
	rm Silva.seed_v123.tgz silva.seed_v123.*

$(REFS)/silva.v4.align : $(REFS)/silva.seed.align
	mothur "#pcr.seqs(fasta=$(REFS)/silva.seed.align, start=11894, end=25319, keepdots=F, processors=8)"
	mv $(REFS)/silva.seed.pcr.align $(REFS)/silva.v4.align

# Next, we want the RDP reference taxonomy. The current version is v10 and we
# use a "special" pds version of the database files, which are described at
# http://blog.mothur.org/2014/10/28/RDP-v10-reference-files/

$(REFS)/trainset14_032015.% :
	wget -N http://www.mothur.org/w/images/8/88/Trainset14_032015.pds.tgz
	tar xvzf Trainset14_032015.pds.tgz trainset14_032015.pds/trainset14_032015.pds.*
	mv trainset14_032015.pds/* $(REFS)/
	rmdir trainset14_032015.pds
	rm Trainset14_032015.pds.tgz

# Grab greengenes reference files for picrust

$(REFS)/gg_13_5_99.align : 
	wget -N http://www.mothur.org/w/images/b/be/GG_13_5_otuMapTable.zip
	unzip GG_13_5_otuMapTable.zip
	mv GG_13_5_otuMapTable/*.txt $(REFS)/
	rmdir GG_13_5_otuMapTable
	rm GG_13_5_otuMapTable.zip
	wget -N http://www.mothur.org/w/images/9/9d/Gg_13_5_99.taxonomy.tgz
	tar xvzf Gg_13_5_99.taxonomy.tgz
	mv gg_13_5_99* $(REFS)/
	rm -r Gg_13_5_99.taxonomy.tgz ._gg_13_5_99.pds.tax ._pds.notes __MACOSX/ pds.notes
	wget -N http://www.mothur.org/w/images/c/cd/Gg_13_5_99.refalign.tgz
	tar xvzf Gg_13_5_99.refalign.tgz
	mv gg_13_5_99.align $(REFS)/
	rm pds.notes Gg_13_5_99.refalign.tgz



################################################################################
#
# Part 2: Run data through mothur
#
#	Process fastq data through the generation of files that will be used in the
# overall analysis.
#
################################################################################


# Run initial mothur 16S Sequecning and OTU Clustering
$(PROC)/final.0.03.subsample.shared\
$(PROC)/final.groups.ave-std.summary\
$(PROC)/final.groups.summary\
$(PROC)/final.rep.count_table\
$(PROC)/final.rep.seqs\
$(PROC)/final.sharedsobs.0.03.lt.ave.dist\
$(PROC)/final.sharedsobs.0.03.lt.dist\
$(PROC)/final.sharedsobs.0.03.lt.std.dist\
$(PROC)/final.taxonomy\
$(PROC)/final.thetayc.0.03.lt.ave.dist\
$(PROC)/final.thetayc.0.03.lt.dist\
$(PROC)/final.thetayc.0.03.lt.std.dist\
$(PROC)/final.shared : code/mothurCluster.batch code/mothur.batch\
					$(PROC)/unmatched.file $(PROC)/unmatched.dist\
					$(PROC)/unmatched.taxonomy $(PROC)/unmatched.fasta\
					data/references/silva.v4.align\
					data/references/trainset14_032015.pds.fasta\
					data/references/trainset14_032015.pds.tax\
					data/process/stability.files
	bash code/mothur.batch
	bash code/mothurCluster.batch


# Run the greengenes formated file mothur 16S Sequecning and OTU Clustering
$(PROC)/gg_final.0.03.biom\
$(PROC)/gg_final.0.03.biom_shared\
$(PROC)/gg_final.0.03.subsample.shared\
$(PROC)/gg_final.groups.ave-std.summary\
$(PROC)/gg_final.groups.summary\
$(PROC)/gg_final.rep.count_table\
$(PROC)/gg_final.rep.seqs\
$(PROC)/gg_final.sharedsobs.0.03.lt.ave.dist\
$(PROC)/gg_final.sharedsobs.0.03.lt.dist\
$(PROC)/gg_final.sharedsobs.0.03.lt.std.dist\
$(PROC)/gg_final.taxonomy\
$(PROC)/gg_final.thetayc.0.03.lt.ave.dist\
$(PROC)/gg_final.thetayc.0.03.lt.dist\
$(PROC)/gg_final.thetayc.0.03.lt.std.dist\
$(PROC)/gg_final.shared : code/gg_mothurCluster.batch code/gg_mothur.batch\
					code/align_metadata_picrust.R\
					$(PROC)/unmatched.file $(PROC)/unmatched.dist\
					$(PROC)/unmatched.taxonomy $(PROC)/unmatched.fasta\
					data/references/silva.v4.align\
					data/references/trainset14_032015.pds.fasta\
					data/references/trainset14_032015.pds.tax\
					data/process/stability.files
	bash code/gg_mothur.batch
	R -e "source('code/align_metadata_picrust.R')"
	bash code/gg_mothurCluster.batch


# Run the picrust prediction algorithm on the gg generated file
$(PROC)/nsti.txt\
$(PROC)/gg_corr_OTU_table.biom\
$(PROC)/gg_corr_normalized_otus.biom\
$(PROC)/predicted_metagenomes.biom : code/run_picrust.sh
	bash code/run_picrust.sh


# Create the needed contig fasta file
$(RAW)/all_contigs.fasta : code/module.sh\
		code/qual_trim.py code/rm_human_seqs.py code/create_contigs.py\
		code/create_master_contig_file.py code/create_contig_length_table.py
	source code/module.sh
	python code/qual_trim.py -s "data/raw/whole_metagenome_samples.txt"
	python code/rm_human_seqs.py -s "data/raw/whole_metagenome_samples.txt"
	python code/create_contigs.py -s "data/raw/whole_metagenome_samples.txt"
	python code/create_master_contig_file.py -s "data/raw/whole_metagenome_samples.txt"
	python code/create_contig_length_table.py -cf "data/raw/all_contigs.fasta"


# Create the needed OPFs with contigs 1kb or greater
$(RAW)/all_contigs_1kbto10kb.fasta\
$(PROC)/orf_gene_alignment.tsv\
$(RAW)/diamond_analysis/orf_abund.tsv : code/remove_short_contigs.py\
		code/concoct_scripts/cut_up_fasta.py code/prodigal_wrap.py\
		code/diamond_wrap.py code/run_create_OPF_abund_table.R\
		$(RAW)/mmseq2_opf_run/resultDB.m8
	python code/remove_short_contigs.py -s "data/raw/whole_metagenome_samples.txt"
	python2 code/concoct_scripts/cut_up_fasta.py -c 10000 -o 0 -m $(RAW)/all_contigs_greater1kb.fasta > $(RAW)/all_contigs_1kbto10kb.fasta
	python code/prodigal_wrap.py
	cp data/raw/mmseq2_opf_run/resultDB.m8 data/process/orf_gene_alignment.tsv
	python code/diamond_wrap.py
	Rscript code/run_create_OPF_abund_table.R $(RAW)/mmseq2_opf_run/clu.tsv $(RAW)/diamond_analysis/orf_abund.tsv


# Set up the variables for SCFA processing
ACETATE = $(foreach S, $(PLATES), $(RAW)/Raw_hplc_files/acetate/$(S))
ACETATE_PLATES = $(addsuffix _scfa_crc_acetate.txt,$(ACETATE))
TR_ACETATE = $(foreach S, $(PLATES), $(RAW)/Raw_hplc_files/acetate/transformed_$(S))
TR_ACETATE_PLATES = $(addsuffix _scfa_crc_acetate.csv,$(TR_ACETATE))

BUTYRATE = $(foreach S, $(PLATES), $(RAW)/Raw_hplc_files/butyrate/$(S))
BUTYRATE_PLATES = $(addsuffix _scfa_crc_acetate.txt,$(BUTYRATE))
TR_BUTYRATE = $(foreach S, $(PLATES), $(RAW)/Raw_hplc_files/acetate/transformed_$(S))
TR_BUTYRATE_PLATES = $(addsuffix _scfa_crc_acetate.csv,$(TR_BUTYRATE))

PROPIONATE = $(foreach S, $(PLATES), $(RAW)/Raw_hplc_files/propionate/$(S))
PROPIONATE_PLATES = $(addsuffix _scfa_crc_acetate.txt,$(PROPIONATE))
TR_PROPIONATE = $(foreach S, $(PLATES), $(RAW)/Raw_hplc_files/acetate/transformed_$(S))
TR_PROPIONATE_PLATES = $(addsuffix _scfa_crc_acetate.csv,$(TR_PROPIONATE))

SCFA_PROCESSING_CODE = $(foreach S, $(PLATES), code/run_$(S)_convert_scfa_concentrations.R)

FINAL_SCFA_DATA = $(foreach S, $(SCFA), $(TABLES)/$(S)_final_data.csv)


# Process the SCFA data
$(TR_ACETATE) $(TR_BUTYRATE)\
$(TR_PROPIONATE) $(FINAL_SCFA_DATA) : $(ACETATE_PLATES) $(BUTYRATE_PLATES)\
$(PROPIONATE_PLATES) $(SCFA_PROCESSING_CODE) code/combine_plates.R
	Rscript code/run_plate1_convert_scfa_concentrations.R
	Rscript code/run_plate2_convert_scfa_concentrations.R
	Rscript code/run_plate3_convert_scfa_concentrations.R
	Rscript code/run_plate4_convert_scfa_concentrations.R
	Rscript code/run_plate5_convert_scfa_concentrations.R
	Rscript code/run_plate6_convert_scfa_concentrations.R
	Rscript code/run_plate7_convert_scfa_concentrations.R
	Rscript code/run_plate8_convert_scfa_concentrations.R
	Rscript code/combine_plates.R
		

# Run the SCFA analysis for cross-sectional and pre/post-treatment
SCFA_KRUSKAL_DATA = $(foreach S, $(SCFA), $(TABLES)/$(S)_kruskal_crc_groups.csv)

$(SCFA_KRUSKAL_DATA)\
$(TABLES)/scfa_treatment_comparison_pvalues.csv : $(FINAL_SCFA_DATA)\
		$(RAW)/metadata/metaI_final.csv $(RAW)/metadata/good_metaf_final.csv\
		code/Run_scfa_analysis.R code/Run_treatment_scfa_analysis.R
	Rscript code/Run_scfa_analysis.R
	Rscript code/Run_treatment_scfa_analysis.R


# Run the random forest analysis for combined Adenoma classifications
$(TABLES)/adn_full_eighty_twenty_splits.csv\
$(TABLES)/adn_full_test_data.csv\
$(TABLES)/adn_full_AUC_model_summary.csv\
$(TABLES)/adn_full_raw_mda_values.csv\
$(TABLES)/adn_full_MDA_Summary.csv : $(PROC)/final.0.03.subsample.shared\
		$(RAW)/metadata/metaI_final.csv $(RAW)/good_metaf_final.csv $(FINAL_SCFA_DATA)\
		code/Run_full_adn_RF_data_setup.R code/Run_full_adn_RF_reference.R\
		code/Run_aggregate_full_adn_model.R
	Rscript code/Run_full_adn_RF_data_setup.R
	for ((i=1;i<=100;i++));
	do
		Rscript code/Run_full_adn_RF_reference.R $i
	done
	Rscript code/Run_aggregate_full_adn_model.R


# Run the random forest analysis for combined Carcinoma classifications
$(TABLES)/crc_full_eighty_twenty_splits.csv\
$(TABLES)/crc_full_test_data.csv\
$(TABLES)/crc_full_AUC_model_summary.csv\
$(TABLES)/crc_full_raw_mda_values.csv\
$(TABLES)/crc_full_MDA_Summary.csv : $(PROC)/final.0.03.subsample.shared\
		$(RAW)/metadata/metaI_final.csv $(RAW)/good_metaf_final.csv $(FINAL_SCFA_DATA)\
		code/Run_full_crc_RF_data_setup.R code/Run_full_crc_RF_reference.R\
		code/Run_aggregate_full_crc_model.R
	Rscript code/Run_full_crc_RF_data_setup.R
	for ((i=1;i<=100;i++));
	do
		Rscript code/Run_full_crc_RF_reference.R $i
	done
	Rscript code/Run_aggregate_full_crc_model.R



# Run the random forest analysis for OTU only Adenoma classifications
$(TABLES)/adn_OTU_only_eighty_twenty_splits.csv\
$(TABLES)/adn_OTU_only_full_test_data.csv\
$(TABLES)/adn_otu_only_AUC_model_summary.csv\
$(TABLES)/adn_otu_only_raw_mda_values.csv\
$(TABLES)/adn_otu_only_MDA_Summary.csv : $(PROC)/final.0.03.subsample.shared\
		$(RAW)/metadata/metaI_final.csv $(RAW)/good_metaf_final.csv $(FINAL_SCFA_DATA)\
		code/Run_OTU_only_adn_RF_data_setup.R code/Run_OTU_only_adn_RF_reference.R\
		code/Run_aggregate_otu_only_adn_model.R
	Rscript code/Run_OTU_only_adn_RF_data_setup.R
	for ((i=1;i<=100;i++));
	do
		Rscript code/Run_OTU_only_adn_RF_reference.R $i
	done
	Rscript code/Run_aggregate_otu_only_adn_model.R


# Run the random forest analysis for OTU only Carcinoma classifications
$(TABLES)/crc_full_eighty_twenty_splits.csv\
$(TABLES)/crc_full_test_data.csv\
$(TABLES)/crc_otu_only_AUC_model_summary.csv\
$(TABLES)/crc_otu_only_raw_mda_values.csv\
$(TABLES)/crc_otu_only_MDA_Summary.csv : $(PROC)/final.0.03.subsample.shared\
		$(RAW)/metadata/metaI_final.csv $(RAW)/good_metaf_final.csv $(FINAL_SCFA_DATA)\
		code/Run_OTU_only_crc_RF_data_setup.R code/Run_OTU_only_crc_RF_reference.R\
		code/Run_aggregate_otu_only_crc_model.R
	Rscript code/Run_OTU_only_crc_RF_data_setup.R
	for ((i=1;i<=100;i++));
	do
		Rscript code/Run_OTU_only_crc_RF_reference.R $i
	done
	Rscript code/Run_aggregate_otu_only_crc_model.R

# Run the AUC model analysis
$(TABLES)/rf_model_ttest_comparisons.csv : $(TABLES)/adn_full_AUC_model_summary.csv\
		$(TABLES)/crc_full_AUC_model_summary.csv\
		$(TABLES)/adn_otu_only_AUC_model_summary.csv\
		$(TABLES)/crc_otu_only_AUC_model_summary.csv code/Run_rf_model_analysis.R
	Rscript code/Run_rf_model_analysis.R

# Run the PICRUSt analysis
$(PROC)/picrust_metadata\
$(TABLES)/specific_scfa_kruskal_picrust_summary.csv\
$(TABLES)/selected_scfa_gene_data.csv\
$(TABLES)/specific_scfa_wilcox_picrust_treatment_summary.csv\
$(TABLES)/selected_scfa_treatment_picrust_gene_data.csv : $(RAW)/metadata/metaI_final.csv\
		$(RAW)/metadata/good_metaf_final.csv $(PROC)/final.shared\
		$(PROC)/predicted_metagenomes.biom\
		code/align_metadata_picrust.R code/run_targeted_scfa_picrust.R
	Rscript code/align_metadata_picrust.R
	Rscript code/run_targeted_scfa_picrust.R


# Get specific OPF SCFA KEGG matches
$(PROC)/select_scfa_opf_matches.tsv : code/kegg_parse.py $(PROC)/scfa_kegg_ids.txt\
		$(PROC)/orf_gene_alignment.tsv
	python code/kegg_parse.py -iko "data/process/scfa_kegg_ids.txt" -gad "data/process/orf_gene_alignment.tsv" -o "data/process/select_scfa_opf_matches.tsv" 


# Run specific SCFA OPF analysis
$(PROC)/select_scfa_opf_matches.tsv\
$(PROC)/opf_shared.tsv\
$(PROC)/sra_meta_conversion.txt : $(TABLES)/select_scfa_opf_kruskal_summary.csv\
		$(TABLES)/select_scfa_opf_data.csv code/run_opf_select_scfa_analysis.R
	Rscript code/run_opf_select_scfa_analysis.R


# Run RF model classifications of SCFA high/low
SCFA_RF_DATA = $(foreach S, $(SCFA), $(TABLES)/$(S)_classification_RF_summary.csv)
SCFA_RF_IMP = $(foreach S, $(SCFA), $(TABLES)/$(S)_imp_otus_classification_RF_summary.csv)

$(SCFA_RF_DATA)\
$(SCFA_RF_IMP) : $(PROC)/picrust_metadata $(RAW)/metadata/good_metaf_final.csv\
		$(PROC)/final.0.03.subsample.shared $(FINAL_SCFA_DATA)\
		code/run_rf_scfa_predictions.R
	Rscript code/run_rf_scfa_predictions.R



################################################################################
#
# Part 3: Figure and table generation
#
#	Run scripts to generate figures and tables
#
################################################################################







################################################################################
#
# Part 4: Pull it all together
#
# Render the manuscript
#
################################################################################


$(FINAL)/study.% : 			\ #include data files that are needed for paper
						$(FINAL)/peerj.csl\
						$(FINAL)/references.bib\
						$(FINAL)/study.Rmd
	R -e 'render("$(FINAL)/study.Rmd", clean=FALSE)'
	mv $(FINAL)/study.knit.md $@
	rm $(FINAL)/study.utf8.md

write.paper : $(TABLES)/table_1.pdf $(TABLES)/table_2.pdf\ #customize to include
				$(FIGS)/figure_1.pdf $(FIGS)/figure_2.pdf\	# appropriate tables and
				$(FIGS)/figure_3.pdf $(FIGS)/figure_4.pdf\	# figures
				$(FINAL)/study.Rmd $(FINAL)/study.md\
				$(FINAL)/study.tex $(FINAL)/study.pdf
