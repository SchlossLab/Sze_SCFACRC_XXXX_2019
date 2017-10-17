REFS = data/references
FIGS = results/figures
TABLES = results/tables
PROC = data/process
FINAL = submission/

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
