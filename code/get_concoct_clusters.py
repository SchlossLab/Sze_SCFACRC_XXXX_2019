#!python

#This code will run concoct and create contig clusters based on 
# both the coverage and linkage tables provided. I have decided to use
# a combination approach of both Matt and Geoff. Mostly the concoct part
# will be a merger of the two approaches while the linkage and binning will
# mostly draw from Matt.

# Changes from Geof protocol is using default 500 iterations

############# Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, argparse

# Import other code with useful functions that will be used
from qual_trim import create_samples_to_download

# Set up working directory
workdir = "data/raw/"
refdir = "data/references/"


############################################################################################

# needed files
# "data/raw/all_contigs_1kbto10kb.fasta" 
# "data/raw/overall_coverage_table.tsv"


def create_concoct_clusters(contigfile, coveragefile, outputName):

	os.system("concoct \
--coverage_file %s \
--composition_file %s \
--clusters 500 \
--kmer_length 4 \
--length_threshold 1000 \
--read_length 150 \
--basename %s" % (coveragefile, contigfile, outputName))




# Matt way
#/share/apps/rhel6/concoct/0.4.0/bin/concoct 
#-c 400 --coverage_file ${metagenome}.coveragetable.tsv 
#--composition_file ${metagenome}.final.contigs.1k.10k.fa 
#-b ${metagenome}.concoct_output/

# Incorporate linkage with hierarchical clustering
#perl /share/apps/rhel6/concoct/0.4.0/scripts/ClusterLinkNOverlap.pl 
#--cfile=${metagenome}.concoct_output/clustering_gt1000.csv 
#--lfile=${metagenome}.linkagetable.tsv 
#--covfile=${metagenome}.coveragetable.tsv 
#--ofile=${metagenome}.concoct_output/clustering_gt1000_link.csv

# Bin clustered contigs into separate fasta files
#cd ${metagenome}.concoct_output
#python /home/mljenior/scripts/bin_contigs.py ../{metagenome}.final.contigs.1k.10k.fa clustering_gt1000.csv
#cd binned_contigs
#find . -size -1000k -delete
#for x in *.fasta; do python /mnt/EXT/Schloss-data/bin/seq_stats.py $x > $x.summary; done




# Runs the overall program 
def main(contigfileName, coveragetableName, outputDirectoryName):

	create_concoct_clusters(contigfileName, coveragetableName, outputDirectoryName)

	


# Initializes at the start of the program
if __name__ == '__main__': 

	# Command line argument parser with tags for componenets within the program
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("-ct", "--contig_table", 
		default="%sall_contigs_1kbto10kb.fasta" % (workdir), 
		type=str, help="Overall combined contig file. \
		The default is set to 'all_contigs_1kbto10kb.fasta'.\n")
	parser.add_argument("-cov", "--coverage_table", 
		default="%soverall_coverage_table.tsv" % (workdir), 
		type=str, help="Overall coverage table. \
		The default is set to 'overall_coverage_table.tsv'.\n")
	parser.add_argument("-of", "--output_file", 
		default="%sconcoct_output/" % (workdir), 
		type=str, help="File where concoct ouput goes. \
		The default is set to 'data/raw/concoct_output/'.\n")
	args = parser.parse_args()

	# Runs the main function with the following cmd line arguments ported into it
	main(args.contig_table, args.coverage_table, args.output_file)