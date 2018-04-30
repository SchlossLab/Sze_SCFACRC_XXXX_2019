#!python

#This code will run concoct and create contig clusters based on 
# both the coverage and linkage tables provided. I have decided to use
# a combination approach of both Matt and Geoff. Mostly the concoct part
# will be a merger of the two approaches while the linkage and binning will
# mostly draw from Matt.

# Changes from Geof protocol is using default 500 iterations
# Without threading it takes 372 hours 22 minutes and 59 seconds
	# This is with 1 processor and 46GB


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

# Create a function that generates the contig clusters
def create_concoct_clusters(contigfile, coveragefile, outputName, threading):

	if threading == "Yes":

		os.system("concoct \
--coverage_file %s \
--composition_file %s \
--clusters 500 \
--kmer_length 4 \
--length_threshold 1000 \
--read_length 150 \
--num-threads 10 \
--basename %s" % (coveragefile, contigfile, outputName))

	else:

		os.system("concoct \
--coverage_file %s \
--composition_file %s \
--clusters 500 \
--kmer_length 4 \
--length_threshold 1000 \
--read_length 150 \
--basename %s" % (coveragefile, contigfile, outputName))


# Updates the clusters based on linkage analysis (if available)
def run_linkage_analysis(concoct_file_path, linkage_table, coverage_table):

	# THis throws a divide by 0 error and might be due to only using a subset of the data

	os.system("perl /sw/med/centos7/concoct/0.4.1/bin/ClusterLinkNOverlap.pl \
--cfile=%sclustering_gt1000.csv \
--lfile=%s \
--covfile=%s \
--ofile=%s/clustering_gt1000_link.csv" % 
		(concoct_file_path, linkage_table, coverage_table, concoct_file_path))



# Runs the overall program 
def main(contigfileName, coveragetableName, 
	linkagetableName, outputDirectoryName, threadingCall):

	create_concoct_clusters(contigfileName, coveragetableName, 
		outputDirectoryName, threadingCall)

	#run_linkage_analysis(outputDirectoryName, linkagetableName, coveragetableName)


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
	parser.add_argument("-lin", "--linkage_table", 
		default="%soverall_linkage_table.tsv" % (workdir), 
		type=str, help="Overall linkage table. \
		The default is set to 'overall_linkage_table.tsv'.\n")
	parser.add_argument("-of", "--output_file", 
		default="%sconcoct_output/" % (workdir), 
		type=str, help="File where concoct ouput goes. \
		The default is set to 'data/raw/concoct_output/'.\n")
	parser.add_argument("-t", "--threading", 
		default="Yes", 
		type=str, help="Logical as to whether or not to \
		include threading or not. \
		The default is set to 'Yes'.\n")
	args = parser.parse_args()

	# Runs the main function with the following cmd line arguments ported into it
	main(args.contig_table, args.coverage_table, 
		args.linkage_table, args.output_file, args.threading)