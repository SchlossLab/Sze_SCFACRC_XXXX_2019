#!python

#This code will run concoct and create contig clusters based on 
# both the coverage and linkage tables provided

############# Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, argparse

# Import other code with useful functions that will be used
from qual_trim import create_samples_to_download

# Set up working directory
workdir = "data/raw/"
refdir = "data/references/"


############################################################################################


# Geof way
#concoct \
#		--coverage_file ./data/ContigAbundForConcoctVirus.tsv \
#		--composition_file ./data/totalcontigsvirus.fa \
#		--clusters 500 \
#		--kmer_length 4 \
#		--length_threshold 1000 \
#		--read_length 150 \
#		--basename ./data/ContigClustersVirus/ \
#		--no_total_coverage \
#		--iterations 50


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
def main(sampleListFile):

	samples_to_be_used = create_samples_to_download(sampleListFile)
	make_contigs(samples_to_be_used)


# Initializes at the start of the program
if __name__ == '__main__': 

	# Command line argument parser with tags for componenets within the program
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("-s", "--sample_list", 
		default="%swhole_metagenome_samples.txt" % (workdir), 
		type=str, help="Text file with list of samples to run through the program.\n")
	args = parser.parse_args()

	# Runs the main function with the following cmd line arguments ported into it
	main(args.sample_list)