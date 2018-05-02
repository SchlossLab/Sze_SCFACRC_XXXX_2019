#!python

# This code does a diamond blast on a select set of sequences
# It uses the the uniprot database


############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, argparse

# Import other code with useful functions that will be used
from qual_trim import create_samples_to_download

# Set up working directory
workdir = "data/raw/"

# Set up refrence directory
refdir = "data/references/"


############################################################################################



















# Run the actual program
def main(inputFilePath, curated_uniprot, uncurated_uniprot):
	



# Initialized when program starts
if __name__ == '__main__': 
	# Command line argument parser with tags for componenets within the program
	parser = argparse.ArgumentParser(description=__doc__)
	# Add arguments that will be called from the command line
	parser.add_argument("-i", "--input_fasta", 
		default="data/process/sig_dunn_protein.fasta", 
		type=str, help="Fasta file of select proteins to search. \
		The default file used is the one after a Dunn's post hoc test. \
		The file is 'data/process/sig_dunn_protein.fasta.\n" %(workdir))
	parser.add_argument("-cu", "--curated_uni", 
		default="%suniprot_sprot.fasta" % (refdir), 
		type=str, help="Fasta file of uniprot swiss-prot proteins fasta. \
		The default file used was downloaded on May 2, 2018 . \
		The file is '%suniprot_sprot.fasta'.\n" %(refdir))
	parser.add_argument("-uu", "--uncurated_uni", 
		default="%suniprot_trembl.fasta" % (refdir), 
		type=str, help="Fasta file of uniprot swiss-prot proteins fasta. \
		The default file used was downloaded on May 2, 2018 . \
		The file is '%suniprot_trembl.fasta'.\n" %(refdir))

	args = parser.parse_args()
	# Runs the main function with the following cmd line arguments ported into it
	main(args.input_fasta, args.curated_uni, args.uncurated_uni)