#!python

# This code does a diamond blast on a select set of sequences
# Designed to be used after the opf_seq_picker.py program
# It uses the the kegg database


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


# Function to make a new empty directory if it doesn't exist
def make_new_dir():

	# Create a directory to hold diamond analysis
	if os.path.exists("%skegg_match" % (workdir)):

		print("kegg_match directory exists")

	else:

		print("No kegg_match directory found, creating...")

		os.system("mkdir %skegg_match" % (workdir))


# Function to make uniprot diamond database:
def make_db(ref_fasta):

	print("I am using the database with the path %s." %(ref_fasta))


	os.system("diamond makedb --in %s -d %skegg_match/kegg_reference_db" % 
	 	(ref_fasta, workdir))


	print("Completed creating a kegg_reference_db.")


# Run diamond of each fastq against the created database 
def run_diamond_alignment(input_fasta_file):

	print("Running diamond analysis on file with path %s." % (input_fasta_file))

	# Use blast to get hits of sequences to Uniprot genes
	os.system("diamond blastp -q %s -d %skegg_match/kegg_reference_db.dmnd \
-a %skegg_match/match_output.daa -t %skegg_match/ \
--max-target-seqs 1 --evalue 1e-15 --id 0.90 --threads 8" % 
			(input_fasta_file, workdir, workdir, workdir))

	os.system("diamond view -a %skegg_match/match_output.daa \
-o %skegg_match/match_output.diamondout" % (workdir, workdir))

		# Diamond output file has the following default columns (left to right)
		# qseqid = Query Seq - id, sseqid = Subject Seq - id, 
		# pident = Percentage of identical matches,  length = Alignment length, 
		# mismatch = Number of mismatches, gapopen = Number of gap openings, 
		# qstart = Start of alignment in query, qend = End of alignment in query, 
		# sstart = Start of alignment in subject, send = End of alignment in subject, 
		# evalue = Expect value, bitscore = Bit score

# Send a copy to the data/process directory to be GitHub tracked
def make_match_copy(output_file_name):

	print("Making a copy of results in data/process directory.")

	os.system("cp %skegg_match/match_output.diamondout %s" % (workdir, output_file_name))




# Run the actual program
def main(inputFilePath, curated_kegg, outputName):

	make_new_dir()

	make_db(curated_kegg)

	run_diamond_alignment(inputFilePath)
	
	make_match_copy(outputName)


# Initialized when program starts
if __name__ == '__main__': 
	# Command line argument parser with tags for componenets within the program
	parser = argparse.ArgumentParser(description=__doc__)
	# Add arguments that will be called from the command line
	parser.add_argument("-i", "--input_fasta", 
		default="data/process/sig_dunn_protein.fasta", 
		type=str, help="Fasta file of select proteins to search. \
		The default file used is the one after a Dunn's post hoc test. \
		The file is 'data/process/sig_dunn_protein.fasta'.\n")
	parser.add_argument("-ck", "--curated_kegg", 
		default="%skegg_protein.fasta" % (refdir), 
		type=str, help="Fasta file of KEGG protein fasta. \
		The default file used was downloaded in 2016. \
		The file is '%suniprot_sprot.fasta'.\n" %(refdir))
	parser.add_argument("-o", "--output_file", 
		default="data/process/sig_dunn_kegg_protein_matches.tsv", 
		type=str, help="A tsv file of protein matches from the uniprot database. \
		The default file used for the output to be saved under is \
		'data/process/sig_dunn_protein_matches.tsv'.\n")

	args = parser.parse_args()
	# Runs the main function with the following cmd line arguments ported into it
	main(args.input_fasta, args.curated_kegg, args.output_file)