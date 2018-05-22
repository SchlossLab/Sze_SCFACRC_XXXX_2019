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

# Function to create the bacterial genome download directory
def make_bacterial_genome_directory(directory_name):

	# Create a directory to hold diamond analysis
	if os.path.exists(directory_name):

		print("%s directory exists" % (directory_name))

	else:

		print("No %s directory found, creating..." % (directory_name))

		os.system("mkdir %s" % (directory_name))







# Run the actual program
def main(bacteria, data_base, type_of_file, av_genomes, outputPath):

	make_bacterial_genome_directory(outputPath)





# Initialized when program starts
if __name__ == '__main__': 
	# Command line argument parser with tags for componenets within the program
	parser = argparse.ArgumentParser(description=__doc__)
	# Add arguments that will be called from the command line
	parser.add_argument("-b", "--bacterium", 
		default="ALL", 
		type=str, help="Sets the specific bacterium genome to download. \
		The default is to download all bacteria and the argument is set. \
		to 'ALL'. Otherwise entries should be in the form of eg. Klebsiella_pneumoniae.\n")
	parser.add_argument("-d", "--database_used", 
		default="genbank", choices = ["refseq", "genbank"],  
		type=str, help="The database to download the genome(s) from. \
		The default database that is used is set to 'genbank'. \
		The input must be one of 'genbank' or 'refseq'.")
	parser.add_argument("-t", "--file_type", 
		default="fasta",choices = ["genbank", "fasta", "gff", "feature_table"], type=str, 
		help="The file type that will be downloaded from the database. \
		The default is set to fasta. There are total of 4 options: \
		'genbank', 'fasta', 'gff', and 'feature_table'.")
	parser.add_argument("-g", "--genomes", 
		choices=["genbank", "refseq"], type=str, 
		help="This will list all available genomes in specified database. \
		The argument is used in isolation (alone).")
	parser.add_argument("-o", "--output_path", 
		default="%sbacterial_genomes/" % (workdir), type=str, 
		help="output directory")


	args = parser.parse_args()
	# Runs the main function with the following cmd line arguments ported into it
	main(args.bacterium, args.database_used, args.file_type, 
		args.genomes, args.output_path)