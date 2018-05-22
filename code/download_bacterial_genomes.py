#!python

# Code inspired by:
	# github repo https://github.com/ryjohnson09/bacteria_genome_pull.git

# The main role for this code is to download full genomes from NCBI to be used to search 
# subsequent OGUs for closest taxonomy


############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, argparse, subprocess, fnmatch

# Import other code with useful functions that will be used
from ftplib import FTP

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



# Function to find all genomes in the bacteria directory
def find_genomes(data_base_to_use, output_dir):

	#initialize FTP server 
	ftp_site = 'ftp.ncbi.nlm.nih.gov'
	ftp = FTP(ftp_site)

	ftp.login()
	ftp.cwd("genomes/%s/bacteria" % (data_base_to_use))
	dirs = ftp.nlst()

	for genome in dirs:

		print(genome)

	sys.exit()


# Function to run the actual download
def download_from_ncbi(bacteria_needed, data_base_to_use, file_needed, output_dir):

	#initialize FTP server 
	ftp_site = 'ftp.ncbi.nlm.nih.gov'
	ftp = FTP(ftp_site)

	# set up system calls for downloads of specific file types
	if file_needed == "fasta":
		file_type = "--exclude='*cds_from*' --exclude='*rna_from*' --include='*genomic.fna.gz' --exclude='*'"
	elif file_needed == "genbank":
		file_type = "--include='*genomic.gbff.gz' --exclude='*'"
	elif file_needed == "gff":
		file_type = "--include='*genomic.gff.gz' --exclude='*'"
	elif args.type == "feature_table":
		file_type = "--include='*feature_table.txt.gz' --exclude='*'"

	if bacteria_needed == "ALL":

		try:
			subprocess.call("rsync -Lrtv %s rsync://ftp.ncbi.nlm.nih.gov/genomes/%s/bacteria/*/latest_assembly_versions/*/ %s" % 
				(file_type, data_base_to_use, output_dir), shell=True)
		except:
			raise Exception("Failed to download files")

	else:

		try:
			subprocess.call("rsync -Lrtv %s rsync://ftp.ncbi.nlm.nih.gov/genomes/%s/bacteria/%s/latest_assembly_versions/*/ %s" % 
				(file_type, data_base_to_use, bacteria_needed, output_dir), shell=True)
		except:
			raise Exception("Failed to download files")

# Need to implement this as a fail safe 
#Check that bacterium is in NCBI genbank/refseq
# if args.bacterium:
# 	ftp.login()
# 	ftp.cwd('genomes/%s/bacteria' % NCBI_database)
# 	dirs = ftp.nlst()
# 	pattern = args.bacterium
# 	genome_match = fnmatch.filter(dirs, pattern)
# 	if len(genome_match) == 0:
# 		print("Bacterium %s not found in %s -----> Exiting program" % (args.bacterium, NCBI_database))
# 		sys.exit(1)
# else:
# 	print("")
# 	print("Please provide bacterium name '-b'")
# 	print("Exiting program")
# 	print("")
# 	sys.exit(1)
		


# Run the actual program
def main(bacteria, data_base, type_of_file, av_genomes, outputPath):

	make_bacterial_genome_directory(outputPath)

	if av_genomes != None:

		find_genomes(data_base, av_genomes)

	download_from_ncbi(bacteria, data_base, type_of_file, outputPath)





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