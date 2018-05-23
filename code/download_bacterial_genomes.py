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
	# gets all bacteria that are contained in the current database
	ftp.login()
	ftp.cwd("genomes/%s/bacteria" % (data_base_to_use))
	dirs = ftp.nlst()
	# Open txt file to write to
	temp_file = open("%sbacteria_in_database.txt" % (output_dir), 'w')
	# Output message to stdout
	print("Identifying all bacterial genomes in %s database..." % 
		(data_base_to_use))
	# Create the needed list of all bacteria present
	for genome in dirs:

		temp_file.write(genome+'\n')

	sys.exit()


# Function to run the actual download
def download_from_ncbi(bacteria_needed, data_base_to_use, file_needed, output_dir):

	#initialize FTP server 
	ftp_site = 'ftp.ncbi.nlm.nih.gov'
	ftp = FTP(ftp_site)

	# set up system calls for downloads of specific file types
	if file_needed == "fasta":
		file_type = "--reject-regex='*cds_from*' --reject-regex='*rna_from*' --accept-regex='*genomic.fna.gz'"
	elif file_needed == "genbank":
		file_type = "--accept-regex='*genomic.gbff.gz'"
	elif file_needed == "gff":
		file_type = "--accept-regex='*genomic.gff.gz'"
	elif args.type == "feature_table":
		file_type = "--accept-regex='*feature_table.txt.gz'"
	# Checks if the default to download all bacterial genomes is set
	if bacteria_needed == "ALL":

		print("Downloading all bacterial genomes from %s database..." % 
			(data_base_to_use))
		# Attempts to dowload all bacterial genomes and throws exception if it fails
		try:
			os.system("wget -r -np %s ftp://ftp.ncbi.nlm.nih.gov/genomes/%s/bacteria/*/latest_assembly_versions/*/ %s" % 
				(file_type, data_base_to_use, output_dir))
		except:
			raise Exception("Failed to download files")
	# If a specific bacterium is needed downloads this from given database
	else:
		# login into the ftp server
		ftp.login()
		ftp.cwd('genomes/%s/bacteria' % (data_base_to_use))
		dirs_list = ftp.nlist()
		# finds any matches within the given database
		genome_matches = fnmatch.filter(dirs_list, bacteria_needed)
		# Exits program if there is no matches found
		if len(genome_matches) == 0:

			print("%s not found in %s database. Exiting." % 
				(bacteria_needed, data_base_to_use))

			sys.exit(1)

		print("Downloading %s genome from %s database..." % 
			(bacteria_needed, data_base_to_use))
		# If there is a match then trys to downloads genome from desired database or throws an error
		try:
			subprocess.call("wget -r -np %s ftp://ftp.ncbi.nlm.nih.gov/genomes/%s/bacteria/%s/latest_assembly_versions/*/ %s" % 
				(file_type, data_base_to_use, bacteria_needed, output_dir), shell=True)
		except:
			raise Exception("Failed to download files")
	


# Run the actual program
def main(bacteria, data_base, type_of_file, av_genomes, outputPath):
	# Creates a download directory if it doesn't exist
	make_bacterial_genome_directory(outputPath)
	# Runs a search of all bacterial genomes and returns a txt file of them if 
	# argument is not blank
	if av_genomes != None:

		find_genomes(data_base, av_genomes)

	# Runs the download from the NCBI database
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