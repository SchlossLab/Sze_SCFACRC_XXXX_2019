#!python

# Code inspired by:
	# github repo https://github.com/ryjohnson09/bacteria_genome_pull.git

# The main role for this code is to download full genomes from NCBI to be used to search 
# subsequent OGUs for closest taxonomy


############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, argparse, subprocess, fnmatch, glob

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
	# Download data on all genomes in database
	os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/%s/bacteria/assembly_summary.txt -P %s" % 
		(data_base_to_use, output_dir))
	# Read in the downloaded file
	temp_file = open("%sassembly_summary.txt" % (output_dir), 'r')
	# Create an empty storage dictionary
	temp_link_list = {}
	# Create a counter
	x = 1

	print("Creating dictionary of %s database to download..." % (data_base_to_use))
	# Parse the text file and populate the dictionary
	for line in temp_file:
		# separate the line by tab
		temp_list = line.split('\t')
		# iterate through each part of the list
		for value in temp_list:
			# Find values that are ftp and not in the first two lines
			if "ftp" in value and x > 2:
				# Add to dictionary key = bacterium, value = link
				temp_link_list[temp_list[7]] = value
		# Increase the counter
		x += 1
	# close the file
	temp_file.close()

	# Write the file paths out to at txt file
	temp_path_file = open("%scomplete_genomes.txt" % (output_dir), 'w')

	print("Creating a key txt file...")
	# Create a master txt file of the key value pairs
	for bacterium in temp_link_list:

		temp_path_file.write(bacterium+'\t'+temp_link_list[bacterium]+'\n')
	# close the writing to file
	temp_path_file.close()

	# return the dictionary for later use
	return(temp_link_list)



# Function to search output directory for current genomes and remove those that 
# do not need updating
def get_current_files(data_base_to_use, output_dir):

	print("Getting existing files in %s." %(output_dir))
	# Copy the existing stored dictionary -- to delete from later
	temp_dict = dict(data_base_to_use)
	# create a list of all the files with gz in the output directory
	print("Generating query list...")
	all_files_in_dir = glob.glob("%s*.gz" % (output_dir))

	print("Comparing and removing already downloaded files from list.")
	# loop through the current bacterial database and get the link
	for link in data_base_to_use:

		temp_string = data_base_to_use[link]
		# find and remove the ftp link from string
		if temp_string.startswith('ftp'):

			temp_string = temp_string[55:]
		# loop through all gz files in directory
		for file in all_files_in_dir:
			# remove from the temp dictionary if their is a match
			if temp_string in file:

				del temp_dict[link]
	# return the reduced dictionary
	return(temp_dict)



# Function to run the actual download
def download_from_ncbi(bacteria_needed, data_base_to_use, bacterial_list, output_dir):
		# Execute full download if default ALL is used
		if bacteria_needed == "ALL":

			# Execute the download
			for link in bacterial_list:

				os.system("wget -P %s --reject='*cds_from_genomic.fna.gz' \
--reject='rna_from_genomic.fna.gz' \
--accept='*genomic.fna.gz' %s/*genomic.fna.gz" % 
			(output_dir, bacterial_list[link]))

		# If it is not ALL try to download specific fasta
		else:

			try:

				os.system("wget -P %s --reject='*cds_from_genomic.fna.gz' \
--reject='rna_from_genomic.fna.gz' \
--accept='*genomic.fna.gz' %s/*genomic.fna.gz" % 
			(output_dir, bacterial_list[bacteria_needed]))
			# If it fails pump out this error
			except:

				print("Bacterium not in database. Cannot download...")





# Run the actual program
def main(bacteria, data_base, outputPath):
	# Creates a download directory if it doesn't exist
	make_bacterial_genome_directory(outputPath)
	
	# Runs a search of all bacterial genomes and returns a txt file of them
	# based on database selected
	genome_list = find_genomes(data_base, outputPath)

	# Check to find previously downloaded files and ouput a new genome_list
	temp_file = os.listdir(outputPath)

	if any("gz" in s for s in os.listdir(outputPath)):

		print("Zipped files already exist. Checking downloaded files...")

		genome_list = get_current_files(genome_list, outputPath)

	# Runs the download from the NCBI database
	download_from_ncbi(bacteria, data_base, genome_list, outputPath)





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
	parser.add_argument("-o", "--output_path", 
		default="%sbacterial_genomes/" % (workdir), type=str, 
		help="output directory")


	args = parser.parse_args()
	# Runs the main function with the following cmd line arguments ported into it
	main(args.bacterium, args.database_used, args.output_path)