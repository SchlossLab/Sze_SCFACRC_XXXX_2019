#!python

# This code is used to create a master bowtie reference file, aligning 
# sequences to contigs, and then using this reference generate 
# contig relative abundnace.
# This task was completed using BowTie2


############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, re

# Import other code with useful functions that will be used
from qual_trim import command_line, create_samples_to_download

# Set up working directory
workdir = "data/raw/"
refdir = "data/references/"
all_contigs = "all_contigs.fasta"
reference = "bowtieReference"

############################################################################################

# Function that creates a bowtie reference file from all contigs
def make_ref_file():

	print("Starting to build reference database")

	os.system("bowtie2-build -q %s%s %s/bowtieReference" % 
		(workdir, all_contigs, workdir))

	print("Completed building reference database")


# Function to align reverse sequences to the master contig reference database
def run_alignment(sampleList):

	for r2_fasta in sampleList:

		os.system("bowtie2 -x %s%s -q %s%s_qf_2.fastq -S %s%s_bowtie.sam -p 8 -L 25 -N 1" % 
			(workdir, reference, workdir, r2_fasta, workdir, r2_fasta))


# Function to create a dictionary with 0 counts for all the contigs 
def create_count_dict():

	print("Generating blank contig count dictionary")

	tempDict = {}

	temp_file = open("%s%s" % (workdir, all_contigs), 'r')

	x = 0

	for line in temp_file:

		if ">" in line:

			test = line.split(" ")[0]

			tempName = test.strip('>')

			tempDict[tempName] = 0

	temp_file.close()

	print("Completed generation of blank contig count dictionary")

	return tempDict



# Function to count specific abundances within the generated sam files
def get_contig_abundance(sampleList, countDict):

	for sample_id in sampleList:

		tempDict = countDict

		temp_file = open("%s%s_bowtie.sam" % (workdir, sample_id), 'r')

		

# Needs to take in a sam file and then write an output file (tsv)
	# Read in sam file
	# Parse file and find /^@[A-Z]{2}/ 
	# separate on the \t and take the second component 
	# find matching contig and add a count to it.
	# Read out as a tsv
# Need to be able to do this for every sample
# Keep as individual tsv files for now


# Runs the overall program 
def main():

	meta_genome_file_name = command_line()
	samples_to_be_used = create_samples_to_download(meta_genome_file_name)
	#make_ref_file()
	#run_alignment(samples_to_be_used)	
	blank_count_dict = create_count_dict()

if __name__ == '__main__': main()