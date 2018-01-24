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


# Runs the overall program 
def main():

	meta_genome_file_name = command_line()
	samples_to_be_used = create_samples_to_download(meta_genome_file_name)
	#make_ref_file()
	run_alignment(samples_to_be_used)	


if __name__ == '__main__': main()