#!python

# This code is used on aligned sam files to generate the necessary
# files used with concoct. Uses samtools, bedtools, and picard.


############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, re, argparse

# Import other code with useful functions that will be used
from qual_trim import create_samples_to_download

# Set up working directory
workdir = "data/raw/"
refdir = "data/references/"
#all_contigs = "all_contigs.fasta"
#reference = "bowtieReference"
#summarydir = "data/process/tables/"

############################################################################################

# Function to create bam files from sam files
def create_bam_files(sampleList, contigName, outputName):

	# create fai file
	print("samtools faidx %s.fasta" % (contigName))
	
	# Iterate through each individual sam file and convert to bam
	for sam_file in sampleList:
		# convert to bam 
		print("samtools view -bt %s.fasta.fai %s%s_%s_bowtie.sam > %s%s_%s.bam" % 
			(contigName, workdir, sam_file, outputName, workdir, sam_file, outputName))

		# sort the bam file
		print("samtools sort %s%s_%s.bam %s%s_%s_sort.bam" % 
			(workdir, sam_file, outputName, workdir, sam_file, outputName))

		#index the bam file
		print("samtools index %s%s_%s_sort.bam" % 
			(workdir, sam_file, outputName))



# mark duplicates and sort
# get coverage data
# Call specific code for genome coverage table
	#python2 /sw/med/centos7/concoct/0.4.1/bin/gen_input_table.py
# get linkage table

# Runs the overall program 
def main(sampleListFile, contigFile, outputEnding):
	# read in file list from -s call
	samples_to_be_used = create_samples_to_download(sampleListFile)
	create_bam_files(samples_to_be_used, contigFile, outputEnding)



# Upon program call executes these commands automatically
if __name__ == '__main__': 
	# Command line argument parser with tags for componenets within the program
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("-l", "--sample_list", default="%swhole_metagenome_samples.txt" % (workdir), type=str, help="Text file with list of samples\n")
	parser.add_argument("-c", "--contig_file", default="%sall_contigs_1kbto10kb" % (workdir), type=str, help="Combined cut contig fasta fil name\n")
	parser.add_argument("-s", "--sam_file_end", default="contig_1to10kb", type=str, help="Unique Ending of Samfiles\n")
	args = parser.parse_args()
	# Runs the main function with the following cmd line arguments ported into it
	main(args.sample_list, args.contig_file, args.sam_file_end)