#!python

# This code is used to generate the total sequences within a given fasta file.
# It assumes that the programs that have been used maintain sequence
# pairing when processing the R1 (forward) and R2 (reverse) fastq files


############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys, re

# Import other code with useful functions that will be used
from qual_trim import command_line, create_samples_to_download

# Set up working directory
workdir = "data/raw/"
refdir = "data/references/"
summarydir = "data/process/tables/"


############################################################################################


# Function to read in a fastq and count unique sequences
# Default would be to collect total and then 
# Sequences after quality filtering steps
def get_seq_counts(sampleList, seq_ending):
	# sampleList is a list of samples to be analyzed
	# seq_ending is the unique ending attached to the fastq file

	# initializes an empty dictionary
	temp_dict = {}
	# iteratively goes through each sample
	for fastq_name in sampleList:
		# Prints a notification to stdoutput
		print("Counting sequences in %s%s" % (fastq_name, seq_ending))
		# Counter
		x = 0
		# Opens the respective fastq sample
		temp_file_qf = open("%s%s%s" % 
			(workdir, fastq_name, seq_ending), 'r')
		# goes through each line and checks for the unique ID 
		for line in temp_file_qf:
			# Looks if the unique ID is present on a line
			if "@SRR" in line:
				# adds to count 
				x += 1
		# Stores the total counts for the respective fastq sample
		temp_dict[fastq_name] = x
		# Closes the opened file
		temp_file_qf.close()

	# Returns the completed dictionary to the working environment
	return(temp_dict)


# Function to write out a summary file
def write_summary(full_dict, qf_dict, hrm_dict):
	# full_dict is the seq counts before any filtering
	# qf_dict is the seq counts after quality filtering
	# hrm_dict is the seq counts after human filtering removal 

	# Opens the necessary summary file
	test = open("%sseq_filter_summary.tsv" % (summarydir),'w')
	# Counter
	x = 0
	# Iterates through every sample 
	for sample_name in full_dict:
		# Checks if it is the first sample
		if x == 0:
			# Adds header and sample counts
			test.write("sample"+'\t'+"all_seqs"+'\t'+ \
				"quality_filtered_seqs"+'\t'+"human_filtered_seqs"+'\n'+ \
				sample_name+'\t'+ str(full_dict[sample_name])+'\t'+ \
				str(qf_dict[sample_name])+'\t'+str(hrm_dict[sample_name])+'\n')

		else:
			# Adds only sample counts
			test.write(sample_name+'\t'+ str(full_dict[sample_name])+'\t'+ \
				str(qf_dict[sample_name])+'\t'+str(hrm_dict[sample_name])+'\n')

		# adds to count
		x += 1
	# Clost the file
	test.close()


# Runs the overall program 
def main():

	meta_genome_file_name = command_line()
	samples_to_be_used = create_samples_to_download(meta_genome_file_name)
	test_full = get_seq_counts(samples_to_be_used, "_1.fastq")
	test_qf = get_seq_counts(samples_to_be_used, "_qf_1.fastq")
	test_hrm = get_seq_counts(samples_to_be_used, "_hrm_r1.fastq")

	write_summary(test_full, test_qf, test_hrm)

if __name__ == '__main__': main()