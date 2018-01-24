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

	temp_dict = {}

	for fastq_name in sampleList:

		print("Counting sequences in %s%s" % (fastq_name, seq_ending))
		
		x = 0

		temp_file_qf = open("%s%s%s" % 
			(workdir, fastq_name, seq_ending), 'r')

		for line in temp_file_qf:

			if "@SRR" in line:

				x += 1

		temp_dict[fastq_name] = x

		temp_file_qf.close()


	return(temp_dict)


# Function to write out a summary file
def write_summary(full_dict, qf_dict, hrm_dict):

	test = open("%sseq_filter_summary.tsv" % (summarydir),'w')

	x = 0

	for sample_name in full_dict:
		
		if x == 0:

			test.write("sample"+'\t'+"all_seqs"+'\t'+ \
				"quality_filtered_seqs"+'\t'+"human_filtered_seqs"+'\n'+ \
				sample_name+'\t'+ str(full_dict[sample_name])+'\t'+ \
				str(qf_dict[sample_name])+'\t'+str(hrm_dict[sample_name])+'\n')

		else:

			test.write(sample_name+'\t'+ str(full_dict[sample_name])+'\t'+ \
				str(qf_dict[sample_name])+'\t'+str(hrm_dict[sample_name])+'\n')

		
		x += 1

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