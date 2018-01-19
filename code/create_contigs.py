#!python

# This code is used to create contigs for each sample.
# It parts are taken from the perl code of Geoffrey Hannigan and 
# they can be found here if needed 
#(https://github.com/SchlossLab/Hannigan_CRCVirome_mBio_2017/blob/master/bin/ContigAssembly.sh)
# I use megahit like Geoff but feel free to sub-in your algorithm of choice


############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys

# Import other code with useful functions that will be used
from qual_trim import command_line, create_samples_to_download

# Set up working directory
workdir = "data/raw/"
refdir = "data/references/"


############################################################################################

# Function to run the assemble contig command of megahit
def make_contigs(sampleList):

	for i in sampleList:

		os.system("megahit -1 %s%s_hrm_r1.fastq -2 %s%s_hrm_r2.fastq \
			--min-contig-len 1000 --k-min 21 --k-max 101 --k-step 20 -t 4 -o %s%s" % 
		(workdir, i, workdir, i, workdir, i))

	








# Runs the overall program 
def main():

	meta_genome_file_name = command_line()
	samples_to_be_used = create_samples_to_download(meta_genome_file_name)
	make_contigs(samples_to_be_used)




if __name__ == '__main__': main()


