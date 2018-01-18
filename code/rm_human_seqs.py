#!python

# This code is used to remove human seqs from fastq files.
# It uses samtools to complete this task but any program that IDs human sequences can be used
# and it can be subsituted in that part of the code.
# The code also used hg19 but if you really want to you can use the newer hg23 instead.

############## Internal parameters used by all functions is program ########################

# Import needed libraries
import os, sys

# Import other code with useful functions that will be used
from qual_trim import command_line, create_samples_to_download

# Set up working directory
workdir = "data/raw/"
refdir = "data/references/"


############################################################################################

# Function used to download human reference genome or to skip the step if present already
def human_database_download():

	file_dirs = os.listdir(refdir)

	hg19_present = False

	x = 0

	while hg19_present == False:

		if "hg19" in file_dirs[x]:

			hg19_present = True

		else:

			hg19_present = False

		x += 1

	if hg19_present == True:

		print("Human database present, skipping download step.")

	else:

		os.system("wget -P %s ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip" % 
			(refdir))

		os.system("unzip %shg19.zip -d %s" % (refdir, refdir))


def main():

	meta_genome_file_name = command_line()
	samples_to_be_used = create_samples_to_download(meta_genome_file_name)
	human_database_download()

	#print(samples_to_be_used)




if __name__ == '__main__': main()