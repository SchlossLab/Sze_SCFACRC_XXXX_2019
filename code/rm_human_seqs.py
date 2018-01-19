#!python

# This code is used to remove human seqs from fastq files.
# It uses samtools to complete this task but any program that IDs human sequences can be used
# and it can be subsituted in that part of the code
# Workflow taken from http://www.metagenomics.wiki/tools/short-read/remove-host-sequences.
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
	# creates a list of all files in a directory
	file_dirs = os.listdir(refdir)
	# Original assumption is database is not in reference directory
	hg19_present = False
	# Set initial counter to 0
	x = 0
	# iterate through the directory list and only stop if the file end is reached
	while hg19_present == False & x <= len(file_dirs):

		# check for hg19 in directory
		if "hg19" in file_dirs[x]:
			# change to true if present
			hg19_present = True

		else:

			hg19_present = False
		# add one to the counter
		x += 1
	
	# check to see if hg19_present has changed from false to true
	if hg19_present == True:

		print("Human database present, skipping download step.")

	else:
		# Downloads the hg19 reference from bowtie
		os.system("wget -P %s ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip" % 
			(refdir))
		# Unzips the downloaded reference database in reference directory
		os.system("unzip %shg19.zip -d %s" % (refdir, refdir))


# Function to map and remove human sequences
def remove_human_seqs(samplesList):

	# Remove host sequences using bowtie2 and samtools
	for fastq_sample in samplesList:
		# Map sequences to human (hg19)
		print("Mapping reads to human geneome for %s." % (fastq_sample))

		os.system("bowtie2 -x hg19 -1 %s%s_qf_1.fastq -2 %s%s_qf_2.fastq -S %s%s_map_unmap.sam" % 
			(workdir, fastq_sample, workdir, fastq_sample, workdir, fastq_sample))

		# Convert files to bam
		print("Converting files to BAM format for %s." % (fastq_sample))
		os.system("samtools view -bS %s%s_map_unmap.sam > %s%s.map_unmap.bam" % (workdir, fastq_sample, workdir, fastq_sample))
		
		# Filter unmapped pairs
		print("Removing unmapped sequence pairs for %s." % (fastq_sample))
		os.system("samtools view -b -f 12 -F 256 %s%s.map_unmap.bam > %s%s.bothends_unmap.bam" % 
    		(workdir, fastq_sample, workdir, fastq_sample))
				# -f 12 = Extract only (-f) alignments with both reads unmapped: <read unmapped><mate unmapped>
        		# -F 256 = Do not(-F) extract alignments which are: <not primary alignment>

        # split paired-ends reads into separated fastq files		
		print("Splitting paired-end reads into F and R fastq for %s." % (fastq_sample))
		os.system("samtools sort -n %s%s.bothends_unmap.bam -o %s%s_bothends_unmap_sorted.bam" % 
    		(workdir, fastq_sample, workdir, fastq_sample))

		os.system("bedtools bamtofastq -i %s%s_bothends_unmap_sorted.bam -fq %s%s_hrm_r1.fastq -fq2 %s%s_hrm_r2.fastq" % 
    		(workdir, fastq_sample, workdir, fastq_sample, workdir, fastq_sample))

		# Remove intermediary files
		print("Removing uneeded intermediate files for %s." % (fastq_sample))
		os.system("rm %s%s.map_unmap.bam %s%s.bothends_unmap.bam %s%s_bothends_unmap_sorted.bam" % 
    		(workdir, fastq_sample, workdir, fastq_sample, workdir, fastq_sample))

    	
    	
# Runs the overall program 
def main():

	meta_genome_file_name = command_line()
	samples_to_be_used = create_samples_to_download(meta_genome_file_name)
	human_database_download()
	remove_human_seqs(samples_to_be_used)
	#print(samples_to_be_used)




if __name__ == '__main__': main()