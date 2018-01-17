#!python

# This is re-written python code to accomplish two tasks
	# Convert SRA to fastq
	# Run quality control over each read per sample


# Import needed libraries
import os

# Set up working directory
workdir = "data/raw/"

# Set up shortcuts for keys and other download components
asp = "/home/marcsze/.aspera/cli/etc/asperaweb_id_dsa.openssh"
sra_nih = "anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByStudy/sra"
database = "ERP"
seq_dir = "ERP005"
study = "ERP005534"

# Set up variables to be used in analysis
d_files = ["ERR478975", "ERR478971", "ERR478968", "ERR478984", 
"ERR478989", "ERR478976", "ERR478980", "ERR478964", "ERR478962", "ERR478960"]


# Cycles through and downloads each sample in d_files and moves it to appropriate directory
#for files in d_files:
	# Runs the actual download
#	os.system("ascp -QTr -k 1 -l 1G -i %s %s/%s/%s/%s/%s/ %s" % 
#		(asp, sra_nih, database, seq_dir, study, files, workdir))
	# Moves file back up a directory
#	os.system("mv %s%s/%s.sra %s" % (workdir, files, files, workdir))
	# Removes extra empty directory
#	os.system("rm -r %s%s" % (workdir, files))


# Convert the sra file to fastq files (one forward (F) and one reverse (R))
#for i in os.listdir(workdir):
	
#	if "sra" in i:

#		os.system("fastq-dump --split-files %s%s -O %s" % (workdir, i, workdir))


# Quality filter and toss out reads that don't make the cut 
# -t quality filter of nucleotide (trimmed from the end of sequence)
# -Q average quality cutoff of sequence
# -l minimum length of sequence if shorter after trimming it is tossed
for fastq_sample in d_files:
    
    os.system("sickle pe -f %s%s_1.fastq -r %s%s_2.fastq \
-t sanger -o %s%s_qf_1.fastq -p %s%s_qf_2.fastq \
-s %s%s_qf.orphan.fastq -q 30 -l 75" % 
    	(workdir, fastq_sample, workdir, fastq_sample, 
    		workdir, fastq_sample, workdir, fastq_sample, 
    		workdir, fastq_sample))
    
    print("Completed quality filtering of sample %s" % (fastq_sample))
