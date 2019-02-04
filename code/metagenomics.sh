#!bash

module load sratoolkit/2.8.2-1
module load sickle/1.33.6
module load megahit/1.0.6
module load prodigal/2.6.3
module load diamond/0.8.34


#Set local variables
RAWDIR=data/raw/metagenome
mkdir -p $RAWDIR

WORKDIR=data/metagenome
mkdir -p $WORKDIR


#Get list of files to download
wget -O $RAWDIR/SRP108915_info.csv 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=SRP108915+AND+whole'

#Split sra files into fastq files
URLS=`cut -f 10 -d , $RAWDIR/SRP108915_info.csv | grep "http"`
for sample in $URLS
do
	fastq-dump --split-3 $sample -O $RAWDIR
	echo $sample
done

#Extract the list of SRR files
SRR_FILES=`cut -f 1 -d , $RAWDIR/SRP108915_info.csv | grep "SRR"`


#Quality filtering:	files taken from SRA so additional adapter cutting and removal of human reads
#										is not needed
for file in $SRR_FILES
do
	sickle pe -f $RAWDIR/${file}_1.fastq -r $RAWDIR/${file}_2.fastq -t sanger -o $WORKDIR/${file}_qf_1.fastq -p $WORKDIR/${file}_qf_2.fastq -s $WORKDIR/${file}_qf.orphan.fastq -q 30 -l 75
	echo ${file}
done


#Sequence assembly - t=12 threads
for file in $SRR_FILES
do
	megahit -1 $WORKDIR/${file}_qf_1.fastq -2 $WORKDIR/${file}_qf_2.fastq --min-contig-len 1000 --k-min 21 --k-max 101 --k-step 20 -t 12 -o $WORKDIR/${file}
	echo $file
done


#Create master contig file
> data/metagenome/all_contigs.fasta
for file in $SRR_FILES
do
sed "s/>/>${file}_/" $WORKDIR/${file}/final.contigs.fa >> $WORKDIR/all_contigs.fasta
done


# Marc had this, but the code did ">1000", not ">=1000". MegaHit, above already screens for length
#	python code/remove_short_contigs.py -s "data/raw/whole_metagenome_samples.txt"


#Generate table listing lengths of contigs
grep ">" $WORKDIR/all_contigs.fasta | sed -E "s/>(.*) flag.*len=/\1\t/" > $WORKDIR/all_contigs.lengths.tsv


#not sure this is actually ever used or why it is needed
#python2 code/concoct_scripts/cut_up_fasta.py -c 10000 -o 0 -m $(RAW)/all_contigs_greater1kb.fasta > $(RAW)/all_contigs_1kbto10kb.fasta


#python code/prodigal_wrap.py
# Call genes using prodigal
# the prodigal wiki is way ahead of its time and the options it lists aren't the same as what's
# available with version 2.6.3, which is still the latest release as of 2019-02-02
rm -rf $WORKDIR/prodigal
mkdir -p $WORKDIR/prodigal
split -l 1000 --suffix-length=5 $WORKDIR/all_contigs.fasta $WORKDIR/prodigal/
ls $WORKDIR/prodigal/* | xargs -I {} --max-procs=12 prodigal -q -c -i {} -o {}.genes -a {}.out -d {}.nucl -p meta
cat $WORKDIR/prodigal/*.nucl > $WORKDIR/all_contigs.genes.fasta
cat $WORKDIR/prodigal/*.out > $WORKDIR/all_contigs.aa.fasta
rm -rf $WORKDIR/prodigal



# Assign genes to OPFs using mmseq2
rm -rf $WORKDIR/mmseq2
mkdir -p $WORKDIR/mmseq2/temp
mkdir -p $WORKDIR/mmseq2/result

# De novo cluster genes into OPFs
## Set up target DB
mmseqs createdb $WORKDIR/all_contigs.aa.fasta $WORKDIR/mmseq2/target_DB

## Create index
mmseqs createindex $WORKDIR/mmseq2/target_DB $WORKDIR/mmseq2/tmp

## Cluster the target DB
mmseqs cluster $WORKDIR/mmseq2/target_DB $WORKDIR/mmseq2/target_clu $WORKDIR/mmseq2/tmp -e 0.001 --min-seq-id 0.4

## Create a target clustered sequence file
mmseqs createseqfiledb $WORKDIR/mmseq2/target_DB $WORKDIR/mmseq2/target_clu $WORKDIR/mmseq2/clu_seq

## Create target clustered fasta file
mmseqs result2flat $WORKDIR/mmseq2/target_DB $WORKDIR/mmseq2/target_DB $WORKDIR/mmseq2/clu_seq $WORKDIR/mmseq2/clu_seq.fasta

## Create a clustered target sequence tsv file
mmseqs createtsv $WORKDIR/mmseq2/target_DB $WORKDIR/mmseq2/target_DB $WORKDIR/mmseq2/target_clu $WORKDIR/mmseq2/clu.tsv #this could have a better name

#result2repseq ???


# Reference based-clustering
## Set up query DB
mmseqs createdb data/references/genes.pep.format.fasta $WORKDIR/mmseq2/query_DB

## Create alignment
mmseqs search $WORKDIR/mmseq2/query_DB $WORKDIR/mmseq2/target_DB $WORKDIR/mmseq2/resultDB $WORKDIR/mmseq2/tmp

## Create BLAST formatted file of Result DB (tsv based file)
mmseqs convertalis $WORKDIR/mmseq2/query_DB $WORKDIR/mmseq2/target_DB $WORKDIR/mmseq2/resultDB $WORKDIR/mmseq2/resultDB.m8

# cp $WORKDIR/mmseq2/resultDB.m8 data/process/orf_gene_alignment.tsv



rm -rf $WORKDIR/diamond
mkdir -p $WORKDIR/diamond

# the first grep is necessary because mmseq2 generates fasta headers without sequences and these
# need to be removed
grep -E "^>.*start|^[^>]" $WORKDIR/mmseq2/clu_seq.fasta > $WORKDIR/diamond/opf_proteins.fasta
diamond makedb --in $WORKDIR/diamond/opf_proteins.fasta -d $WORKDIR/diamond/reference_db

>$WORKDIR/diamond/orf_abund.tsv

for file in $SRR_FILES
do

#read 1
diamond blastx -q $WORKDIR/${file}_qf_1.fastq -d $WORKDIR/diamond/reference_db.dmnd -a $WORKDIR/diamond/${file}_1_output.daa -t $WORKDIR/diamond/ --max-target-seqs 1 --evalue 1e-15 --id 0.90 --threads 12

diamond view -a $WORKDIR/diamond/${file}_1_output.daa -o $WORKDIR/diamond/${file}_1.diamondout

#read 2
diamond blastx -q $WORKDIR/${file}_qf_2.fastq -d $WORKDIR/diamond/reference_db.dmnd -a $WORKDIR/diamond/${file}_2_output.daa -t $WORKDIR/diamond/ --max-target-seqs 1 --evalue 1e-15 --id 0.90 --threads 12

diamond view -a $WORKDIR/diamond/${file}_2_output.daa -o $WORKDIR/diamond/${file}_2.diamondout

#orphan read
diamond blastx -q $WORKDIR/${file}_qf.orphan.fastq -d $WORKDIR/diamond/reference_db.dmnd -a $WORKDIR/diamond/${file}_o_output.daa -t $WORKDIR/diamond/ --max-target-seqs 1 --evalue 1e-15 --id 0.90 --threads 12

diamond view -a $WORKDIR/diamond/${file}_o_output.daa -o $WORKDIR/diamond/${file}_o.diamondout

#pool results
cut -f 2 $WORKDIR/diamond/${file}_?.diamondout | sort | uniq -c | sed 's/^ *//' | sed "s/ /\\t/" | sed "s/$/\\t${file}/" >> $WORKDIR/diamond/orf_abund.tsv
done


#synthesize to make shared file
Rscript code/metagenomics_get_shared.R $WORKDIR/mmseq2/clu.tsv $WORKDIR/diamond/orf_abund.tsv
