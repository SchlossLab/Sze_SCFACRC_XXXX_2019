#!bash

# Load needed modules
module load sratoolkit/2.8.2-1

#Set local variables
# mothurRv=/nfs/turbo/schloss-lab/bin/mothur_src/mothur
RAWDIR=data/raw/16S
rm -rf $RAWDIR
mkdir -p $RAWDIR

WORKDIR=data/mothur
rm -rf $WORKDIR
mkdir -p $WORKDIR

REF=data/references

#Get list of files to download
wget -O $RAWDIR/SRP062005_info.csv 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=SRP062005'

wget -O $RAWDIR/SRP096978_info.csv 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=SRP096978'

#Split sra files into fastq files
URLS=`cut -f 10 -d , $RAWDIR/*_info.csv | grep "http"`
for sample in $URLS
do
	echo $sample
	SRR=`echo $sample | sed -e "s_.*/__"`
	wget -O $RAWDIR/$SRR $sample
	fastq-dump --split-files $RAWDIR/$SRR -O $RAWDIR
	rm $RAWDIR/$SRR
done

# make files file
cut -d , -f 1,12 data/raw/16S/SRP0* | grep "SRR" | sed -E "s/(^SRR.*),([^_]*)_.*/\2\t\\1_1.fastq\t\\1_2.fastq/" > $RAWDIR/crc.files

#Run mothur process
mothur "#set.seed(seed=19760620);
	make.contigs(file=crc.files, inputdir=$RAWDIR, outputdir=$WORKDIR, processors=6);
	screen.seqs(fasta=current, group=current, maxambig=0, maxlength=275);
	unique.seqs(fasta=current);
	count.seqs(name=current, group=current);
	align.seqs(fasta=current, reference=$REF/silva.v4.align);
	screen.seqs(fasta=current, count=current, start=1968, end=11550, maxhomop=8);
	filter.seqs(fasta=current, vertical=T, trump=.);
	unique.seqs(fasta=current, count=current);
	pre.cluster(fasta=current, count=current, diffs=2);
	chimera.uchime(fasta=current, count=current, dereplicate=t);
	remove.seqs(fasta=current, accnos=current);
	classify.seqs(fasta=current, count=current, reference=$REF/trainset16_022016.pds.fasta, taxonomy=$REF/trainset16_022016.pds.tax, cutoff=80);
	remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);
	remove.groups(fasta=current, count=current, taxonomy=current, groups=mock1-mock2-mock5-mock6-mock7);
	cluster.split(fasta=current, count=current, taxonomy=current, taxlevel=4, cutoff=0.03);
	make.shared(list=current, count=current, label=0.03);
	sub.sample(shared=current);
	classify.otu(list=current, count=current, taxonomy=current, label=0.03);
	system(mv data/mothur/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons.taxonomy data/mothur/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons_rdp.taxonomy);
	classify.seqs(fasta=current, count=current, reference=$REF/gg_13_5_99.gg.fasta, taxonomy=$REF/gg_13_5_99.gg.tax, cutoff=80);
	classify.otu(list=current, count=current, taxonomy=current, label=0.03);
	system(mv data/mothur/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons.taxonomy data/mothur/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons_gg.taxonomy)"

cp data/mothur/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.shared data/mothur/crc.otu.shared
