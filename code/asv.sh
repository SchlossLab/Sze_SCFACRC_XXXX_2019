WORKDIR=data/asv
mkdir -p $WORKDIR

#dependencies
FASTA=data/mothur/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta
COUNT=data/mothur/crc.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.pick.count_table
TAXONOMY=data/mothur/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.taxonomy

mothur "#make.shared(count=$COUNT, label=ASV, outputdir=$WORKDIR);rename.seqs(fasta=$FASTA, taxonomy=$TAXONOMY, map=$WORKDIR/crc.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.pick.map);sub.sample(shared=current)"

mv $WORKDIR/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.renamed.fasta $WORKDIR/crc.asv.fasta
mv $WORKDIR/crc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.renamed.taxonomy $WORKDIR/crc.asv.taxonomy
mv $WORKDIR/crc.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.pick.ASV.subsample.shared $WORKDIR/crc.asv.shared


rm $WORKDIR/crc*map
rm $WORKDIR/crc.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.pick.shared
