WORKDIR=data/picrust2
mkdir -p $WORKDIR

#dependencies
FASTA=data/asv/crc.asv.fasta
SHARED=data/asv/crc.asv.shared

BIOM=$WORKDIR/crc.asv.biom


mothur "#make.biom(shared=$SHARED, outputdir=$WORKDIR)"

mv $WORKDIR/crc.asv.ASV.biom $BIOM

# Comments on whether this is an acceptable biom format
biom validate-table -i $BIOM

# Convert to an acceptable biom format for picrust
# http://biom-format.org/documentation/biom_conversion.html
biom convert -i $BIOM -o $WORKDIR/temp_OTU_table.txt --table-type="OTU table" --to-tsv
biom convert -i $WORKDIR/temp_OTU_table.txt -o $BIOM.json --table-type="OTU table" --to-json
biom convert -i $WORKDIR/temp_OTU_table.txt -o $BIOM.hdf5 --table-type="OTU table" --to-hdf5

# Double check that it is now the right format
biom validate-table -i $BIOM.json
biom validate-table -i $BIOM.hdf5

rm $WORKDIR/temp_OTU_table.txt
