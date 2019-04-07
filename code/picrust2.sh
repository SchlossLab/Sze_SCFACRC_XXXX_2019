THREADS=6

module load python-anaconda3/latest
module load R/3.5.1

WORKDIR=data/picrust2
rm -rf $WORKDIR
mkdir -p $WORKDIR

#dependencies
FASTA=$WORKDIR/crc.asv.fasta
BIOM=$WORKDIR/crc.asv.biom

R -e "source('code/picrust2_utilities.R'); sort_and_unalign_sequences()"

mothur "#make.biom(shared=data/asv/crc.asv.shared, outputdir=$WORKDIR)"

mv $WORKDIR/crc.asv.ASV.biom $BIOM

conda activate
conda activate picrust2

# Comments on whether this is an acceptable biom format
biom validate-table -i $BIOM

# Convert to an acceptable biom format for picrust
# http://biom-format.org/documentation/biom_conversion.html
biom convert -i $BIOM -o $WORKDIR/temp_OTU_table.txt --table-type="OTU table" --to-tsv
biom convert -i $WORKDIR/temp_OTU_table.txt -o $BIOM --table-type="OTU table" --to-json #--to-hdf5

# Double check that it is now the correct format
biom validate-table -i $BIOM

rm $WORKDIR/temp_OTU_table.txt

picrust2_pipeline.py -s $FASTA -i $BIOM -o $WORKDIR/output --threads $THREADS --verbose

conda deactivate

PC_OUTPUT=$WORKDIR/output

R -e "source('code/picrust2_utilities.R');convert_tsv_to_shared('$PC_OUTPUT/pathways_out/path_abun_unstrat.tsv')"
R -e "source('code/picrust2_utilities.R');convert_tsv_to_shared('$PC_OUTPUT/KO_metagenome_out/pred_metagenome_unstrat.tsv')"
R -e "source('code/picrust2_utilities.R');convert_tsv_to_shared('$PC_OUTPUT/EC_metagenome_out/pred_metagenome_unstrat.tsv')"
