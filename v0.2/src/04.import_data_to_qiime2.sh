#!/usr/bin/bash
## DianaOaxaca
## Import data to QIIME2

#Run in qiime conda environment
#conda activate qiime2-2023.5

#import rep seqs
qiime tools import --input-path results/03.Dada2/ASVs.fasta --type 'FeatureData[Sequence]' --output-path results/04.qiime/ASV_rep_seq.qza

# append missing header to the table for import
cat <(echo -n "#OTU Table") results/03.Dada2/ASV_to_seqs-nochim.tsv > temp.txt

# convert to biom
biom convert -i temp.txt -o temp.biom --table-type="OTU table" --to-hdf5

# and create table-type qza
qiime tools import --input-path temp.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path results/04.qiime/ASV_table.qza

# remove temporal files
rm temp.*
