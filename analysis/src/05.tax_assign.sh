#!/usr/bin/bash
#DianaOaxaca
#Taxonomy assignment with sklearn and silva database in qiime2

##Train the Silva database taking V3-V4 amplicon to increase the taxonomic assignment resolution.

#mkdir -p data/dbs

#get the database and check its integrity

#wget https://data.qiime2.org/2023.5/common/silva-138-99-seqs.qza

#md5sum data/dbs/silva-138-99-seqs.qza 
#de8886bb2c059b1e8752255d271f3010  data/dbs/silva-138-99-seqs.qza

#wget https://data.qiime2.org/2023.5/common/silva-138-99-tax.qza

#md5sum data/dbs/silva-138-99-tax.qza 
#f12d5b78bf4b1519721fe52803581c3d  data/dbs/silva-138-99-tax.qza

#extract specific fragments
#qiime feature-classifier extract-reads \
#--i-sequences data/dbs/silva-138-99-seqs.qza \
#--p-f-primer CCTACGGGNGGCWGCAG --p-r-primer GACTACHVGGGTATCTAATCC \
#--p-min-length 250 --p-max-length 450 \
#--o-reads data/dbs/silva-138-99-seqs-extracted.qza --p-n-jobs 40

#train the database
#qiime feature-classifier fit-classifier-naive-bayes \
#--i-reference-reads data/dbs/silva-138-99-seqs-extracted.qza \
#--i-reference-taxonomy data/dbs/silva-138-99-tax.qza \
#--o-classifier data/dbs/classifier_silva_138_trained.qza

#taxonomy assignment
date
echo " Start taxonomy classification: "

qiime feature-classifier classify-sklearn \
  --i-classifier data/dbs/classifier_silva_138_trained.qza \
  --i-reads results/04.qiime/ASV_rep_seq.qza \
  --o-classification results/04.qiime/taxonomy.qza --p-n-jobs 40

date
echo " Get tax visualization "
#get visualization
qiime metadata tabulate \
  --m-input-file results/04.qiime/taxonomy.qza \
  --o-visualization results/04.qiime/taxonomy.qzv

echo " Get rep seqs visualization "
#get visual fasta to compare the taxonomic assignments with the top BLASTn hits for certain ASVs
qiime feature-table tabulate-seqs \
--i-data results/04.qiime/ASV_rep_seq.qza \
--o-visualization results/04.qiime/ASV_rep_seq.qzv

echo "Finish tax process"
date
