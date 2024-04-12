#!/usr/bin/bash
#Get iqtree phylogeny

date
echo "Start phylogeny"

qiime phylogeny align-to-tree-mafft-iqtree \
 --p-n-threads auto --i-sequences results/04.qiime/ASV_rep_seq_filters.qza \
 --o-alignment results/04.qiime/align.qza \
 --o-masked-alignment results/04.qiime/masked-align.qza \
 --o-tree results/04.qiime/unrooted-tree-iqtree.qza \
 --o-rooted-tree results/04.qiime/rooted-tree-iqtree.qza --verbose

echo "finish phylogeny!"

date

