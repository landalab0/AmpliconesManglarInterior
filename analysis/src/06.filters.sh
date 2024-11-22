#!/usr/bin/bash
#DianaOaxaca
#Filters

#Summary of the qza table imported from R
qiime feature-table summarize \
--i-table results/04.qiime/ASV_table.qza \
--o-visualization results/04.qiime/ASV_table.qzv

#Filter by frequency
#Here I removed all ASVs with a frequency of less than 0.1% of the mean sample depth. 
#This cut-off excludes ASVs that are likely due to MiSeq bleed-through between runs (reported by Illumina to be 0.1% of reads). 
#To calculate this cut-off I identified the mean sample depth, multiplied it by 0.001, and rounded to the nearest integer. 
#This step are describe in [this paper](https://journals.asm.org/doi/pdf/10.1128/msystems.00127-16)

qiime feature-table filter-features --i-table  results/04.qiime/ASV_table.qza \
 --p-min-samples 1 --p-min-frequency 218 --o-filtered-table results/04.qiime/ASV_table_filter_freq218.qza

qiime feature-table summarize --i-table results/04.qiime/ASV_table_filter_freq218.qza \
 --o-visualization results/04.qiime/ASV_table_filter_freq218.qzv

#Filter Mitochondria, chloroplast and Eukaryota

qiime taxa filter-table --i-table results/04.qiime/ASV_table_filter_freq218.qza \
 --i-taxonomy results/04.qiime/taxonomy.qza --p-exclude Eukaryota,Mitochondria,Chloroplast \
 --p-include p__ --o-filtered-table results/04.qiime/ASV_table_filter_freq218_emc.qza

qiime feature-table summarize --i-table results/04.qiime/ASV_table_filter_freq218_emc.qza \
 --o-visualization results/04.qiime/ASV_table_filter_freq218_emc.qzv

#remove in fasta sequences
qiime feature-table filter-seqs  --i-table results/04.qiime/ASV_table_filter_freq218_emc.qza \
 --i-data results/04.qiime/ASV_rep_seq.qza --o-filtered-data results/04.qiime/ASV_rep_seq_filters.qza
