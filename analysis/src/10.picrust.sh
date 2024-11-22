# picrust
nohup qiime picrust2 full-pipeline --i-table results/04.qiime/ASV_table_filter_freq218_emc.qza --i-seq results/04.qiime/ASV_rep_seq_filters.qza --output-dir results/04.qiime/picrust2 --p-threads 40 --verbose > outs/picrust2.nohup &

# visual
qiime feature-table summarize --i-table results/04.qiime/picrust2/pathway_abundance.qza --o-visualization results/04.qiime/picrust2/pathway_abundance.qzv
qiime feature-table summarize --i-table results/04.qiime/picrust2/ec_metagenome.qza --o-visualization results/04.qiime/picrust2/ec_metagenome.qzv
qiime feature-table summarize --i-table results/04.qiime/picrust2/pathway_abundance.qza --o-visualization results/04.qiime/picrust2/pathway_abundance.qzv

# export
qiime tools export --input-path results/04.qiime/picrust2/ko_metagenome.qza --output-path results/04.qiime/picrust2/ko_metagenome-exported

# ancom 
qiime composition ancombc --i-table results/04.qiime/picrust2/pathway_abundance.qza --m-metadata-file data/metadata3.tsv --p-formula Depths --o-differentials results/04.qiime/picrust2/ancombc_DA_pathways.qza

qiime composition da-barplot --i-data results/04.qiime/picrust2/ancombc_DA_pathways.qza --o-visualization results/04.qiime/picrust2/ancombc_DA_pathways.qzv
qiime tools export --input-path results/04.qiime/picrust2/ancombc_DA_pathways.qza --output-path results/04.qiime/picrust2/ancombc_DA_pathways-exported

qiime tools export --input-path results/04.qiime/picrust2/pathway_abundance.qza --output-path results/04.qiime/picrust2/pathway_abundance-exported

biom convert -i results/04.qiime/picrust2/pathway_abundance-exported/feature-table.biom -o results/04.qiime/picrust2/pathway_abundance-exported/pathway_abundance.txt --to-tsv

Rscript pathways_descriptions.R -i results/04.qiime/picrust2/pathway_abundance-exported/pathway_abundance.txt -o results/04.qiime/picrust2/pathway_abundance-exported/pathway_abundance_descriptions.txt
