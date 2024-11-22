#!/usr/bin/bash
#DianaOaxaca
# Get info per sample


#taxonomy
qiime tools export --input-path results/04.qiime/taxonomy.qza --output-path results/04.qiime/exports/

#feature table
qiime tools export --input-path results/04.qiime/ASV_table_filter_freq218_emc.qza --output-path results/04.qiime/exports

#reformat taxonomy tsv
sed -i -e '1 s/Feature/#Feature/' -e '1 s/Taxon/taxonomy/' results/04.qiime/exports/taxonomy.tsv

#Add taxonomy to feature table
biom add-metadata -i results/04.qiime/exports/feature-table.biom -o results/04.qiime/exports/feature-table_tax.biom --observation-metadata-fp results/04.qiime/exports/taxonomy.tsv --sc-separated results/04.qiime/exports/taxonomy.tsv

#Convert to tsv from biom format
biom convert -i results/04.qiime/exports/feature-table_tax.biom -o results/04.qiime/exports/feature-table_tax.tsv --to-tsv --header-key taxonomy

#Get effective reads per sample
biom summarize-table -i results/04.qiime/exports/feature-table_tax.biom -o results/04.qiime/exports/feature-table_tax_reads_per_sample_summary.txt

#Get ASVs per sample
biom summarize-table -i results/04.qiime/exports/feature-table_tax.biom --qualitative -o results/04.qiime/exports/feature-table_tax_ASVs_per_sample_summary.txt

