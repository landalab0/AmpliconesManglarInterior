#!/usr/bin/bash
#Accumulative curves

date

#Define package vectors
cran_packages <- c("knitr", "qtl", "bookdown", "magrittr", "plyr", "ggplot2",
                   "grid","gridExtra", "tidyverse", "devtools", "dplyr",
                   "pheatmap", "xtable",
                   "kableExtra", "remotes", "Rtsne", "vegan", "RColorBrewer",
                   "PoiClaClu",
                   "gtools", "gplots", "reshape2", "MASS", "usethis",
                   "indicspecies", "Polychrome")

bioc_packages <- c("airway", "phyloseq", "dada2", "DECIPHER", "phangorn",
                   "ggpubr","DESeq2",
                   "genefilter", "philr", "GenomeInfoDb", "microbiome",
                   "metagenomeSeq", "mia", "curatedMetagenomicData",
                   "ANCOMBC","microbiomeMarker")

git_packages <- c("btools", "fantaxtic", "ampvis2", "tsnemicrobiota",
                  "qiime2R", "ranacapa")
                  

#Load libraries
sapply(c(cran_packages, bioc_packages, git_packages), require, character.only = TRUE)

# 01. Load data
physeq_qiime2 <- qza_to_phyloseq(
  features = "../results/04.qiime/ASV_table_filter_freq218_aemc.qza",
  #tree = "results/04.qiime/",
  taxonomy = "../results/04.qiime/taxonomy.qza",
  metadata = "../data/metadata.tsv")

#Acumulation curves
#library(ranacapa)
## 03.2 Get accumulation curves
acumulation_curves <- ggrare(physeq_qiime2, step = 100, color = "Site", label = "Sample")

#custom plot
acumulation_curves_plot <- acumulation_curves + facet_wrap(~Site) +
  labs(title="Accumulative curves") + theme_bw() + scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")

#save plot
pdf("../results/plots/Rarefaction_curves_ranacapa.pdf")
acumulation_curves_plot
dev.off()

acumulation_curves_plot
ggsave("../results/plots/Rarefaction_curves_ranacapa.png")

date
