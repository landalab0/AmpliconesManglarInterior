# Diversity

## Explore data

#load packages

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
                   "metagenomeSeq", "mia",
                   "ANCOMBC","microbiomeMarker")

git_packages <- c("btools", "fantaxtic", "ampvis2", "tsnemicrobiota",
                  "qiime2R", "ranacapa")

#Load libraries
sapply(c(cran_packages, bioc_packages, git_packages), require, character.only = TRUE)

#Import data to R from qiime and check prevalence to filter data 

# 01. Load data
physeq_qiime2 <- qza_to_phyloseq(
  features = "results/04.qiime/ASV_table_filter_aemc_freq94_1minsamp.qza",
  tree = "results/04.qiime/rooted-tree-fasttree2.qza",
  taxonomy = "results/04.qiime/taxonomy.qza",
  metadata = "data/metadata.tsv")

physeq_qiime2

# 02. Explore prevalence
## 02.1 Get prevalence
prevdf = apply(X = otu_table(physeq_qiime2),
               MARGIN = ifelse(taxa_are_rows(physeq_qiime2), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
## 02.2 Add taxonomy
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(physeq_qiime2),
                    tax_table(physeq_qiime2))
## 02.3 Check prevalence at Phylum level
dfprev <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

## 02.4 Get genus prevalence plot
prevalence_genus = subset(prevdf, Genus %in% get_taxa_unique(physeq_qiime2, "Genus"))
prev_genus <- ggplot(prevalence_genus, aes(TotalAbundance, 
                                           Prevalence /nsamples(physeq_qiime2),color=Genus)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#save prevalence plot
pdf("results/plots/03.diversity/01.Prevalence_Phylum_and_genus.pdf")
prev_genus
dev.off()

## 3.1 FASTER rarefaction curves with Vegan
mat <- as(t(otu_table(physeq_qiime2)), "matrix")
raremax <- min(rowSums(mat))
pdf("results/plots/03.diversity/02.Rarefaction_curves_vegan.pdf")
vegan_rarefaction_curves <- system.time(rarecurve(mat, step = 100, sample = raremax, 
                                                  col = "blue", label = FALSE))
dev.off()

vegan_rarefaction_curves

#Acumulation curves
#library(ranacapa)
## 03.2 Get accumulation curves
acumulation_curves <- ggrare(physeq_qiime2, step = 100, color = "Location", label = "Sample")

acumulation_curves_plot <- acumulation_curves + facet_wrap(~Location) +
  labs(title="Accumulative curves") + theme_bw()

pdf("results/plots/03.diversity/03.Rarefaction_curves_ranacapa.pdf")
acumulation_curves_plot
dev.off()

##png for summary
png("results/plots/summary/06.acumulation_curves.png")
acumulation_curves_plot
dev.off()

print(acumulation_curves_plot)


################################################################################
# Bray NMDS 
nmds_bray <- ordinate(physeq_qiime2, method = "NMDS", distance = "bray")
# Get stress value
var_stress_nmds_bray <- round(nmds_bray$stress, 5)

#checks that the fit is good with shepard plot
stressplot(nmds_bray)

nmds_bray_plot <- plot_ordination(physeq_qiime2, nmds_bray, label = "Sample",
                           color = "Location", shape = "Location") + theme_bw() + 
  labs(col = "Location") + labs(title="NMDS, Bray-Curtis distance") +
  geom_point(size=3) + theme_bw() 

nmds_bray_plot <- nmds_bray_plot +
  annotate("text", x = Inf, y = -Inf, label = paste("Stress:", var_stress_nmds_bray),
           hjust = 2.1, vjust = -1.9, size = 4)
nmds_bray_plot

pdf("results/plots/03.diversity/04.NMDS_Shepard_Bray_Fit.pdf")
stressplot(nmds_bray)
dev.off()

pdf("results/plots/03.diversity/05.NMDS_Bray_Location.pdf")
nmds_bray_plot
dev.off()

##png for summary
png("results/plots/summary/07.nmds_bray_location.png")
nmds_bray_plot
dev.off()


#Weigthed UniFrac take relative abundance and it is less sensitive to sample size 

nmds_wunifrac <- ordinate(physeq_qiime2, method = "NMDS", distance = "wunifrac")
# stress variable
var_stress_nmds_wu <- round(nmds_wunifrac$stress, 5)
var_stress_nmds_wu

stressplot(nmds_wunifrac)# checks that the fit is good
#nmds_wunifrac$points

# Weigthed UniFrac NMDS
nmds_wu <- plot_ordination(physeq_qiime2, nmds_wunifrac, label = "Sample",
                           color = "Location", shape = "Location") + theme_bw() + 
  labs(col = "Location") + labs(title="NMDS, Weighted UniFrac distance") +
  geom_point(size=3) + theme_bw() 

nmds_wu <- nmds_wu +
  annotate("text", x = Inf, y = -Inf, label = paste("Stress:", var_stress_nmds_wu),
           hjust = 1.1, vjust = -1.1, size = 4)
nmds_wu

pdf("results/plots/03.diversity/06.NMDS_Shepard_WUniFrac_Fit.pdf")
stressplot(nmds_wunifrac)
dev.off()

pdf("results/plots/03.diversity/07.NMDS_WUniFrac_Location.pdf")
nmds_wu
dev.off()

#### ANOSIM
#extract metadata
metadata <- data.frame(phyloseq::sample_data(physeq_qiime2),
           check.names = FALSE)
#get significant difference between location with anosim
anosim_location <- anosim(x= as.data.frame(t(otu_table(physeq_qiime2))),
       grouping = metadata$Location,
       permutations = 9999, distance = "bray")
#get values
anosim_significance <- anosim_location$signif
anosim_statistic  <- anosim_location$statistic
anosim_location
####
