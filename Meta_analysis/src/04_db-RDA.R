
library(phyloseq)
library(qiime2R)
library(tydiverse)
library(dplyr)

#PCoA with environmental parameters
#Extract data with environmental parameters associated
metadatatonga <- data.frame(phyloseq::sample_data(physeq_qiime3), check.names = FALSE)

#Filter data by the Parametros column with the value "si"
PCoA_metadatatonga <- metadatatonga[metadatatonga$Parametros %in% c("si"), ]

#Combine phyloseq object with the filtered metadata
PCoA_physeqonga <- prune_samples(rownames(PCoA_metadatatonga), physeq_qiime2)

# Extract metadata
PCoA_metadata <- as.data.frame(sample_data(PCoA_physeqonga))

# select environmental parameters you want to use
# for example, "pH", "Temperature", "Salinity", etc.
env_parameters <- PCoA_metadata[, c("Temperatura", "Salinidad", "pH", "Redox.mV", "S.2", "SO4")]

# Verify values with NA
env_parameters <- na.omit(env_parameters)

#calculate Bray Curtis distance
bray_dist <- phyloseq::distance(PCoA_physeqonga, method = "bray")

# make an db-RDA (using capscale in vegan)
dbrda_result <- capscale(bray_dist ~ Temperatura + Salinidad + pH + Redox.mV + S.2 + SO4, data = PCoA_metadatatonga)

# see results
summary(dbrda_result)

# Graph results of db-RDA
plot(dbrda_result, scaling = 2, main = "db-RDA: Community vs Environmental Parameters")

# add arrows that represents environmental parameters
envfit_result <- envfit(dbrda_result, env_parameters, perm = 999)
plot(envfit_result, col = "red")

# Permutation test
anova(dbrda_result, permutations = 999)

#Significance of every axis
anova(dbrda_result, by = "terms")

#Kruskal-Wallis test to see singificance in no parametric data
# for every parameter in env_parameters, apply Kruskal-Wallis between Study_zone
kruskal_results <- lapply(env_parameters, function(param) {
  kruskal.test(param ~ PCoA_metadatatonga$Study_zone)
})

# show results
kruskal_results
