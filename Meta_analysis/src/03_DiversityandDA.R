physeq_qiime3 <- qza_to_phyloseq(
  features = "cluster_table_filter_freq218_emcEstero.qza",
  #tree = "results/04.qiime/",
  taxonomy = "taxonomyEstero.qza",
  metadata = "metadata14.tsv")

#Filter metadata
# Extract metadata from phyloseq object
metadat <- data.frame(phyloseq::sample_data(physeq_qiime3), check.names = FALSE)

#Filter metadata by Study_zone column
filtered_metadat <- metadat[metadat$Study_zone %in% c("Estero_pargo", "Celestún", "Rio San Pedro"), ]

#Filter again by Ecological Type column, selecting Fringe parameter
filtered_metadat2 <- filtered_metadat[filtered_metadat$Ecological_type %in% c("Fringe"), ]

#Filter again by season column to include all samples that were collected form a flood season
filtered_metadat3 <- filtered_metadat2[filtered_metadat2$season %in%c("flood"),]

#Combine filtered metadata with phyloseq object to create a filtered phyloseq object
filtered_physeq2 <- prune_samples(rownames(filtered_metadat), physeq_qiime3)

#Calculate Bray-Curtis distance for filtered phyloseq
nmds_bray_filtered1 <- ordinate(filtered_physeq2, method = "NMDS", distance = "bray")

#Graph results in a NMDS plot with Study_zone and depth as parameters
nmds_bray_filtered_plot <- plot_ordination(filtered_physeq2, nmds_bray_filtered1,
                                           label = "SRA", color = "Study_zone") +
  theme_bw() +
  geom_point(aes(shape = depth), size = 3) +  # Agrega forma según profundidad
  labs(title = "NMDS, Bray-Curtis distance") +
  scale_color_npg() +
  scale_fill_npg() +
  scale_shape_manual(values = c(16, 17, 18, 15, 8, 3, 6, 2))  # Personaliza formas geométricas

# Show plot
nmds_bray_filtered_plot

#                   ANCOMBC2
################################################
# Load package
library(qiime2R)
library(ANCOMBC)
library(phyloseq)
library(tidyverse)
library(DT)
library(ggplot2)
library(vegan)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color':
  '#000', 'color': '#fff'});","}")))


# phyloseq to TreeSummarizedExperiment
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(filtered_physeq2)

# ANCOMBC2
output <- ancombc2(
  data = tse,
  assay_name = "counts",
  tax_level = "Phylum",
  fix_formula = "Study_zone",
  #rand_formula = "(1 | ID)",
  p_adj_method = "BH",
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 1000,
  s0_perc = 0.05,
  group = "Study_zone",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 2,
  verbose = TRUE,
  global = TRUE,
  pairwise = TRUE,
  dunnet = TRUE,
  trend = TRUE,
  iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
  em_control = list(tol = 1e-5, max_iter = 100),
  lme_control = lme4::lmerControl(),
  mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                              nrow = 2,
                                              byrow = TRUE),
                                       matrix(c(-1, 0, 1, -1),
                                              nrow = 2,
                                              byrow = TRUE),
                                       matrix(c(1, 0, 1, -1),
                                              nrow = 2,
                                              byrow = TRUE)),
                       node = list(2, 2, 1),
                       solver = "ECOS",
                       B = 100
  )
)

##################################################
#          LOCATION
##################################################
#To manually change the reference level, not alphabdetic order
tse$Location = factor(tse$Location, levels = c("El Cacahuate", "La Piedad", "Dique Miguelito"))


#variable
tse$Location = recode(as.character(tse$Location),
                      `El Cacahuate` = "EC",
                      `La Piedad` = "LP",
                      `Dique Miguelito` = "DM",
                      .missing = "unknown")

# Ejecutar ANCOMBC2 para comparar por Location
output_location <- ancombc2(
  data = tse,
  assay_name = "counts",
  tax_level = "Phylum",
  fix_formula = "Depth_cm + Location",  # Cambiamos la fórmula fija
  rand_formula = "(1 | ID)",
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 1000,
  s0_perc = 0.05,
  group = "Location",  # Cambiamos el grupo de interés
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 2,
  verbose = TRUE,
  global = TRUE,
  pairwise = TRUE,
  dunnet = TRUE,
  trend = TRUE,
  iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
  em_control = list(tol = 1e-5, max_iter = 100),
  lme_control = lme4::lmerControl(),
  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                              nrow = 2,
                                              byrow = TRUE),
                                       matrix(c(-1, 0, 1, -1),
                                              nrow = 2,
                                              byrow = TRUE),
                                       matrix(c(1, 0, 1, -1),
                                              nrow = 2,
                                              byrow = TRUE)),
                       node = list(2, 2, 1),
                       solver = "ECOS",
                       B = 100))

# Ver los resultados
res_pair_location <- output$res_pair

# Visualizar los resultados
head(res_pair_location)

# Crear el dataframe df_fig_pair1 para los valores log-fold change (LFC)
df_fig_pair1_location <- res_pair_location %>%
  dplyr::filter(`diff_Study_zoneLaguna Cacahuate` == 1 |
                  `diff_Study_zoneRio San Pedro` == 1 |
                  `diff_Study_zoneRio San Pedro_Study_zoneLaguna Cacahuate` == 1) %>%
  dplyr::mutate(lfc1 = ifelse(`diff_Study_zoneLaguna Cacahuate` == 1, round(`lfc_Study_zoneLaguna Cacahuate`, 2), NA),
                lfc2 = ifelse(`diff_Study_zoneRio San Pedro` == 1, round(`lfc_Study_zoneRio San Pedro`, 2), NA),
                lfc3 = ifelse(`diff_Study_zoneRio San Pedro_Study_zoneLaguna Cacahuate` == 1, round(`lfc_Study_zoneRio San Pedro_Study_zoneLaguna Cacahuate`, 2), NA)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3, names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

# Crear el dataframe df_fig_pair2 para los valores significativos
df_fig_pair2_location <- res_pair_location %>%
  dplyr::filter(`diff_Study_zoneLaguna Cacahuate` == 1 |
                  `diff_Study_zoneRio San Pedro` == 1 |
                  `diff_Study_zoneRio San Pedro_Study_zoneLaguna Cacahuate` == 1) %>%
  dplyr::mutate(lfc1 = ifelse(`passed_ss_Study_zoneLaguna Cacahuate` == 1 & `diff_Study_zoneLaguna Cacahuate` == 1, "#00CCFF", "black"),
                lfc2 = ifelse(`passed_ss_Study_zoneRio San Pedro` == 1 & `diff_Study_zoneRio San Pedro` == 1, "#00CCFF", "black"),
                lfc3 = ifelse(`passed_ss_Study_zoneRio San Pedro_Study_zoneLaguna Cacahuate` == 1 & `diff_Study_zoneRio San Pedro_Study_zoneLaguna Cacahuate` == 1, "#00CCFF", "black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3, names_to = "group", values_to = "color") %>%
  dplyr::arrange(taxon)

# Combinar los dataframes df_fig_pair1 y df_fig_pair2
df_fig_pair_location <- df_fig_pair1_location %>%
  dplyr::left_join(df_fig_pair2_location, by = c("taxon", "group"))

# Renombrar los grupos para hacerlos más legibles
df_fig_pair_location$group <- recode(df_fig_pair_location$group,
                                     `lfc1` = "LC vs E",
                                     `lfc2` = "RSP vs E",
                                     `lfc3` = "RSP vs LC")

# #df_fig_pair_location$group <- factor(df_fig_pair_location$group,
#                                      levels = c("LC vs ",
#                                                 "LP vs DM",
#                                                 "LP vs EC"))

# Filtrar las comparaciones que no tienen valores significativos
df_fig_pair_location <- df_fig_pair_location %>%
  group_by(group) %>%
  filter(any(!is.na(value))) %>%
  ungroup()

# Calcular los límites y el punto medio para la escala del heatmap
lo <- floor(min(df_fig_pair_location$value, na.rm = TRUE))
up <- ceiling(max(df_fig_pair_location$value, na.rm = TRUE))
mid <- (lo + up) / 2

# Crear el heatmap
#fig_pair_location <- df_fig_pair_location %>%
  ggplot(aes(x = group, y = taxon, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#4cbfb7", high = "#c866a3", mid = "gray90",#c866a3
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = "LFC") +
  geom_text(aes(label = round(value, 2), color = color), size = 5, na.rm = TRUE) +
  scale_color_identity(guide = FALSE) +
  labs(x = "Location Comparison", y = NULL, title = "LogFold Change between locations") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9))

fig_pair_location <- df_fig_pair_location %>%
    ggplot(aes(x = group, y = taxon, fill = value)) +
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "#4cbfb7", high = "#c866a3", mid = "gray90",
                         na.value = "white", midpoint = mid, limit = c(lo, up),
                         name = "LFC") +
    geom_text(aes(label = round(value, 2)), color = "black", size = 5, na.rm = TRUE) +  # Números en negro
    labs(x = "Location Comparison", y = NULL, title = "LogFold Change between locations") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 9))

# Mostrar el heatmap
print(fig_pair_location)
