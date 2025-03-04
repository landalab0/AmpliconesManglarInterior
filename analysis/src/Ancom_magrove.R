################################################
#                   ANCOMBC2
################################################
# Load package
library(qiime2R)
library(ANCOMBC)
library(phyloseq)
library(tidyverse)
library(DT)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color':
  '#000', 'color': '#fff'});","}")))

#load("~/Documents/Manglares/update/src/ancom_mangrove.RData")
# phyloseq object
ps <- qza_to_phyloseq(
  features = "ASV_table_filter_freq218_emc.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv")

# phyloseq to TreeSummarizedExperiment
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(ps)

# contrast matrix
contrast_matrices <- list(
  matrix(c(1, -1, 0,
           0, 1, -1,
           0, 0, 1),
         nrow = 3, byrow = TRUE),
  matrix(c(1, 0, -1,
           0, 1, -1,
           0, 0, 1),
         nrow = 3, byrow = TRUE),
  matrix(c(1, 0, 0,
           0, 1, -1,
           0, 0, 1),
         nrow = 3, byrow = TRUE),
  matrix(c(1, -1, 0,
           0, 1, 0,
           0, 0, 1),
         nrow = 3, byrow = TRUE))

# nodes
nodes <- list(1, 1, 1, 1)

# ANCOMBC2
output <- ancombc2(
  data = tse,
  assay_name = "counts",
  tax_level = "Phylum",
  fix_formula = "Location + Depth_cm",
  rand_formula = "(1 | ID)",
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 1000,
  s0_perc = 0.05,
  group = "Depth_cm",
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
  trend_control = list(
    contrast = contrast_matrices,
    node = nodes,
    solver = "ECOS",
    B = 100
  )
)

library(dplyr)
library(tidyr)
library(ggplot2)

# Extraer los resultados de comparaciones pareadas del objeto output
res_pair <- output$res_pair

# Crear el dataframe df_fig_pair1
df_fig_pair1 <- res_pair %>%
  dplyr::filter(`diff_Depth_cm16-30` == 1 |
                  `diff_Depth_cm31-45` == 1 |
                  `diff_Depth_cm50-75` == 1 |
                  `diff_Depth_cm31-45_Depth_cm16-30` == 1 |
                  `diff_Depth_cm50-75_Depth_cm16-30` == 1 |
                  `diff_Depth_cm50-75_Depth_cm31-45` == 1) %>%
  dplyr::mutate(lfc1 = ifelse(`diff_Depth_cm16-30` == 1, round(`lfc_Depth_cm16-30`, 2), NA),
                lfc2 = ifelse(`diff_Depth_cm31-45` == 1, round(`lfc_Depth_cm31-45`, 2), NA),
                lfc3 = ifelse(`diff_Depth_cm50-75` == 1, round(`lfc_Depth_cm50-75`, 2), NA),
                lfc4 = ifelse(`diff_Depth_cm31-45_Depth_cm16-30` == 1, round(`lfc_Depth_cm31-45_Depth_cm16-30`, 2), NA),
                lfc5 = ifelse(`diff_Depth_cm50-75_Depth_cm16-30` == 1, round(`lfc_Depth_cm50-75_Depth_cm16-30`, 2), NA),
                lfc6 = ifelse(`diff_Depth_cm50-75_Depth_cm31-45` == 1, round(`lfc_Depth_cm50-75_Depth_cm31-45`, 2), NA)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc6, names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

# Crear el dataframe df_fig_pair2
df_fig_pair2 <- res_pair %>%
  dplyr::filter(`diff_Depth_cm16-30` == 1 |
                  `diff_Depth_cm31-45` == 1 |
                  `diff_Depth_cm50-75` == 1 |
                  `diff_Depth_cm31-45_Depth_cm16-30` == 1 |
                  `diff_Depth_cm50-75_Depth_cm16-30` == 1 |
                  `diff_Depth_cm50-75_Depth_cm31-45` == 1) %>%
  dplyr::mutate(lfc1 = ifelse(`passed_ss_Depth_cm16-30` == 1 & `diff_Depth_cm16-30` == 1, "blue3", "black"),
                lfc2 = ifelse(`passed_ss_Depth_cm31-45` == 1 & `diff_Depth_cm31-45` == 1, "blue3", "black"),
                lfc3 = ifelse(`passed_ss_Depth_cm50-75` == 1 & `diff_Depth_cm50-75` == 1, "blue3", "black"),
                lfc4 = ifelse(`passed_ss_Depth_cm31-45_Depth_cm16-30` == 1 & `diff_Depth_cm31-45_Depth_cm16-30` == 1, "#00CCFF", "black"),
                lfc5 = ifelse(`passed_ss_Depth_cm50-75_Depth_cm16-30` == 1 & `diff_Depth_cm50-75_Depth_cm16-30` == 1, "#00CCFF", "black"),
                lfc6 = ifelse(`passed_ss_Depth_cm50-75_Depth_cm31-45` == 1 & `diff_Depth_cm50-75_Depth_cm31-45` == 1, "#00CCFF", "black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc6, names_to = "group", values_to = "color") %>%
  dplyr::arrange(taxon)

# Combinar los dataframes df_fig_pair1 y df_fig_pair2
df_fig_pair <- df_fig_pair1 %>%
  dplyr::left_join(df_fig_pair2, by = c("taxon", "group"))

# Renombrar los grupos para hacerlos más legibles
df_fig_pair$group <- recode(df_fig_pair$group,
                            `lfc1` = "16-30 vs 0-15",
                            `lfc2` = "31-45 vs 0-15",
                            `lfc3` = "50-75 vs 0-15",
                            `lfc4` = "31-45 vs 16-30",
                            `lfc5` = "50-75 vs 16-30",
                            `lfc6` = "50-75 vs 31-45")

df_fig_pair$group <- factor(df_fig_pair$group,
                            levels = c("16-30 vs 0-15",
                                       "31-45 vs 0-15",
                                       "50-75 vs 0-15",
                                       "31-45 vs 16-30",
                                       "50-75 vs 16-30",
                                       "50-75 vs 31-45"))

# Filtrar las comparaciones que no tienen valores significativos
df_fig_pair <- df_fig_pair %>%
  group_by(group) %>%
  filter(any(!is.na(value))) %>%
  ungroup()

# Calcular los límites y el punto medio para la escala del heatmap
lo <- floor(min(df_fig_pair$value, na.rm = TRUE))
up <- ceiling(max(df_fig_pair$value, na.rm = TRUE))
mid <- (lo + up) / 2

# Crear el heatmap
fig_pair <- df_fig_pair %>%
  ggplot(aes(x = group, y = taxon, fill = value)) +
  geom_tile(color = "gray50") +
  scale_fill_gradient2(low = "#c05098", high = "#1BB6AFFF", mid = "gray90", #"#c04fa5"paletteer_d("MoMAColors::Doughton")
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = "LFC") +
  geom_text(aes(label = round(value, 2), color = color), size = 4, na.rm = TRUE) +
  scale_color_identity(guide = FALSE) +
  labs(x = "Depth Comparison", y = NULL) + #, title = "LogFold changes between depth levels") +
  theme_minimal() +
  theme(#plot.title = element_text(hjust = 0.5, size = 11),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 9))

# Mostrar el heatmap
print(fig_pair)








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
res_pair_location <- output_location$res_pair

# Visualizar los resultados
head(res_pair_location)

# Crear el dataframe df_fig_pair1 para los valores log-fold change (LFC)
df_fig_pair1_location <- res_pair_location %>%
  dplyr::filter(`diff_LocationEC` == 1 |
                  `diff_LocationLP` == 1 |
                  `diff_LocationLP_LocationEC` == 1) %>%
  dplyr::mutate(lfc1 = ifelse(`diff_LocationEC` == 1, round(`lfc_LocationEC`, 2), NA),
                lfc2 = ifelse(`diff_LocationLP` == 1, round(`lfc_LocationLP`, 2), NA),
                lfc3 = ifelse(`diff_LocationLP_LocationEC` == 1, round(`lfc_LocationLP_LocationEC`, 2), NA)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3, names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

# Crear el dataframe df_fig_pair2 para los valores significativos
df_fig_pair2_location <- res_pair_location %>%
  dplyr::filter(`diff_LocationEC` == 1 |
                  `diff_LocationLP` == 1 |
                  `diff_LocationLP_LocationEC` == 1) %>%
  dplyr::mutate(lfc1 = ifelse(`passed_ss_LocationEC` == 1 & `diff_LocationEC` == 1, "#00CCFF", "black"),
                lfc2 = ifelse(`passed_ss_LocationLP` == 1 & `diff_LocationLP` == 1, "#00CCFF", "black"),
                lfc3 = ifelse(`passed_ss_LocationLP_LocationEC` == 1 & `diff_LocationLP_LocationEC` == 1, "#00CCFF", "black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3, names_to = "group", values_to = "color") %>%
  dplyr::arrange(taxon)

# Combinar los dataframes df_fig_pair1 y df_fig_pair2
df_fig_pair_location <- df_fig_pair1_location %>%
  dplyr::left_join(df_fig_pair2_location, by = c("taxon", "group"))

# Renombrar los grupos para hacerlos más legibles
df_fig_pair_location$group <- recode(df_fig_pair_location$group,
                                     `lfc1` = "EC vs DM",
                                     `lfc2` = "LP vs DM",
                                     `lfc3` = "LP vs EC")

df_fig_pair_location$group <- factor(df_fig_pair_location$group,
                                     levels = c("EC vs DM",
                                                "LP vs DM",
                                                "LP vs EC"))

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
fig_pair_location <- df_fig_pair_location %>%
  ggplot(aes(x = group, y = taxon, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "purple3", high = "green4", mid = "gray90",
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = "LFC") +
  geom_text(aes(label = round(value, 2), color = color), size = 5, na.rm = TRUE) +
  scale_color_identity(guide = FALSE) +
  labs(x = "Location Comparison", y = NULL, title = "LogFold Change between locations") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9))

# Mostrar el heatmap
print(fig_pair_location)

# Gráfico de Barras Apiladas para Location
df_fig_pair_location %>%
  ggplot(aes(x = taxon, y = value, fill = group, label = round(value, 2))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(position = position_dodge(width = 0.9), vjust = -0.5) +
  scale_fill_manual(values = c("EC vs DM" = "#1f78b4", "LP vs DM" = "#33a02c", "LP vs EC" = "#e31a1c")) +
  labs(title = "Differential Abundance by Location",
       x = "Taxon",
       y = "Log-Fold Change",
       fill = "Comparison") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 16))

# Gráfico de Puntos para Location
df_fig_pair_location %>%
  ggplot(aes(x = group, y = taxon, size = abs(value), color = color)) +
  geom_point() +
  scale_color_identity() +
  labs(title = "Differential Abundance by Location",
       x = "Location Comparison",
       y = "Taxon",
       size = "Log-Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16))

# Como localidad no fue muy informativo hacer el grafico solo con depth

###############
#heatmap con buble plot
# Cargar las librerías necesarias
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)

loc_colors <- c("El Cacahuate" = "#2E8B57",
                "La Piedad" = "#74C476",
                "Dique Miguelito" =  "#628B2E")


# transformar a relativa
psrel <- transform_sample_counts(ps, function(x) (x / sum(x)) * 100)
# Filtrar los phyla diferencialmente abundantes
differential_phyla <- c("Zixibacteria", "WS2", "uncultured", "Thermoplasmatota", "TA06",
                        "Sva0485", "Spirochaetota", "RCP2-54", "Nitrospinota", "NB1-j",
                        "Myxococcota", "Modulibacteria", "Methylomirabilota", "MBNT15",
                        "LCP-89", "Latescibacterota", "Fibrobacterota", "Crenarchaeota",
                        "Chloroflexi", "Calditrichota", "Armatimonadota")

# Filtrar el objeto phyloseq para incluir solo los phyla diferencialmente abundantes
ps_filtered <- subset_taxa(psrel, Phylum %in% differential_phyla)

# Calcular la abundancia relativa
#ps_rel_abundance <- transform_sample_counts(ps_filtered, function(x) x / sum(x))

# Extraer datos para el bubble plot
# df_abundance <- psmelt(ps_rel_abundance)
df_abundance <- psmelt(ps_filtered)

# Resumir por profundidad y phylum
df_abundance_summary <- df_abundance %>%
  group_by(Phylum, Depth_cm) %>%
  summarise(Mean_Abundance = mean(Abundance * 100, na.rm = TRUE))  # Convertir a porcentaje

depth_colors <- c("0-15" = "#f5e5c4",
                  "16-30" = "#baa179",
                  "31-45" = "#a67451",
                  "50-75" = "#80673b")

bubble_plot <- ggplot(df_abundance_summary, aes(x = Depth_cm, y = Phylum, size = Mean_Abundance, fill = Depth_cm)) +
  geom_point(shape = 21, color = "black") +
  scale_size_continuous(range = c(1, 12), guide = "none") +  # Oculta la leyenda de tamaño
  scale_fill_manual(values = depth_colors) +  # Aplica la paleta de colores
  labs(#title = "Relative Abundance of Differentially Abundant Phyla by Depth",
       #x = "Depth (cm)",
       y = "Phylum",
       fill = "Depth (cm)") +  # Mantener solo la leyenda de fill
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    #legend.title = element_text(size = 10),  # Opcional: ajusta el tamaño del título de la leyenda
    legend.text = element_text(size = 9)    # Opcional: ajusta el tamaño del texto de la leyenda
  )


# Mostrar el bubble plot
print(bubble_plot)

# save
ggsave("results/final_plots/DA_phy_buble.pdf", bubble_plot)
ggsave("results/final_plots/DA_phy_ancomheatmap.pdf", fig_pair)


# Extra Crear un heatmap
heatmap_plot <- ggplot(df_abundance_summary, aes(x = Depth_cm, y = Phylum, fill = Mean_Abundance)) +
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Heatmap of Differentially Abundant Phyla by Depth",
       x = "Depth (cm)",
       y = "Phylum",
       fill = "Mean Abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Mostrar el heatmap
print(heatmap_plot)
df_abundance_summary
