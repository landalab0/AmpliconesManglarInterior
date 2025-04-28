# Alpha diversity

## Full

### Get Hill numbers

```{r}
library(hilldiv)
library(tidyverse)
library(kableExtra)

#
otu_table_raw <- otu_table(filtered_physeq2)
print(otu_table_omfg)

otu_table_df <- as.data.frame(otu_table_raw)

filtered_metadata <- as.data.frame(filtered_metadatatonga3)
#Get taxonomy table
tax_table_raw <- tax_table(filtered_physeq2)
print(tax_table_raw)

# Get hill numbers
q0 <- hill_div(otu_table_df, qvalue = 0)
q1 <- hill_div(otu_table_df, qvalue = 1)
q2 <- hill_div(otu_table_df, qvalue = 2)
#qf <- hill_div(otu_table_raw, qvalue = 0, tree = tree_umMPL)
# Merge metadata with Hill numbers
q012f_all <- cbind(q0, q1, q2) %>% as.data.frame() %>% rownames_to_column(var = "SampleID")
library(dplyr)

filtered_metadata2 <- filtered_metadata %>%
  rename(SampleID = BioSample)


metadata_with_hill <- q012f_all %>%
  inner_join(filtered_metadata, by = c("SampleID"="SampleID"))



#save table
write.table(metadata_with_hill, "Metadata_with_hill.tsv", quote = FALSE, sep = "\t", row.names = TRUE, col.names=TRUE)

# show
metadata_with_hill %>% head() %>% kable()
```

```{r}
metadata_with_hill %>% head() %>% kable()
```



### Explore Hill distribution

**Get means**

  ```{r}
# Reorder q values
meta_qs <- metadata_with_hill %>%
  pivot_longer(cols = q0:q2, names_to = "q", values_to = "value") %>%
  filter(q %in% c("q0", "q1", "q2")) %>%
  mutate(
    qs = case_when(
      q == "q0" ~ "q0=Observed",
      q == "q1" ~ "q1=Exp Shannon",
      q == "q2" ~ "q2=Inv Simpson",
    ))

#Get means of Hill numbers
means <- meta_qs %>% group_by(Study_zone) %>%
  summarise(means = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            .groups = 'drop')

#save table
write.table(means, "Hill_means_sd.tsv", quote = FALSE, sep = "\t", row.names = TRUE, col.names=TRUE)

print(means)
```

### Depth vs Hill Correlation

**Alpha diversity depth correlation to order q=0**

  ```{r}
#Get effective reads
Effective_reads <- colSums(otu_table_df[-1, ])
metadata_with_hill$Effective_reads <- Effective_reads
```


```{r}
#install.packages("ggpubr")
library(ggpubr)

q0_vs_depth <- ggscatter(metadata_with_hill,
                         x = "Effective_reads",
                         y = "q0",
                         xlab= "Sequencing depth (number of reads)",
                         ylab = "q=0",
                         #ylab="Alpha diversity q=0 (effective number of total ASVs)",
                         add = "reg.line",  # Add regression line
                         add.params = list(color = "#114e9d", fill = "lightgray"), # Customize reg. line
                         conf.int = TRUE, # Add confidence interval
                         cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                         cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)+
  theme(legend.title = element_blank(), legend.position = "none")

```


**Alpha diversity depth correlation to order q=1**

  ```{r}
q1_vs_depth <- ggscatter(metadata_with_hill,
                         x = "Effective_reads",
                         y = "q1",
                         xlab= "Sequencing depth (number of reads)",
                         ylab = "q=1",
                         #ylab="Alpha diversity q=1 (effective number of total ASVs)",
                         add = "reg.line",  # Add regression line
                         add.params = list(color = "#114e9d", fill = "lightgray"), # Customize reg. line
                         conf.int = TRUE, # Add confidence interval
                         cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                         cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)+
  theme(legend.title = element_blank(), legend.position = "none")

```


**Alpha diversity depth correlation to order q=2**

  ```{r}
q2_vs_depth <- ggscatter(metadata_with_hill,
                         x = "Effective_reads",
                         y = "q2",
                         xlab= "Sequencing depth (number of reads)",
                         ylab = "q=2",
                         #ylab="Alpha diversity q=2 (effective number of total ASVs)",
                         add = "reg.line",  # Add regression line
                         add.params = list(color = "#114e9d", fill = "lightgray"), # Customize reg. line
                         conf.int = TRUE, # Add confidence interval
                         cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                         cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)+
  theme(legend.title = element_blank(), legend.position = "none")

```

```{r}
library(ggplot2)
library(cowplot)

#save plots

#Combine plot
title_corr_plot <- ggdraw() + draw_label("Alpha diversity depth correlation with fltr samples")
correlation_plot_q012 <- plot_grid(title_corr_plot, q0_vs_depth, q1_vs_depth, q2_vs_depth, labels = c(" ","A", "B", "C"), ncol = 1, rel_heights = c(0.15, 1, 1, 1))
correlation_plot_q012
#save plot
#ggsave("../results/plots/alpha_depth_correlation_samples.png", correlation_plot_q012_, width = 8, height = 8)
```


### Check Hill numbers Normality

Shapiro test to Hill numbers


```{r}
library(ggplot2)
library(cowplot)

#Shapiro test
shapiro_q0 <- shapiro.test(metadata_with_hill$q0)
shapiro_q0_pvalue <- round(shapiro_q0$p.value, 5)
shapiro_q1 <- shapiro.test(metadata_with_hill$q1)
shapiro_q1_pvalue <- round(shapiro_q1$p.value, 5)
shapiro_q2 <- shapiro.test(metadata_with_hill$q2)
shapiro_q2_pvalue <- round(shapiro_q2$p.value, 5)

#Histograms
histplot_q0 <- ggplot(metadata_with_hill, aes(x = q0, xlab="q=0")) +
  geom_histogram(fill = "lightblue", bins = 15) +
  ggtitle(paste("Shapiro, p-value:", shapiro_q0_pvalue)) +
  theme_bw() + xlab("q=0") + ylab("Frequency")

histplot_q1 <- ggplot(metadata_with_hill, aes(x = q1)) +
  geom_histogram(fill = "lightblue", bins = 15) +
  ggtitle(paste("Shapiro, p-value:", shapiro_q1_pvalue)) +
  theme_bw() + xlab("q=1") + ylab("Frequency")

histplot_q2 <- ggplot(metadata_with_hill, aes(x = q2)) +
  geom_histogram(fill = "lightblue", bins = 15) +
  ggtitle(paste("Shapiro, p-value:", shapiro_q2_pvalue)) +
  theme_bw() + xlab("q=2") + ylab("Frequency")

#Combine plot
title_plot <- ggdraw() + draw_label("Histogram of Hill diversity") #, fontface = 'bold', x = 0.5, hjust = 0.5)
histplot_q012 <- plot_grid(title_plot, histplot_q0, histplot_q1, histplot_q2, labels = c(" ","A", "B", "C"), ncol = 1, rel_heights = c(0.15, 1, 1, 1))
histplot_q012
#save plot
#ggsave("../results/plots/hill_normality_test_histogram_samples.png", histplot_q012)
```


A pesar de que tienen una distribución normal, dado que los datos del proyecto tienen medidas repetidas y desbalanceo, lo adecuado es hacer comparaciones basadas en modelos. En este caso modelos lineales mixtos

### Plot Hill distribution


```{r}
library(ggpubr)
library(tidyverse)

# factor to reorder plot
metadata_with_hill$Study_zone <- factor(metadata_with_hill$Study_zone, levels = c("Rio San Pedro", "Celestún", "Laguna Cacahuate"))
metadata_with_hill$Depth  <- factor(metadata_with_hill$Depth_cm, levels = c("0-15", "16-30", "31-45", "50-75"))
#1f77b4

#plot
Studyzone_colors <- c("Rio San Pedro" = "#2ca02c", "Celestún" = "#1f77b4", "Laguna Cacahuate" = "#5cd194")

hill_barplot_explore <- metadata_with_hill %>%
  pivot_longer(cols = q0:q2, names_to = "q", values_to = "value") %>%
  filter(q %in% c("q0", "q1", "q2")) %>%
  mutate(
    qs = case_when(
      q == "q0" ~ "q0=Observed",
      q == "q1" ~ "q1=Exp Shannon",
      q == "q2" ~ "q2=Inv Simpson",
    )) %>%
  ggbarplot(
    x = "Depth",
    y = "value",
    add = "mean_se",
    facet.by = c("qs", "Study_zone"),
    fill = "Study_zone"
  ) +
  scale_fill_manual(values = Studyzone_colors) +
  geom_jitter(size = 0.8, position = position_jitter(width = 0.1)) +
  facet_grid(rows = vars(qs), cols = vars(Study_zone), scales = "free_y") +
  theme(
    strip.background =  element_blank(), #element_rect(fill = "")
    strip.text.x = element_text(face= "bold", size = 13, margin = margin(0.5, 0, 0.5, 0)),
    strip.text.y = element_text(face= "bold", size = 14, angle = 0, margin = margin(0, 0.5, 0, 0.5)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(colour = "black"),
    #panel.border = element_blank(),
    panel.spacing.x = unit(0.5, "lines"),
    panel.background = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(colour = "black", size = 11),
    legend.position = "bottom"
  ) +
  labs(
    fill = "Study_zone",
    y = "",
    x = "",
    title = "Biodiversity across Sites and Depths"
  ) +
  theme(legend.text = element_text(size = 12))

hill_barplot_explore

ggsave("Hill_diversity_across_Sites_and_Depths_fltr.pdf", hill_barplot_explore, width = 11, height = 11)

```

### Plot Alpha by depth with lmer

```{r}
library(lme4)
library(nlme)
library(lmerTest)
library(pgirmess)
install.packages("pgirmess")
library(emmeans)
library(easystats)
library(ggpubr)
library(multcomp)
```



```{r}
# Reorganizar los niveles de Location
metadata_with_hill$Study_zone <- relevel(metadata_with_hill$Study_zone, ref = "Rio San Pedro")

library(lme4)

# Asegúrate de que Study_zone es un factor
metadata_with_hill$Study_zone <- as.factor(metadata_with_hill$Study_zone)

# Ajustar el modelo lineal mixto
q0lmer <- lmer(q0 ~ Study_zone + (1 | ID), data = metadata_with_hill)

# Resumen del modelo
summary(q0lmer)

# Ajustar el modelo lineal mixto
q0lmer <- lmer(q0 ~ Study_zone + Depth + (1 | ), data = metadata_with_hill)
summary(q0lmer)

# Realizar la prueba de permutación con lmerTest
perm_test <- rand(q0lmer, nsim = 1000)

# Resumen de los resultados de la prueba de permutación
perm_test

#efectos
anova(q0lmer)
```

```{r}
library(easystats)
check_model(q0lmer)
```

```{r}
library(easystats)
report(q0lmer)
```

q0lmer_simple <- lmer(q0 ~ Study_zone + (1 | ID), data = metadata_with_hill)
q0_lmer_means_simple <- emmeans(q0lmer_simple, pairwise ~ Study_zone)
cld_results_simple <- cld(object = q0_lmer_means_simple$emmeans, Letters = letters)


```{r}
# Obtener las diferencias significativas usando emmeans
q0_lmer_means <- emmeans(q0lmer, pairwise ~ Study_zone)

# Obtener las letras de significancia
cld_results <- cld(object = q0_lmer_means$emmeans, Letters = letters)

# Convertir a data frame
emmeans_df <- as.data.frame(cld_results)
```

```{r}
# Crear el gráfico de barras con letras de significancia
q0plot <- ggplot(metadata_with_hill, aes(x = Study_zone, y = q0)) +
  geom_bar(data = emmeans_df, aes(y = emmean, fill = Study_zone), stat = "identity", position = position_dodge(), color = "black", width = 0.7) +
  geom_errorbar(data = emmeans_df, aes(y = emmean, ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.9)) +
  geom_jitter(aes(fill = Study_zone), width = 0.1, alpha = 0.4) +
  geom_text(data = emmeans_df, aes(y = emmean, label = .group), vjust = -9, color = "black", fontface = "bold", position = position_dodge(0.9)) +
  labs(title = "Richness", y = NULL, x = NULL ) + #,
  #caption = "letters by lmer emmeans") +
  theme_bw() + scale_fill_manual(values = Studyzone_colors) +
  theme(legend.position = "none")
q0plot
```

**q1**

  ```{r}

# Reorder Location
metadata_with_hill$Location <- relevel(metadata_with_hill$Location, ref = "Rio San Pedro")

# Fitted linear mixed model
q1lmer <- lmer(q1 ~ Study_zone + (1 | ID), data = metadata_with_hill)
summary(q1lmer)

# lmerTest
perm_testq1 <- rand(q1lmer, nsim = 1000)

# Summary
perm_testq1

# Efects
anova(q1lmer)

# Check model
library(easystats)

check_model(q1lmer)
report(q1lmer)

# Sig dif emmeans
q1_lmer_means <- emmeans(q1lmer, pairwise ~ Study_zone)

# letters
cld_resultsq1 <- cld(object = q1_lmer_means$emmeans, Letters = letters)

# Convert df
emmeans_dfq1 <- as.data.frame(cld_resultsq1)

# plot
q1plot <- ggplot(metadata_with_hill, aes(x = Study_zone, y = q1)) +
  geom_bar(data = emmeans_dfq1, aes(y = emmean, fill = Study_zone), stat = "identity", position = position_dodge(), color = "black", width = 0.7) +
  geom_errorbar(data = emmeans_dfq1, aes(y = emmean, ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.9)) +
  geom_jitter(aes(fill = Study_zone), width = 0.1, alpha = 0.4) +
  geom_text(data = emmeans_dfq1, aes(y = emmean, label = .group), vjust = -9, color = "black", fontface = "bold",position = position_dodge(0.9)) +
  labs(title = "Exp Shannon", y = NULL, x = NULL ) + #,
  #caption = "letters by lmer emmeans") +
  theme_bw() + scale_fill_manual(values = Studyzone_colors) +
  theme(legend.position = "none")
q1plot
```


**q2**

  ```{r}

# Reorder Location
metadata_with_hill$Study_zone <- relevel(metadata_with_hill$Study_zone, ref = "Laguna Cacahuate")

# Fitted linear mixed model
q2lmer <- lmer(q2 ~ Study_zone + (1 | ID), data = metadata_with_hill)
summary(q2lmer)

# lmerTest
perm_testq2 <- rand(q2lmer, nsim = 1000)

# Summary
perm_testq2

# Efects
anova(q2lmer)

# Check model
library(easystats)

check_model(q2lmer)
report(q2lmer)

# Sig dif emmeans
q2_lmer_means <- emmeans(q2lmer, pairwise ~ Study_zone)

# letters
cld_resultsq2 <- cld(object = q2_lmer_means$emmeans, Letters = letters)

# Convert df
emmeans_dfq2 <- as.data.frame(cld_resultsq2)

# plot
q2plot <- ggplot(metadata_with_hill, aes(x = Study_zone, y = q2)) +
  geom_bar(data = emmeans_dfq2, aes(y = emmean, fill = Study_zone), stat = "identity", position = position_dodge(), color = "black", width = 0.7) +
  geom_errorbar(data = emmeans_dfq2, aes(y = emmean, ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.9)) +
  geom_jitter(aes(fill = Study_zone), width = 0.1, alpha = 0.4) +
  geom_text(data = emmeans_dfq2, aes(y = emmean, label = .group), vjust = -9, color = "black", fontface = "bold",position = position_dodge(0.9)) +
  labs(title = "Inv Simpson", y = NULL, x = NULL ) + #,
  #caption = "letters by lmer emmeans") +
  theme_bw() + scale_fill_manual(values = Studyzone_colors) +
  theme(legend.position = "none")
q2plot
```
library (cowplot)
plot_grid(q0plot, q1plot, q2plot, labels = c("A", "B", "C"), ncol = 2)

**qf**

  ```{r}

# Reorder Location
metadata_with_hill$Location <- relevel(metadata_with_hill$Location, ref = "El Cacahuate")

# Fitted linear mixed model
qflmer <- lmer(qf ~ Depth_cm + Location + Site + (1 | ID), data = metadata_with_hill)
summary(qflmer)

# lmerTest
perm_testqf <- rand(qflmer, nsim = 1000)

# Summary
perm_testqf

# Efects
anova(qflmer)

# Check model
library(easystats)

check_model(qflmer)
report(qflmer)

# Sig dif emmeans
qf_lmer_means <- emmeans(qflmer, pairwise ~ Depth_cm)

# letters
cld_resultsqf <- cld(object = qf_lmer_means$emmeans, Letters = letters)

# Convert df
emmeans_dfqf <- as.data.frame(cld_resultsqf)

# plot
qfplot <- ggplot(metadata_with_hill, aes(x = Depth_cm, y = qf)) +
  geom_bar(data = emmeans_dfqf, aes(y = emmean, fill = Depth_cm), stat = "identity", position = position_dodge(), color = "black", width = 0.7) +
  geom_errorbar(data = emmeans_dfqf, aes(y = emmean, ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.9)) +
  geom_jitter(aes(fill = Depth_cm), width = 0.1, alpha = 0.4) +
  geom_text(data = emmeans_dfqf, aes(y = emmean, label = .group), vjust = -8, color = "black", fontface = "bold", position = position_dodge(0.9)) +
  labs(title = "PD Faith", y = NULL, x = NULL ) + #,
  #caption = "letters by lmer emmeans") +
  theme_bw() + scale_fill_manual(values = depth_colors) +
  theme(legend.position = "none")
qfplot
```

```{r}
library(ggplot2)
library(cowplot)

#Combine plot
#title_alpha_plot <- ggdraw() + draw_label("Alpha Diversity by Depth" , fontface = "bold")
ytitle <- ggdraw() + draw_label("Effective number of ASVs", angle = 90, size = 14 )
q01plot <- plot_grid(q0plot, q1plot, labels = c("A", "B"), ncol = 2)
q2fplot <- plot_grid(q2plot, qfplot, labels = c("C", "D"), ncol = 2)
alphaplots <- plot_grid(q01plot, q2fplot, ncol = 1, rel_heights = c(1, 1))
alphaplots_y <- plot_grid(ytitle, alphaplots, ncol = 2, rel_widths = c(0.05, 1))
alphaplots_y
#save plot
#ggsave("../results/plots/alpha_depth_correlation_samples.png", correlation_plot_q012_, width = 8, height = 8)
```


### Alpha by Location with lmer

**q0**


  ```{r}
# Obtener las diferencias significativas usando emmeans
q0_lmer_means_loc <- emmeans(q0lmer, pairwise ~ Location)

# Obtener las letras de significancia
q0cld_results_loc <- cld(object = q0_lmer_means_loc$emmeans, Letters = letters)

# Convertir a data frame
q0emmeans_df_loc <- as.data.frame(q0cld_results_loc)
q0emmeans_df_loc
```


```{r}
# Obtener las diferencias significativas usando emmeans
q1lmer_means_loc <- emmeans(q1lmer, pairwise ~ Location)

# Obtener las letras de significancia
q1cld_results_loc <- cld(object = q1lmer_means_loc$emmeans, Letters = letters)

# Convertir a data frame
q1emmeans_df_loc <- as.data.frame(q1cld_results_loc)
q1emmeans_df_loc
```

```{r}
# Obtener las diferencias significativas usando emmeans
q2lmer_means_loc <- emmeans(q2lmer, pairwise ~ Location)

# Obtener las letras de significancia
q2cld_results_loc <- cld(object = q2lmer_means_loc$emmeans, Letters = letters)

# Convertir a data frame
q2emmeans_df_loc <- as.data.frame(q2cld_results_loc)
q2emmeans_df_loc
```
```{r}
# Obtener las diferencias significativas usando emmeans
qflmer_means_loc <- emmeans(qflmer, pairwise ~ Location)

# Obtener las letras de significancia
qfcld_results_loc <- cld(object = qflmer_means_loc$emmeans, Letters = letters)

# Convertir a data frame
qfemmeans_df_loc <- as.data.frame(qfcld_results_loc)
qfemmeans_df_loc
```
