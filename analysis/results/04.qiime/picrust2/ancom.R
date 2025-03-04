# Load package
library(qiime2R)
library(ANCOMBC)
library(phyloseq)
library(dplyr)
library(tidyverse)
library(DT)

# phyloseq object
ps <- qza_to_phyloseq(
  features = "ko_metagenome.qza",
  metadata = "../../../data/metadata.tsv")

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
output_KO <- ancombc2(
  data = tse,
  assay_name = "counts",
  tax_lev = NULL,
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


saveRDS(output_KO, file = "output_KO.rds")
