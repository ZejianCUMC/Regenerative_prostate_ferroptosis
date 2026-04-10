# =============================================================================
# 00_setup.R
# Global setup: library loading, color palettes, and parallel settings
# =============================================================================

# --- Core libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(tidyverse)
  library(stringr)
  library(data.table)
  library(ggpubr)
  library(MAST)
  library(rstatix)
  library(harmony)
  library(DropletUtils)
  library(ggnewscale)
  library(tibble)
  library(ggplot2)
})

# --- Parallel processing ------------------------------------------------------
library(future)
plan(multicore, workers = 4)
options(future.globals.maxSize = 199 * 1024^3)
options(future.seed = TRUE)

# --- Global color palette -----------------------------------------------------
# 43-color palette used for cluster visualization throughout the project
col.cluster2.1 <- c(
  '#1C79B7','#F38329','#2DA248','#DC403E','#976BA6',
  '#8D574C','#D07DB0','#BDBF3B','#27BDD0','#B0C7E5',
  '#F9BA79','#A0CC8A','#F09EA1','#C7B0D1','#D4BCB8',
  '#F7CCDD','#DBDD8D','#A3D7E5','#B04E4F','#A38A58',
  '#ED5351','#0C8945','#27251f','#A82764','#F8DAE4',
  '#7A4B1F','#5E6BAE','#8AC997','#DAC9AC','#0F5148',
  '#A0B5DC','#9F858E','#5C181D','#7B8380','#E8DDDA',
  '#264220','#5AB747','#5169AE','#4B3D58','#CD428A',
  '#62615A','#B82129','#66762E'
)

# --- Luminal marker gene sets (Karthaus et al. eLife markers) ----------------
# Used throughout luminal signature scoring
LumA_marker_lab <- c(
  "Gsdma", "Tgm4", "Krt8", "Cd24a", "Pbsn",
  "Hoxb13", "Ceacam1", "Prom1", "Nkx3-1"
)

LumP_marker_lab <- c(
  "Cldn10", "Lrrc26", "Ppp1r1b", "Krt4", "Wfdc2",
  "Krt7", "Tacstd2", "Clu", "Ly6a", "Krt8", "Cd24a", "Ceacam1"
)

# --- MUFA-related ferroptosis genes (human) -----------------------------------
MUFA_gene_human <- c(
  'ACSL3', 'SCD', 'MBOAT2', 'GPX4', 'FASN',
  'ACACA', 'ACSL4', 'ELOVL2', 'ELOVL6', 'FADS1', 'FADS2'
)

# --- PTN receptor gene set ----------------------------------------------------
PTN_receptor <- c('PTPRZ1', 'ALK', 'CSPG5', 'LRP1', 'ITGA3', 'ITGB3', 'SDC3')
