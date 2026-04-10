# =============================================================================
# 04_luminal_signature_scoring.R
# Compute LumA / LumP activity scores and integrate ML classifier predictions
# for all cohorts: Wouter (Science 2020), ML001, ML002, MW3, MW4, MW6, MW7, MW8
#
# Scoring method: expression-based z-score with bin-matched null distribution
# (adapted from Karthaus et al. eLife 2020)
#
# Input:  Luminal Seurat objects + classifier probability CSVs
# Output: Per-cohort signature data frames saved as .RDS / .CSV
# =============================================================================

source("R/00_setup.R")
library(ggnewscale)

# =============================================================================
# Core scoring function
# =============================================================================

#' Calculate LumA and LumP activity scores using a bin-matched null distribution
#'
#' Genes are split into 20 expression bins; for each simulation, a matched
#' null gene is randomly drawn from the same bin. The final score is the
#' z-score of the raw gene-set score relative to the null distribution.
#'
#' @param expression_data Matrix (genes × cells)
#' @param gene_sets        Character vector of marker genes
#' @param num_bins         Number of expression bins (default 20)
#' @param num_simulations  Number of null simulations (default 100; use 2000 for publication)
#' @return Numeric vector of activity scores (one per cell)
calculate_gene_scores <- function(expression_data,
                                  gene_sets,
                                  num_bins        = 20,
                                  num_simulations = 100) {

  # Remove silent genes
  expressed_genes  <- rowSums(expression_data != 0) > 0
  expression_data  <- expression_data[expressed_genes, ]

  # Normalize to [0, 1] via logistic transform on z-scores
  normalized_data  <- expression_data %>%
    scale(center = TRUE, scale = TRUE) %>%
    plogis()

  # Bin genes by mean expression
  gene_means <- rowMeans(expression_data) %>%
    tibble::enframe(name = "gene", value = "mean_expression") %>%
    mutate(bin = ntile(mean_expression, num_bins))
  bins              <- gene_means$bin
  names(bins)       <- gene_means$gene

  # Inner function: compute raw score + null z-score for a gene set
  calculate_z_score <- function(gene_set) {
    gene_set          <- intersect(gene_set, rownames(expression_data))
    bins_of_gene_set  <- bins[gene_set]
    raw_scores        <- colMeans(normalized_data[gene_set, , drop = FALSE])

    null_scores_matrix <- matrix(nrow = ncol(normalized_data), ncol = num_simulations)
    for (i in seq_len(num_simulations)) {
      sampled_genes <- sapply(bins_of_gene_set, function(bin_id) {
        genes_in_bin <- gene_means %>% filter(bin == bin_id)
        if (nrow(genes_in_bin) > 0) sample(genes_in_bin$gene, 1) else NA
      })
      if (any(is.na(sampled_genes))) {
        null_scores_matrix[, i] <- NA
      } else {
        null_scores_matrix[, i] <- colMeans(normalized_data[sampled_genes, , drop = FALSE])
      }
    }

    null_mean <- rowMeans(null_scores_matrix, na.rm = TRUE)
    null_sd   <- apply(null_scores_matrix, 1, sd, na.rm = TRUE)
    null_sd[null_sd == 0] <- 1e-10   # prevent division by zero
    (raw_scores - null_mean) / null_sd
  }

  calculate_z_score(gene_sets)
}


# =============================================================================
# Helper: add classifier probabilities and compute PredictionProb
# =============================================================================

#' Merge classifier output into a signature data frame
#'
#' @param sig_df      Data frame with LumType, LumA/LumP scores
#' @param pred_csv    Path to Python classifier output CSV
#' @param lum_type_col Column name for cell type label (default "LumType")
merge_classifier <- function(sig_df, pred_csv, lum_type_col = "LumType") {
  pred <- read.csv(pred_csv)
  sig_df <- cbind(sig_df, pred)
  sig_df$LumA_pro <- sig_df$Epi_Luminal_1
  sig_df$LumP_pro <- sig_df$Epi_Luminal_2Psca
  sig_df <- sig_df %>%
    mutate(PredictionProb = case_when(
      .data[[lum_type_col]] == "LumA" ~ LumA_pro,
      TRUE                             ~ LumP_pro
    ))
  return(sig_df)
}


# =============================================================================
# Wouter Science cohort (mouse)
# =============================================================================
Wouter_science_luminal <- readRDS("Wouter_science_luminal.filter.RDS")

Luminal_A.scores <- calculate_gene_scores(Wouter_science_luminal@assays$RNA@data,
                                           gene_sets = LumA_marker_lab, num_simulations = 2000)
Luminal_P.scores <- calculate_gene_scores(Wouter_science_luminal@assays$RNA@data,
                                           gene_sets = LumP_marker_lab, num_simulations = 2000)

Wouter_science_luminal_meta           <- Wouter_science_luminal@meta.data
Wouter_science_luminal_meta$Luminal_A <- Luminal_A.scores
Wouter_science_luminal_meta$Luminal_P <- Luminal_P.scores

# Merge ML classifier output
Wouter_science_luminal.intact <- subset(Wouter_science_luminal, timepoint == "T00_Intact")
Wouter_science_luminal.other  <- subset(Wouter_science_luminal, timepoint != "T00_Intact")

Wouter_Intact_pred  <- read.csv("data/Wouter_science_luminal_Intact_prediction.csv")
Wouter_OtherTP_pred <- read.csv("data/Wouter_science_luminal_notIntact_prediction_result.csv")
Wouter_Intact_pred$NAME  <- colnames(Wouter_science_luminal.intact)
Wouter_OtherTP_pred$NAME <- colnames(Wouter_science_luminal.other)

Wouter_luminal_pred <- rbind(Wouter_Intact_pred, Wouter_OtherTP_pred[, -1])
Wouter_science_luminal_meta <- left_join(Wouter_science_luminal_meta, Wouter_luminal_pred)
Wouter_science_luminal_meta$LumA_pro <- Wouter_science_luminal_meta$Epi_Luminal_1
Wouter_science_luminal_meta$LumP_pro <- Wouter_science_luminal_meta$Epi_Luminal_2Psca

Wouter_science_luminal_meta <- Wouter_science_luminal_meta %>%
  mutate(
    LumType = case_when(FullType == "Epi_Luminal_1" ~ "LumA",
                        FullType == "Epi_Luminal_2Psca" ~ "LumP"),
    PredictionProb = case_when(LumType == "LumA" ~ LumA_pro, TRUE ~ LumP_pro)
  )
saveRDS(Wouter_science_luminal_meta, "Wouter_science_luminal_signature.RDS")


# =============================================================================
# ML001
# =============================================================================
ML001_Luminal_clean <- readRDS("ML001_Luminal_clean.RDS")

ML001_meta <- ML001_Luminal_clean@meta.data %>%
  mutate(
    LumType = case_when(Cluster == "LumP" ~ "LumP", TRUE ~ "LumA"),
    timepoint = case_when(HTO_maxID %in% c("B301", "B303") ~ "1w", TRUE ~ "4w"),
    CellType  = case_when(HTO_maxID %in% c("B301", "B302") ~ "Luminal AR KO",
                          TRUE ~ "Stromal AR KO")
  )
ML001_meta$CellName <- rownames(ML001_meta)
ML001_LumA_scores   <- calculate_gene_scores(ML001_Luminal_clean@assays$RNA@data,
                                              LumA_marker_lab, num_simulations = 2000)
ML001_LumP_scores   <- calculate_gene_scores(ML001_Luminal_clean@assays$RNA@data,
                                              LumP_marker_lab, num_simulations = 2000)
ML001_meta$Luminal_A <- ML001_LumA_scores
ML001_meta$Luminal_P <- ML001_LumP_scores

ML001_meta <- merge_classifier(ML001_meta, "data/ML001_luminal_prediction.csv")
saveRDS(ML001_meta, "ML001_Luminal_signature_df_lab.RDS")


# =============================================================================
# ML002
# =============================================================================
ML002_Luminal_clean <- readRDS("ML002_Luminal_clean.RDS")

ML002_meta <- ML002_Luminal_clean@meta.data %>%
  mutate(
    LumType   = case_when(Cluster == "LumP" ~ "LumP", TRUE ~ "LumA"),
    timepoint = case_when(HTO_maxID %in% c("B301", "B303") ~ "1w", TRUE ~ "4w"),
    CellType  = case_when(HTO_maxID %in% c("B301", "B302") ~ "Luminal Acsl4 KO",
                          TRUE ~ "Luminal Atg7 KO")
  )
ML002_meta$CellName <- rownames(ML002_meta)
ML002_LumA_scores   <- calculate_gene_scores(ML002_Luminal_clean@assays$RNA@data,
                                              LumA_marker_lab, num_simulations = 2000)
ML002_LumP_scores   <- calculate_gene_scores(ML002_Luminal_clean@assays$RNA@data,
                                              LumP_marker_lab, num_simulations = 2000)
ML002_meta$Luminal_A <- ML002_LumA_scores
ML002_meta$Luminal_P <- ML002_LumP_scores
ML002_meta <- merge_classifier(ML002_meta, "data/ML002_luminal_prediction.csv")
saveRDS(ML002_meta, "ML002_Luminal_signature_df_lab.RDS")


# =============================================================================
# MW3
# =============================================================================
MW3_luminal <- readRDS("MW3_seurat_object.clean.RDS")

MW3_meta <- MW3_luminal@meta.data %>%
  mutate(
    LumType  = case_when(cluster == "LumP" ~ "LumP", TRUE ~ "LumA"),
    CellType = case_when(
      HTO_maxID == "B301" ~ "AR-intact",
      HTO_maxID == "B302" ~ "AR Castration 4w",
      HTO_maxID == "B303" ~ "Stromal and Luminal AR KO",
      HTO_maxID == "B304" ~ "FVB WT cas 3d"
    )
  )
MW3_LumA_scores <- calculate_gene_scores(GetAssayData(MW3_luminal), LumA_marker_lab,
                                          num_simulations = 2000)
MW3_LumP_scores <- calculate_gene_scores(GetAssayData(MW3_luminal), LumP_marker_lab,
                                          num_simulations = 2000)
MW3_meta$Luminal_A <- MW3_LumA_scores
MW3_meta$Luminal_P <- MW3_LumP_scores
MW3_meta <- merge_classifier(MW3_meta, "data/MW3_luminal_prediction.csv")
saveRDS(MW3_meta, "MW3_Luminal_signature_df.RDS")


# =============================================================================
# MW4
# =============================================================================
MW4_luminal <- readRDS("MW4_luminal.RDS")
MW4_luminal@meta.data <- MW4_luminal@meta.data %>%
  mutate(
    LumType  = case_when(cluster == "LumP" ~ "LumP", TRUE ~ "LumA"),
    CellType = case_when(
      HTO_maxID == "B301" ~ "Atg7-WT intact",
      HTO_maxID == "B302" ~ "Atg7-WT intact Castration 4w",
      HTO_maxID == "B303" ~ "Atg7-Ko intact",
      HTO_maxID == "B304" ~ "NKX3.1+/+ 4w",
      HTO_maxID == "B305" ~ "NKX3.1-/- 4w"
    )
  )
MW4_meta        <- MW4_luminal@meta.data
MW4_LumA_scores <- calculate_gene_scores(GetAssayData(MW4_luminal), LumA_marker_lab,
                                          num_simulations = 2000)
MW4_LumP_scores <- calculate_gene_scores(GetAssayData(MW4_luminal), LumP_marker_lab,
                                          num_simulations = 2000)
MW4_meta$Luminal_A <- MW4_LumA_scores
MW4_meta$Luminal_P <- MW4_LumP_scores
MW4_meta <- merge_classifier(MW4_meta, "data/MW4_luminal_prediction.csv")
saveRDS(MW4_meta, "MW4_Luminal_signature_df.RDS")


# =============================================================================
# MW8
# =============================================================================
MW8_seurat_object <- readRDS("MW8_seurat_object.RDS")
MW8_luminal       <- subset(MW8_seurat_object, ClusterUpdate %in% c("LumA", "LumP"))

MW8_meta        <- MW8_luminal@meta.data
MW8_LumA_scores <- calculate_gene_scores(GetAssayData(MW8_luminal), LumA_marker_lab,
                                          num_simulations = 2000)
MW8_LumP_scores <- calculate_gene_scores(GetAssayData(MW8_luminal), LumP_marker_lab,
                                          num_simulations = 2000)
MW8_meta$Luminal_A <- MW8_LumA_scores
MW8_meta$Luminal_P <- MW8_LumP_scores

MW8_luminal_pred <- read.csv("data/MW8_luminal_prediction.csv")
MW8_meta <- cbind(MW8_meta, MW8_luminal_pred)
MW8_meta$LumA_pro <- MW8_meta$Epi_Luminal_1
MW8_meta$LumP_pro <- MW8_meta$Epi_Luminal_2Psca
MW8_meta$Batch    <- factor(MW8_meta$Batch,
                             levels = c("PTN_KO_4w", "PTN_NKX_DKO_4w",
                                        "RAPTOR_KO_4w", "RAPTOR_AR_DKO_4w"))
MW8_meta <- MW8_meta %>%
  mutate(PredictionProb = case_when(ClusterUpdate == "LumA" ~ LumA_pro, TRUE ~ LumP_pro))
saveRDS(MW8_meta, "MW8_Luminal_signature_df.RDS")


# =============================================================================
# Classifier export helper (for sending to Python ML pipeline)
# =============================================================================

#' Scale and export luminal expression data for ML classification
#'
#' @param seurat_obj  Seurat object
#' @param common_genes Common gene list used across all cohorts
#' @return Data frame with scaled expression + CellType column
export_scaled_seurat_data <- function(seurat_obj, common_genes) {
  data_matrix    <- GetAssayData(
    subset(seurat_obj, features = common_genes), layer = "counts"
  )
  expressed_genes <- rowSums(data_matrix != 0) > 0
  data_matrix     <- data_matrix[expressed_genes, ]
  data_matrix     <- t(t(data_matrix) / (colSums(data_matrix) / 1e6))  # CPM
  scaled_data     <- scale(data_matrix, center = TRUE, scale = TRUE)

  df          <- as.data.frame(t(scaled_data))
  df$CellType <- seurat_obj@meta.data$FullType[
    match(rownames(df), rownames(seurat_obj@meta.data))
  ]
  return(df)
}

# Example usage:
# common_genes <- readRDS("data/ML_MW_Wouter_common_gene_list.RDS")
# Wouter_intact_scaled <- export_scaled_seurat_data(Wouter_science_luminal.intact, common_genes)
# write.csv(Wouter_intact_scaled, "data/Wouter_science_luminal_Intact_scaled.csv", row.names = FALSE)
