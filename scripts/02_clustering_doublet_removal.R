# =============================================================================
# 02_clustering_doublet_removal.R
# Per-sample clustering, DoubletFinder-based doublet removal,
# manual cell-type annotation, and luminal cell extraction.
#
# Samples: ML001, ML002, MW3, MW4, MW8
#
# Input:  Per-sample Seurat objects from 01_data_loading_hto_demux.R
# Output: Cleaned per-sample Seurat objects + luminal subsets
# =============================================================================

source("R/00_setup.R")
suppressPackageStartupMessages(library(DoubletFinder))


# =============================================================================
# Generic helper functions
# =============================================================================

#' Standard Seurat pre-processing: PCA + UMAP + tSNE + neighbor graph
#'
#' @param obj         Seurat object (must contain variable features)
#' @param dims        Integer vector of PCA dims to use (e.g. 1:15)
#' @param seed        Random seed for reproducibility
#' @param assay       Assay to use ("RNA" or "SCT")
run_dim_reduction <- function(obj, dims = 1:15, seed = 1234, assay = "RNA") {
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "pca", dims = dims,
                 verbose = FALSE, seed.use = seed)
  obj <- RunTSNE(obj, reduction = "pca", dims = dims,
                 verbose = FALSE, seed.use = seed)
  obj <- FindNeighbors(obj, reduction = "pca", dims = dims, verbose = FALSE)
  return(obj)
}

#' DoubletFinder wrapper: estimate doublet rate, sweep pK, classify doublets
#'
#' @param obj         Seurat object (pre-processed)
#' @param annotations Character vector of cell-type annotations per cell
#' @param pcs         PCs to use (default 1:40)
#' @param use_sct     Logical: was SCTransform used?
run_doublet_finder <- function(obj, annotations, pcs = 1:40, use_sct = FALSE) {
  n_cells      <- ncol(obj)
  doublet_rate <- n_cells * 8e-6           # 10x Genomics per-cell doublet rate

  homotypic_prop <- modelHomotypic(annotations)
  nExp           <- round(doublet_rate * n_cells)
  nExp_adj       <- round(nExp * (1 - homotypic_prop))

  # Parameter sweep for optimal pK
  sweep_res   <- paramSweep(obj, PCs = pcs, sct = use_sct)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn       <- find.pK(sweep_stats)
  opt_pK      <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  # Run DoubletFinder twice (standard + homotypic-adjusted)
  obj <- doubletFinder(obj, PCs = pcs, pN = 0.25, pK = opt_pK,
                       nExp = nExp,     reuse.pANN = FALSE, sct = use_sct)
  obj <- doubletFinder(obj, PCs = pcs, pN = 0.25, pK = opt_pK,
                       nExp = nExp_adj, reuse.pANN = FALSE, sct = use_sct)

  # Rename auto-generated columns for clarity
  n_col <- ncol(obj@meta.data)
  colnames(obj@meta.data)[n_col]     <- 'DF.classifications.adj'
  colnames(obj@meta.data)[n_col - 2] <- 'DF.classifications'

  # Confidence-based final call
  obj@meta.data <- obj@meta.data %>%
    mutate(
      Ann_final = case_when(
        DF.classifications.adj == "Singlet" & DF.classifications == "Singlet" ~ "Singlet",
        DF.classifications.adj == "Doublet" & DF.classifications == "Doublet" ~ "Doublet High Conf.",
        TRUE ~ "Doublet-Low Conf."
      )
    )

  message("DoubletFinder results: ", table(obj$Ann_final))
  return(obj)
}


# =============================================================================
# ML001 – RNA normalization pipeline (no SCT)
# =============================================================================
ML001_seurat_obj <- readRDS("ML001_seurat_obj_HTO.RDS")
ML001_seurat_obj.subset <- subset(ML001_seurat_obj,
                                   HTO_classification.global == "Singlet")
ML001_seurat_obj.subset <- ML001_seurat_obj.subset %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 4000, verbose = FALSE) %>%
  ScaleData(features = rownames(.))

ML001_seurat_obj.subset <- run_dim_reduction(ML001_seurat_obj.subset, dims = 1:14)
ML001_seurat_obj.subset <- FindClusters(ML001_seurat_obj.subset,
                                         resolution = seq(0.2, 1.0, by = 0.1),
                                         verbose = FALSE)

# Cluster annotation at res 0.3
cluster_mapping_ML001 <- c(
  "0" = "LumA", "1" = "Lum_not_sure", "2" = "LumA", "3" = "Basel",
  "4" = "LumA", "5" = "DeadCell",    "6" = "Lum1_pbsnhigh",
  "7" = "LumP", "8" = "Mesenchymal", "9" = "Macrophages", "10" = "Tcell"
)
ML001_seurat_obj.subset@meta.data$SeuratCluster <-
  as.character(ML001_seurat_obj.subset@meta.data$RNA_snn_res.0.3)
ML001_seurat_obj.subset@meta.data$CellType <-
  cluster_mapping_ML001[ML001_seurat_obj.subset@meta.data$SeuratCluster]

# Doublet removal
ML001_seurat_obj.subset <- run_doublet_finder(
  ML001_seurat_obj.subset,
  annotations = ML001_seurat_obj.subset@meta.data$CellType,
  use_sct = FALSE
)
ML001_seurat_obj.subset <- subset(ML001_seurat_obj.subset, Ann_final == "Singlet")
saveRDS(ML001_seurat_obj.subset, "ML001_clean.RDS")


# =============================================================================
# ML002 – RNA normalization pipeline (no SCT)
# =============================================================================
ML002_seurat_obj <- readRDS("ML002_seurat_obj_HTO.RDS")
ML002_seurat_obj.subset <- subset(ML002_seurat_obj,
                                   HTO_classification.global == "Singlet") %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2500, verbose = FALSE) %>%
  ScaleData(features = rownames(.))

ML002_seurat_obj.subset <- run_dim_reduction(ML002_seurat_obj.subset, dims = 1:15)
ML002_seurat_obj.subset <- FindClusters(ML002_seurat_obj.subset,
                                         resolution = seq(0.1, 1.0, by = 0.1),
                                         verbose = FALSE)

cluster_mapping_ML002 <- c(
  "0" = "LumA",          "1" = "Lum_not_sure_1", "2" = "Lum_not_sure_1",
  "3" = "Basel",         "4" = "LumP",            "5" = "Lum_not_sure_2",
  "6" = "Lum_not_sure_2","7" = "Macrophages",     "8" = "Mesenchymal",
  "9" = "Tcell"
)
ML002_seurat_obj.subset@meta.data$CellType <-
  cluster_mapping_ML002[as.character(ML002_seurat_obj.subset@meta.data$RNA_snn_res.0.3)]

# Find markers and export
ML002.markers <- FindAllMarkers(ML002_seurat_obj.subset, only.pos = TRUE,
                                 min.pct = 0.25, logfc.threshold = 0.25,
                                 group.by = "RNA_snn_res.0.3")
write.csv(ML002.markers, "output/ML002_cluster_markers.csv", row.names = FALSE)

ML002_seurat_obj.subset <- run_doublet_finder(
  ML002_seurat_obj.subset,
  annotations = ML002_seurat_obj.subset@meta.data$CellType,
  use_sct = FALSE
)
ML002_seurat_obj.subset <- subset(ML002_seurat_obj.subset, Ann_final == "Singlet")
saveRDS(ML002_seurat_obj.subset, "ML002_clean.RDS")


# =============================================================================
# Integrate ML001 + ML002 → extract luminal cells
# =============================================================================
ML001_seurat_obj.subset@meta.data$orig.ident   <- "ML001"
ML002_seurat_obj.subset@meta.data$orig.ident   <- "ML002"
ML001_seurat_obj.subset@meta.data$SampleCode   <-
  paste0("ML001_", ML001_seurat_obj.subset@meta.data$HTO_maxID)
ML002_seurat_obj.subset@meta.data$SampleCode   <-
  paste0("ML002_", ML002_seurat_obj.subset@meta.data$HTO_maxID)

merged_ML <- merge(ML001_seurat_obj.subset, y = ML002_seurat_obj.subset,
                   add.cell.ids = c("001", "002"), project = "ShenLab")

# SCTransform + Harmony integration
merged_ML <- merged_ML %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  SCTransform(vars.to.regress = "percent.mt")

merged_ML <- RunPCA(merged_ML, assay = "SCT", npcs = 50)
ML001_and_002_harmonized <- RunHarmony(merged_ML,
                                        group.by.vars = "orig.ident",
                                        reduction    = "pca",
                                        assay.use    = "SCT",
                                        reduction.save = "harmony")
ML001_and_002_harmonized <- RunUMAP(ML001_and_002_harmonized,
                                     reduction = "harmony", assay = "SCT", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = seq(0.1, 1.2, by = 0.1))

# Luminal cell extraction (after manual annotation of integrated object)
# Note: annotation requires inspection of Integrated.markers output
Integrated.markers <- FindAllMarkers(ML001_and_002_harmonized, assay = "SCT",
                                      only.pos = TRUE, min.pct = 0.25,
                                      logfc.threshold = 0.25, group.by = "SCT_snn_res.0.3")
write.csv(Integrated.markers, "output/ML001_ML002_integrated_markers.csv", row.names = FALSE)
# → Manually annotate clusters, add PutativCellType to xlsx, then reload:
# ML001_and_002_harmonized_celltype <- readxl::read_xlsx("Integrated.markers.xlsx", sheet = 2)
# ML001_and_002_harmonized$CellType <- ML001_and_002_harmonized_celltype$PutativCellType[
#   match(ML001_and_002_harmonized$SCT_snn_res.0.3,
#         ML001_and_002_harmonized_celltype$cluster)
# ]


# =============================================================================
# MW3 – RNA normalization pipeline
# =============================================================================
MW3_seurat_object <- readRDS("MW3_seurat_object_HTO.RDS")
MW3_seurat_object.subset <- subset(MW3_seurat_object,
                                    HTO_classification.global == "Singlet") %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2500, verbose = FALSE) %>%
  ScaleData(features = rownames(.))

MW3_seurat_object.subset <- run_dim_reduction(MW3_seurat_object.subset, dims = 1:14)
MW3_seurat_object.subset <- FindClusters(MW3_seurat_object.subset,
                                          resolution = seq(0.1, 1.0, by = 0.1),
                                          verbose = FALSE)

# Export top markers for manual annotation
Luis_marker <- readxl::read_xlsx("data/Luis_marker.xlsx")
gene_clusters <- Luis_marker %>%
  group_by(gene) %>%
  summarize(clusters = paste(unique(cluster), collapse = ",")) %>%
  ungroup()

MW3_marker <- FindAllMarkers(MW3_seurat_object.subset, only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25)
MW3_marker.top100 <- MW3_marker %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 0.5) %>%
  slice_head(n = 100) %>%
  arrange(cluster, -avg_log2FC) %>%
  ungroup()
MW3_marker.top100$MappedLuisClusterTop200 <-
  gene_clusters$clusters[match(MW3_marker.top100$gene, gene_clusters$gene)]
write.csv(MW3_marker.top100, "output/MW3_marker_res0.5_top100.csv", row.names = FALSE)

# → Load manual annotation from xlsx then assign PutativeCellType
MW3.putative.annotation <- readxl::read_xlsx("data/MW3_marker.res0.5.top100.xlsx", sheet = 2)
MW3_seurat_object.subset@meta.data$PutativeCellType <-
  MW3.putative.annotation$PutativeCellType[
    match(MW3_seurat_object.subset@meta.data$RNA_snn_res.0.5,
          MW3.putative.annotation$cluster)
  ]

# DoubletFinder
MW3_seurat_object.subset <- run_doublet_finder(
  MW3_seurat_object.subset,
  annotations = MW3_seurat_object.subset@meta.data$PutativeCellType,
  use_sct = FALSE
)
MW3_seurat_object.subset@meta.data <- MW3_seurat_object.subset@meta.data %>%
  mutate(cluster = case_when(Ann_final == "Singlet" ~ PutativeCellType, TRUE ~ "Doublet"))

# Add HTO→batch label for downstream use
MW3_seurat_object.subset@meta.data <- MW3_seurat_object.subset@meta.data %>%
  mutate(
    batchID = case_when(
      HTO_classification == "B301" ~ "Ar_WT",
      HTO_classification == "B302" ~ "Ar_WT_Cas",
      HTO_classification == "B303" ~ "Lum+StrAr_KO",
      HTO_classification == "B304" ~ "FVB_WT_cas_3d"
    )
  )

# Extract luminal cells
MW3_luminal <- subset(MW3_seurat_object.subset,
                       cluster %in% c("LumA_LumD", "Lum1_not_sure", "LumP", "LumV"))
saveRDS(MW3_luminal,              "MW3_seurat_object.clean.RDS")
saveRDS(MW3_seurat_object.subset, "MW3_seurat_object.subset.RDS")


# =============================================================================
# MW4 – SCTransform pipeline (with cell-cycle regression)
# =============================================================================
MW4_seurat_object <- readRDS("MW4_seurat_object_HTO.RDS")
MW4_seurat_object.subset <- subset(MW4_seurat_object,
                                    HTO_classification.global == "Singlet")

s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(MW4_seurat_object.subset) <- "RNA"

MW4_seurat_object.subset <- MW4_seurat_object.subset %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData()
MW4_seurat_object.subset <- CellCycleScoring(MW4_seurat_object.subset,
                                              s.features   = s.genes,
                                              g2m.features = g2m.genes,
                                              set.ident    = TRUE, verbose = FALSE)
MW4_seurat_object.subset <- SCTransform(MW4_seurat_object.subset,
                                         vars.to.regress = c("percent.mt", "S.Score",
                                                             "G2M.Score", "HTO_maxID"))
MW4_seurat_object.subset <- RunPCA(MW4_seurat_object.subset, assay = "SCT",
                                    npcs = 50, verbose = FALSE) %>%
  RunUMAP(assay = "SCT", dims = 1:40, seed.use = 1234, verbose = FALSE) %>%
  RunTSNE(assay = "SCT", dims = 1:40, seed.use = 1234, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:40, verbose = FALSE) %>%
  FindClusters(resolution = seq(0.1, 1.0, by = 0.1), verbose = FALSE)

Idents(MW4_seurat_object.subset) <- "SCT_snn_res.0.4"
MW4_seurat_object.subset <- PrepSCTFindMarkers(MW4_seurat_object.subset)
MW4_marker <- FindAllMarkers(MW4_seurat_object.subset, only.pos = TRUE,
                              min.pct = 0.4, logfc.threshold = 0.25,
                              group.by = "SCT_snn_res.0.4")
MW4_marker.top100 <- MW4_marker %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 0.5) %>%
  slice_head(n = 100) %>%
  arrange(cluster, -avg_log2FC) %>%
  ungroup()
MW4_marker.top100$MappedLuisClusterTop200 <-
  gene_clusters$clusters[match(MW4_marker.top100$gene, gene_clusters$gene)]
write.csv(MW4_marker.top100, "output/MW4_marker_res0.5_top100.csv", row.names = FALSE)

MW4_putative.annotation <- readxl::read_xlsx("data/MW4_marker.res0.5.top100.xlsx", sheet = 2)
MW4_seurat_object.subset@meta.data$PutativeCellType <-
  MW4_putative.annotation$PutativeCellType[
    match(MW4_seurat_object.subset@meta.data$SCT_snn_res.0.4,
          MW4_putative.annotation$cluster)
  ]
saveRDS(MW4_seurat_object.subset, "MW4_seurat_object.subset.RDS")

MW4_seurat_object.subset <- run_doublet_finder(
  MW4_seurat_object.subset,
  annotations = MW4_seurat_object.subset@meta.data$PutativeCellType,
  use_sct = TRUE
)
MW4_seurat_object.subset@meta.data <- MW4_seurat_object.subset@meta.data %>%
  mutate(
    cluster = case_when(Ann_final == "Singlet" ~ PutativeCellType, TRUE ~ "Doublet"),
    CellType = case_when(
      HTO_maxID == "B301" ~ "Atg7-WT intact",
      HTO_maxID == "B302" ~ "Atg7-WT intact Castration 4w",
      HTO_maxID == "B303" ~ "Atg7-Ko intact",
      HTO_maxID == "B304" ~ "NKX3.1+/+ 4w",
      HTO_maxID == "B305" ~ "NKX3.1-/- 4w"
    )
  )
MW4_luminal <- subset(MW4_seurat_object.subset,
                       cluster %in% c("Lum_not_sure", "LumA", "LumP"))
saveRDS(MW4_luminal,              "MW4_luminal.RDS")
saveRDS(MW4_seurat_object.subset, "MW4_seurat_object.subset.RDS")


# =============================================================================
# MW8 – Harmony-integrated SCTransform pipeline
# =============================================================================
MW8_obj <- readRDS("MW8_seurat_object.RDS")

# Rebuild a clean Seurat object from raw counts to avoid layer conflicts
MW8_counts        <- GetAssayData(MW8_obj, layer = "counts", assay = "RNA")
MW8_seurat_object <- CreateSeuratObject(counts = MW8_counts,
                                         meta.data = MW8_obj@meta.data, assay = "RNA")
rm(MW8_obj)

# Split by batch → SCTransform → Harmony integration
MW8_seurat_object[["RNA"]] <- split(MW8_seurat_object[["RNA"]],
                                     f = MW8_seurat_object$Batch)
MW8_seurat_object <- MW8_seurat_object %>%
  SCTransform() %>%
  RunPCA(verbose = FALSE)

MW8_seurat_object <- IntegrateLayers(
  object               = MW8_seurat_object,
  method               = HarmonyIntegration,
  new.reduction        = "integrated_harmony_SCT",
  normalization.method = "SCT",
  verbose              = FALSE
)
MW8_seurat_object <- RunUMAP(MW8_seurat_object, dims = 1:30,
                               reduction = "integrated_harmony_SCT", verbose = FALSE)
MW8_seurat_object <- RunTSNE(MW8_seurat_object, dims = 1:30,
                               reduction = "integrated_harmony_SCT", verbose = FALSE)
MW8_seurat_object <- FindNeighbors(MW8_seurat_object,
                                    reduction = "integrated_harmony_SCT") %>%
  FindClusters(resolution = seq(0.1, 1, 0.1), verbose = FALSE)

Idents(MW8_seurat_object) <- "SCT_snn_res.0.4"
MW8_seurat_object <- PrepSCTFindMarkers(MW8_seurat_object)
MW8_marker <- FindAllMarkers(MW8_seurat_object, only.pos = TRUE,
                              min.pct = 0.4, logfc.threshold = 0.25,
                              group.by = "SCT_snn_res.0.4")
write.csv(MW8_marker, "output/MW8_marker_res0.4_top100.csv", row.names = FALSE)

# Manual annotation → load from xlsx
MW8_annotation <- readxl::read_xlsx("data/MW7_marker.res0.4.top100.xlsx", sheet = 2)
MW8_seurat_object$cluster <- MW8_annotation$PutativeCellType[
  match(MW8_seurat_object$SCT_snn_res.0.4, MW8_annotation$cluster)
]

# DoubletFinder
MW8_seurat_object <- run_doublet_finder(
  MW8_seurat_object,
  annotations = MW8_seurat_object@meta.data$cluster,
  use_sct = TRUE
)
MW8_seurat_object <- subset(MW8_seurat_object, Ann_final == "Singlet")
saveRDS(MW8_seurat_object, "MW8_seurat_object.RDS")

# Re-run UMAP/tSNE after doublet removal
MW8_seurat_object <- RunUMAP(MW8_seurat_object, dims = 1:30,
                               reduction = "integrated_harmony_SCT", verbose = FALSE)

# Optional cluster update from manual review CSV
# MW8_cluster_update <- read.csv("data/cluster_update.csv")
# MW8_seurat_object$ClusterUpdate <- MW8_cluster_update$cluster[
#   match(rownames(MW8_seurat_object@meta.data), MW8_cluster_update$Barcode)
# ]
MW8_seurat_object$ClusterUpdate <- gsub("LumA_seperate", "CyclingCell",
                                         MW8_seurat_object$ClusterUpdate)
MW8_luminal <- subset(MW8_seurat_object, ClusterUpdate %in% c("LumA", "LumP"))
saveRDS(MW8_seurat_object, "MW8_seurat_object.RDS")
