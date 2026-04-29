# =============================================================================
# 01_data_loading_hto_demux.R
# Load 10X Genomics Cell-Hashing data and perform HTO-based demultiplexing
# for in-house samples: ML001, ML002, MW3, MW4, MW8
#
# Input:  Raw 10X output directories (filtered_feature_bc_matrix)
# Output: Per-sample Seurat objects with HTO assay and demultiplexing labels
#         saved as .RDS files
# =============================================================================

source("R/00_setup.R")

# Helper: load a 10X CellHash sample, apply QC filters, add HTO assay,
# and run HTODemux. Returns a Seurat object with HTO_classification.global set.
load_and_demux_sample <- function(gex_dir,
                                  sample_name,
                                  mt_cutoff    = 15,
                                  min_features = 200,
                                  hto_assay    = "HTO",
                                  clr_margin   = 2) {

  message("Loading: ", sample_name)

  # 1. Load GEX + HTO count matrices
  raw <- Read10X(data.dir = gex_dir)

  # 2. Create Seurat object from GEX
  seurat_obj <- CreateSeuratObject(counts = raw[["Gene Expression"]])
  seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^mt", col.name = "percent.mt")

  # 3. QC filter
  seurat_obj <- subset(seurat_obj,
                       percent.mt < mt_cutoff & nFeature_RNA > min_features)

  # 4. Normalize RNA
  seurat_obj <- NormalizeData(seurat_obj,
                               normalization.method = "LogNormalize",
                               scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "mean.var.plot")
  seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

  # 5. Prepare HTO matrix
  htos_raw <- as.data.frame(raw[["Antibody Capture"]])
  htos_raw <- htos_raw[, colSums(htos_raw) > 0]           # drop empty barcodes
  joint_bcs <- intersect(colnames(htos_raw), colnames(seurat_obj))

  seurat_obj <- seurat_obj[, joint_bcs]
  htos_mat   <- as.matrix(htos_raw[, joint_bcs])

  # 6. Add HTO assay and normalize (CLR)
  seurat_obj[[hto_assay]] <- CreateAssayObject(counts = htos_mat)
  seurat_obj <- NormalizeData(seurat_obj, assay = hto_assay,
                               normalization.method = "CLR", margin = clr_margin)

  # 7. Demultiplex
  seurat_obj <- HTODemux(seurat_obj, assay = hto_assay, positive.quantile = 0.99)

  message("  Singlets: ",
          sum(seurat_obj$HTO_classification.global == "Singlet"),
          " / ", ncol(seurat_obj), " total cells")
  return(seurat_obj)
}


# =============================================================================
# ML001
# =============================================================================
ML001_seurat_obj <- load_and_demux_sample(
  gex_dir     = "/path/to/ML001/filtered_feature_bc_matrix",
  sample_name = "ML001"
)

ML001_seurat_obj$timepoint <- sapply(ML001_seurat_obj@meta.data$HTO_maxID,
          function(x)
              ifelse(
                  x %in% c("B301","B303"),
                  "1w",
                  "4w"
              )
              
          )
ML001_seurat_obj$CellType <- sapply(ML001_seurat_obj@meta.data$HTO_maxID,
          function(x)
              ifelse(
                  x %in% c("B301","B302"),
                  "Luminal AR KO",
                  "Stromal AR KO"
              )
              
          )

saveRDS(ML001_seurat_obj, "ML001_seurat_obj_HTO.RDS")


# =============================================================================
# ML002
# =============================================================================
ML002_seurat_obj <- load_and_demux_sample(
  gex_dir     = "/path/to/ML002/filtered_feature_bc_matrix",
  sample_name = "ML002"
)
ML002_seurat_obj$timepoint <- sapply(ML002_seurat_obj@meta.data$HTO_maxID,
          function(x)
              ifelse(
                  x %in% c("B301","B303"),
                  "1w",
                  "4w"
              )
              
          )

ML002_seurat_obj$CellType <- sapply(ML002_seurat_obj@meta.data$HTO_maxID,
          function(x)
              ifelse(
                  x %in% c("B301","B302"),
                  "Luminal Acsl4 KO",
                  "Luminal Atg7 KO "
              )
              
          )


saveRDS(ML002_seurat_obj, "ML002_seurat_obj_HTO.RDS")


# =============================================================================
# MW3
# =============================================================================
MW3_seurat_object <- load_and_demux_sample(
  gex_dir     = "/path/to/MW3/filtered_feature_bc_matrix",
  sample_name = "MW3"
)

MW3_seurat_object@meta.data <- 
        MW3_seurat_object@meta.data %>%
             mutate( 
                    CellType=case_when(
                        HTO_maxID=="B301" ~ "AR-intact",
                        HTO_maxID=="B302" ~ "AR Castration 4w",
                        HTO_maxID=="B303" ~ "Stromal and Luminal AR KO",
                        HTO_maxID=="B304" ~ "FVB WT cas 3d"
                                        )
        )
saveRDS(MW3_seurat_object, "MW3_seurat_object_HTO.RDS")


# =============================================================================
# MW4
# =============================================================================
MW4_seurat_object <- load_and_demux_sample(
  gex_dir     = "/path/to/MW4/filtered_feature_bc_matrix",
  sample_name = "MW4"
)
MW4_seurat_object@meta.data <-MW4_seurat_object@meta.data %>%
      mutate(
        CellType = case_when(
            HTO_maxID=="B301" ~ "Atg7-WT intact",
            HTO_maxID=="B302" ~ "Atg7-WT intact Castration 4w",
            HTO_maxID=="B303" ~ "Atg7-Ko intact",
            HTO_maxID=="B304" ~ "NKX3.1+/+ 4w",
            HTO_maxID=="B305" ~ "NKX3.1-/- 4w",
        )
    )
saveRDS(MW4_seurat_object, "MW4_seurat_object_HTO.RDS")


# =============================================================================
# MW8
# =============================================================================
MW8_seurat_object <- load_and_demux_sample(
  gex_dir     = "/path/to/MW8/filtered_feature_bc_matrix",
  sample_name = "MW8"
)

# Annotate HTO barcodes -> experimental batch
MW8_seurat_object@meta.data <- MW8_seurat_object@meta.data %>%
  mutate(
    Batch = case_when(
      HTO_classification == 'B301' ~ 'PTN_KO_4w',
      HTO_classification == 'B302' ~ 'PTN_NKX_DKO_4w',
      HTO_classification == 'B303' ~ 'RAPTOR_KO_4w',
      HTO_classification == 'B304' ~ 'RAPTOR_AR_DKO_4w'
    )
  )

# Keep singlets only (demultiplexing QC)
MW8_seurat_object <- subset(MW8_seurat_object,
                             HTO_classification.global == "Singlet")
table(MW8_seurat_object@meta.data$HTO_maxID)
saveRDS(MW8_seurat_object, "MW8_seurat_object.RDS")


# =============================================================================
# HTO QC visualization helper
# Call on any sample to reproduce the standard HTO QC panel
# =============================================================================
plot_hto_qc <- function(seurat_obj, sample_name, hto_pairs) {
  # Ridge plot – HTO signal per barcode
  Idents(seurat_obj) <- "HTO_maxID"
  p_ridge <- RidgePlot(seurat_obj, assay = "HTO",
                        features = rownames(seurat_obj[["HTO"]])[seq_len(4)],
                        ncol = 2)

  # Feature scatter – mutual exclusivity
  scatter_list <- lapply(hto_pairs, function(pair) {
    FeatureScatter(seurat_obj, feature1 = pair[1], feature2 = pair[2])
  })
  p_scatter <- patchwork::wrap_plots(scatter_list, ncol = 2)

  # UMI violin per global class
  Idents(seurat_obj) <- "HTO_classification.global"
  p_vln <- VlnPlot(seurat_obj, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

  list(ridge = p_ridge, scatter = p_scatter, vln = p_vln)
}
