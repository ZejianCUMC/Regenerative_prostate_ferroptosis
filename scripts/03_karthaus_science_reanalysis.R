# =============================================================================
# 03_karthaus_science_reanalysis.R
# Re-analysis of Karthaus et al. 2020 Science scRNA-seq datasets:
#   - Mouse castration time-course (GSE146811)
#   - Human prostate atlas
#
# This script reproduces key figures from the original paper and generates
# PTN/IGF1 expression visualizations used in manuscript figures.
#
# Input:  Public GEO data + metadata files
# Output: Seurat objects for mouse & human datasets; PDF figures
# =============================================================================

source("R/00_setup.R")
library(biomaRt)

# =============================================================================
# Mouse dataset (Karthaus et al. 2020 Science)
# GSE146811: mmProstate10x_timecourse
# =============================================================================

# --- Load count matrix + metadata -------------------------------------------
Wouter_science_obj <- Read10X_h5("data/GSE146811_mmProstate10x_timecourse_rawCount.h5")
Wouter_science_obj <- CreateSeuratObject(counts = Wouter_science_obj)

Wouter_science_meta <- data.table::fread("data/mmProstate10x_scPortal_metadata_modify.txt")
Wouter_science_meta <- data.frame(Wouter_science_meta)

Wouter_science_obj@meta.data$NAME <- rownames(Wouter_science_obj@meta.data)
Wouter_science_obj@meta.data <- left_join(Wouter_science_obj@meta.data,
                                           Wouter_science_meta)
rownames(Wouter_science_obj@meta.data) <- Wouter_science_obj$NAME
Wouter_science_obj@meta.data$NAME     <- NULL

# Remove predicted doublets flagged in the original metadata
Wouter_science_obj <- subset(Wouter_science_obj, IntType != "PredDoublet")

# Basic normalization
Wouter_science_obj <- Wouter_science_obj %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData()

saveRDS(Wouter_science_obj, "Wouter_science_obj.RDS")


# =============================================================================
# Mouse: Intact time point (T00) – reproduce paper tSNE embedding
# =============================================================================
Wouter_science_obj.T00_intact <- subset(Wouter_science_obj, timepoint == "T00_Intact")

# Batch-aware SCTransform + CCA integration
Wouter_science_obj.T00_intact[["RNA"]] <- split(Wouter_science_obj.T00_intact[["RNA"]],
                                                 f = Wouter_science_obj.T00_intact$batchID)
Wouter_science_obj.T00_intact <- Wouter_science_obj.T00_intact %>% SCTransform()

Wouter_science_obj.T00_intact <- IntegrateLayers(
  object               = Wouter_science_obj.T00_intact,
  method               = CCAIntegration,
  new.reduction        = "integrated_cca_SCT",
  normalization.method = "SCT", verbose = FALSE
)
Wouter_science_obj.T00_intact <- RunUMAP(Wouter_science_obj.T00_intact,
                                          dims = 1:30, reduction = "integrated_cca_SCT")
Wouter_science_obj.T00_intact <- RunTSNE(Wouter_science_obj.T00_intact,
                                          dims = 1:30, reduction = "integrated_cca_SCT")

# Embed paper tSNE coordinates (for direct comparison with published figure)
Wouter_science_T00_tsne <- read.delim("data/mmProstate10x_global_T0_tSNE.csv", sep = ",")
DefaultAssay(Wouter_science_obj.T00_intact) <- "RNA"

umap_df <- data.frame(
  UMAP1 = as.numeric(Wouter_science_T00_tsne$X),
  UMAP2 = as.numeric(Wouter_science_T00_tsne$Y),
  row.names = Wouter_science_T00_tsne$NAME
)
umap_df <- na.omit(umap_df[rownames(Wouter_science_obj.T00_intact@meta.data), ])

Wouter_science_obj.T00_intact[["TSNE_paper"]] <- CreateDimReducObject(
  embeddings = as.matrix(umap_df), key = "TSNE_", assay = "RNA"
)

# --- Ext data Fig8k: Ptn expression at Intact timepoint ----------------------
pT00  <- DimPlot(Wouter_science_obj.T00_intact, reduction = "TSNE_paper",
                  group.by = "FullType", cols = col.cluster2.1, label = TRUE,
                  label.size = 5) + theme(legend.position = "bottom")
pPtn  <- FeaturePlot(Wouter_science_obj.T00_intact, reduction = "TSNE_paper",
                      features = "Ptn")
pIgf1 <- FeaturePlot(Wouter_science_obj.T00_intact, reduction = "TSNE_paper",
                      features = "Igf1")

ggsave("output/Wouter_science_T00_tSNE.pdf",       pT00,  width = 12, height = 12, dpi = 300)
ggsave("output/Wouter_science_T00_tSNE_Ptn.pdf",   pPtn,  width = 5,  height = 5,  dpi = 300)
ggsave("output/Wouter_science_T00_tSNE_Igf1.pdf",  pIgf1, width = 5,  height = 5,  dpi = 300)

# Ptn-positive cell fraction by cell type
Wouter_science_obj.T00_intact$Ptn_positive <-
  Wouter_science_obj.T00_intact[["RNA"]]@data["Ptn", ] >= 1
Ptn_percentage <- Wouter_science_obj.T00_intact@meta.data %>%
  group_by(FullType) %>%
  summarise(Ptn_positive_percent = mean(Ptn_positive) * 100)

saveRDS(Wouter_science_obj.T00_intact, "Wouter_science_obj.T00_intact.RDS")


# =============================================================================
# Mouse: Extract luminal cells (all time points)
# =============================================================================
Wouter_science_obj     <- readRDS("Wouter_science_obj.RDS")
Wouter_science_luminal <- subset(Wouter_science_obj, IntType == "Epi_Luminal")
Wouter_science_luminal <- subset(Wouter_science_luminal,
                                  !timepoint %in% c("T00_Epi", "T00_NonEpi"))
Wouter_science_luminal <- subset(Wouter_science_luminal,
                                  !HighLevel %in% c("Imm", "Str"))

Wouter_science_luminal <- Wouter_science_luminal %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(features = rownames(.)) %>%
  RunPCA(verbose = FALSE) %>%
  RunTSNE(reduction = "pca", dims = 1:30, verbose = FALSE, seed.use = 1234) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE, seed.use = 1234) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE)

saveRDS(Wouter_science_luminal, "Wouter_science_luminal_seurat_obj.RDS")

# Export normalized matrix + metadata for external classifier
normalized_counts <- GetAssayData(Wouter_science_luminal, slot = "data", assay = "RNA")
data.table::fwrite(as.data.frame(as.matrix(normalized_counts)),
                   file = "output/Wouter_science_luminal_normalized.csv", row.names = TRUE)
write.csv(Wouter_science_luminal@meta.data,
          "output/Wouter_science_luminal_metainfo.csv", row.names = TRUE)


# =============================================================================
# Mouse: Stromal subset → Ptn / Igf1 expression analysis
# =============================================================================
Wouter_science_obj     <- readRDS("Wouter_science_obj.RDS")
Wouter_science_stromal <- subset(Wouter_science_obj,
                                  HighLevel == "Str" &
                                  !timepoint %in% c("T00_Epi", "T00_NonEpi") &
                                  IntType != "PredDoublet")
rm(Wouter_science_obj)

# Remove residual epithelial cells from stromal subset
Wouter_science_stromal <- subset(Wouter_science_stromal,
                                  !IntType %in% c("Epi_Basal", "Epi_Luminal", "Epi_SV"))

# Batch-integrated normalization (CCA)
Wouter_science_stromal[["RNA"]] <- split(Wouter_science_stromal[["RNA"]],
                                          f = Wouter_science_stromal$batchID)
Wouter_science_stromal <- Wouter_science_stromal %>% SCTransform()

Wouter_science_stromal <- IntegrateLayers(
  object               = Wouter_science_stromal,
  method               = CCAIntegration,
  new.reduction        = "integrated_cca_SCT",
  normalization.method = "SCT", verbose = FALSE
)
Wouter_science_stromal <- RunUMAP(Wouter_science_stromal, dims = 1:30,
                                   reduction = "integrated_cca_SCT")
Wouter_science_stromal <- FindNeighbors(Wouter_science_stromal,
                                         reduction = "integrated_cca_SCT", dims = 1:30) %>%
  FindClusters(resolution = seq(0.1, 1, 0.1), verbose = FALSE)

saveRDS(Wouter_science_stromal, "Wouter_science_stromal.RDS")

# Ptn-positive cluster comparison: Intact vs Cast Day7
Wouter_science_stromal_PtnPos <- subset(Wouter_science_stromal,
                                         !SCT_snn_res.0.3 %in% c("1", "2", "8", "12", "14"))
Wouter_science_stromal_PtnPos_Intact_1w <- subset(Wouter_science_stromal_PtnPos,
                                                    timepoint %in% c("T00_Intact", "T02_Cast_Day7"))
Wouter_science_stromal_PtnPos_Intact_1w <- PrepSCTFindMarkers(
  Wouter_science_stromal_PtnPos_Intact_1w
)
Ptn_comparison <- FindMarkers(Wouter_science_stromal_PtnPos_Intact_1w,
                               ident.1 = "T02_Cast_Day7", ident.2 = "T00_Intact",
                               group.by = "timepoint", min.pct = 0.2, logfc.threshold = 0.2)

# --- Ext Data Fig8l: Stromal Ptn expression ----------------------------------
Idents(Wouter_science_stromal_PtnPos_Intact_1w) <- "timepoint"
p_ptn_stromal <- VlnPlot(Wouter_science_stromal_PtnPos_Intact_1w,
                           group.by = "timepoint", assay = "SCT", features = "Ptn",
                           cols = c("#cf385e", "#468dc1"), pt.size = 0, log = FALSE) +
  labs(y = "Normalized Counts", title = "Ptn Positive Stromal") +
  theme(legend.position = "none")
ggsave("output/Ptn_at_Ptn_positive_stromal.pdf", p_ptn_stromal, width = 4, height = 5)

# Igf1-positive cluster (Ext data Fig8b)
Wouter_science_stromal_Igf1Pos <- subset(Wouter_science_stromal,
                                           !SCT_snn_res.0.3 %in% c("2", "12", "14"))
Wouter_science_stromal_Igf1Pos_Intact_1w <- subset(Wouter_science_stromal_Igf1Pos,
                                                     timepoint %in% c("T00_Intact", "T02_Cast_Day7"))
Wouter_science_stromal_Igf1Pos_Intact_1w <- PrepSCTFindMarkers(
  Wouter_science_stromal_Igf1Pos_Intact_1w
)
Igf1_comparison <- FindMarkers(Wouter_science_stromal_Igf1Pos_Intact_1w,
                                ident.1 = "T02_Cast_Day7", ident.2 = "T00_Intact",
                                group.by = "timepoint", min.pct = 0.2, logfc.threshold = 0.2)

Idents(Wouter_science_stromal_Igf1Pos_Intact_1w) <- "timepoint"
p2 <- VlnPlot(Wouter_science_stromal_Igf1Pos_Intact_1w,
               group.by = "timepoint", assay = "SCT", features = "Igf1",
               cols = c("#cf385e", "#468dc1"), pt.size = 0, log = FALSE) +
  labs(y = "Normalized Counts", title = "Igf1 Positive Stromal") +
  theme(legend.position = "none")
ggsave("output/Igf1_at_Igf1_positive_stromal.pdf", p2, width = 5, height = 5)


# =============================================================================
# Human dataset (Karthaus et al. 2020 Science)
# =============================================================================
Wouter_science_human <- Read10X(
  data.dir = "data/20Science_human/"
)
Wouter_science_human <- CreateSeuratObject(Wouter_science_human)

Wouter_science_human_meta <- data.table::fread("data/Wounter_Science_hsPortalData_metadata.txt")
Wouter_science_human@meta.data$NAME <- rownames(Wouter_science_human@meta.data)
Wouter_science_human@meta.data <- left_join(Wouter_science_human@meta.data,
                                              data.frame(Wouter_science_human_meta))
rownames(Wouter_science_human@meta.data) <- Wouter_science_human@meta.data$NAME
Wouter_science_human@meta.data$NAME <- NULL

# Batch integration across patients
Wouter_science_human[["RNA"]] <- split(Wouter_science_human[["RNA"]],
                                        f = Wouter_science_human$batchID)
Wouter_science_human <- Wouter_science_human %>%
  FindVariableFeatures() %>% ScaleData() %>% SCTransform() %>% RunPCA()

Wouter_science_human <- IntegrateLayers(
  object               = Wouter_science_human,
  method               = HarmonyIntegration,
  new.reduction        = "integrated_harmony_SCT",
  normalization.method = "SCT", verbose = FALSE
)
Wouter_science_human <- RunUMAP(Wouter_science_human, dims = 1:30,
                                 reduction = "integrated_harmony_SCT", verbose = FALSE)
Wouter_science_human <- RunTSNE(Wouter_science_human, dims = 1:30,
                                 reduction = "integrated_harmony_SCT", verbose = FALSE)
Wouter_science_human <- FindNeighbors(Wouter_science_human,
                                       reduction = "integrated_harmony_SCT") %>%
  FindClusters(resolution = seq(0.1, 1, 0.1), verbose = FALSE)

saveRDS(Wouter_science_human, "Wouter_science_human_reanalysis.RDS")

# Marker genes + annotation (see Methods for Luis_marker reference set)
Wouter_science_human <- PrepSCTFindMarkers(Wouter_science_human)
Wouter_science_human_marker <- FindAllMarkers(Wouter_science_human, only.pos = TRUE,
                                               min.pct = 0.4, logfc.threshold = 0.25,
                                               group.by = "SCT_snn_res.0.7")

# Mouse → human ortholog mapping via biomaRt
human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                       host = "https://dec2021.archive.ensembl.org/")
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                       host = "https://dec2021.archive.ensembl.org/")

Luis_marker <- readxl::read_xlsx("data/Luis_marker.xlsx")
gene_clusters_human <- Luis_marker %>%
  group_by(gene) %>%
  summarize(clusters = paste(unique(cluster), collapse = ",")) %>%
  ungroup()

genesV2 <- getLDS(attributes  = "mgi_symbol", filters = "mgi_symbol",
                   values      = gene_clusters_human$gene, mart = mouse_mart,
                   attributesL = "hgnc_symbol", martL = human_mart, uniqueRows = TRUE)
gene_clusters_human$HGNC.symbol <-
  genesV2$HGNC.symbol[match(gene_clusters_human$gene, genesV2$MGI.symbol)]

Wouter_science_human_marker.top100 <- Wouter_science_human_marker %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 0.5) %>%
  slice_head(n = 100) %>%
  arrange(cluster, -avg_log2FC) %>%
  ungroup()
Wouter_science_human_marker.top100$MappedLuisClusterTop200 <-
  gene_clusters_human$clusters[
    match(Wouter_science_human_marker.top100$gene, gene_clusters_human$HGNC.symbol)
  ]
write.csv(Wouter_science_human_marker.top100,
          "output/Wouter_science_human_marker.csv", row.names = FALSE)

# Load manual cell-type annotation
Wouter_science_human_celltype <- readxl::read_xlsx("data/Wouter_science_human_marker.xlsx",
                                                     sheet = 2)
Wouter_science_human$PredCellType <- Wouter_science_human_celltype$PutativeCellType[
  match(Wouter_science_human$SCT_snn_res.0.7,
        Wouter_science_human_celltype$cluster)
]
Wouter_science_human@meta.data <- Wouter_science_human@meta.data %>%
  mutate(PredCellType_v2 = case_when(
    PredCellType %in% c("Basal", "Luminal") ~ "Epithelial",
    TRUE ~ PredCellType
  ))

# PTN & IGF1 feature plots (for manuscript figures)
p1 <- DimPlot(Wouter_science_human, reduction = "umap", group.by = "PredCellType_v2",
              cols = col.cluster2.1, label = TRUE, label.size = 5) +
  theme(legend.position = "right")
p2_ptn  <- FeaturePlot(Wouter_science_human, reduction = "umap",
                        features = "PTN",  max.cutoff = 1.5)
p2_igf1 <- FeaturePlot(Wouter_science_human, reduction = "umap",
                        features = "IGF1", max.cutoff = 1.5)

pdf("output/Wouter_science_human_intact_umap_with_PTN_expression.pdf",
    width = 12, height = 7)
p1 + p2_ptn
dev.off()

pdf("output/Wouter_science_human_intact_umap_with_IGF1_expression.pdf",
    width = 12, height = 7)
p1 + p2_igf1
dev.off()

saveRDS(Wouter_science_human, "Wouter_science_human_reanalysis.RDS")


# =============================================================================
# Human: Extract epithelial cells + embed published tSNE coordinates
# =============================================================================
Wouter_science_human_epithelial <- subset(Wouter_science_human,
                                           PredCellType %in% c("Basal", "Luminal1", "Luminal2"))
hs_epi_tsne <- read.csv("data/hsProst10x_postCNV_Epi_qTPM_tSNE.txt")

Wouter_science_human_epithelial@meta.data$NAME <- rownames(Wouter_science_human_epithelial@meta.data)
Wouter_science_human_epithelial_normal <- subset(Wouter_science_human_epithelial,
                                                  NAME %in% hs_epi_tsne$NAME)
Wouter_science_human_epithelial_normal@meta.data <-
  left_join(Wouter_science_human_epithelial_normal@meta.data, hs_epi_tsne[, -c(2, 3)])
rownames(Wouter_science_human_epithelial_normal@meta.data) <-
  Wouter_science_human_epithelial_normal@meta.data$NAME
Wouter_science_human_epithelial_normal@meta.data$NAME <- NULL

# Embed tSNE coordinates from paper
tsne_df <- data.frame(
  UMAP1 = as.numeric(hs_epi_tsne$X),
  UMAP2 = as.numeric(hs_epi_tsne$Y),
  row.names = hs_epi_tsne$NAME
)
tsne_df <- na.omit(tsne_df[rownames(Wouter_science_human_epithelial_normal@meta.data), ])
Wouter_science_human_epithelial_normal[["TSNE_paper"]] <- CreateDimReducObject(
  embeddings = as.matrix(tsne_df), key = "TSNE_", assay = "RNA"
)
Wouter_science_human_epithelial_NT <- subset(Wouter_science_human_epithelial_normal,
                                              isTumorRegion == "NT")
saveRDS(Wouter_science_human_epithelial_normal, "Wouter_science_human_epithelial.RDS")

# Differential expression: ADT vs INTACT in Luminal 1 (Acinar) and Luminal 2 (Ductal)
Wouter_science_human_lum1 <- subset(Wouter_science_human_epithelial_NT,
                                     FullTypePred == "Epi_Luminal_1" & isTumorRegion == "NT")
Idents(Wouter_science_human_lum1) <- "treatment"
lum1_diff <- FindMarkers(Wouter_science_human_lum1,
                          ident.1 = "ADT", ident.2 = "INTACT",
                          group.by = "treatment", min.pct = 0.1, logfc.threshold = 0)

Wouter_science_human_lum2 <- subset(Wouter_science_human_epithelial_NT,
                                     FullTypePred == "Epi_Luminal_2" & isTumorRegion == "NT")
Idents(Wouter_science_human_lum2) <- "treatment"
Wouter_science_human_lum2$treatment <- factor(Wouter_science_human_lum2$treatment,
                                               levels = c("INTACT", "ADT"))
lum2_diff <- FindMarkers(Wouter_science_human_lum2,
                          ident.1 = "ADT", ident.2 = "INTACT",
                          group.by = "treatment", min.pct = 0.1, logfc.threshold = 0)

# --- ITGA3 comparison violin (Ext data figure) --------------------------------
ITGA3_lum1 <- FetchData(Wouter_science_human_lum1, vars = "ITGA3")
ITGA3_lum2 <- FetchData(Wouter_science_human_lum2, vars = "ITGA3")

ITGA3_df <- rbind(
  data.frame(Cell = rownames(ITGA3_lum1), Exp = ITGA3_lum1$ITGA3,
             CellType = "Luminal Acinar\np=2.56e-34"),
  data.frame(Cell = rownames(ITGA3_lum2), Exp = ITGA3_lum2$ITGA3,
             CellType = "Luminal Ductal\np=3.23e-10")
)
ITGA3_df$Y         <- 2^(ITGA3_df$Exp) - 1
ITGA3_df$treatment <- Wouter_science_human_epithelial_normal@meta.data$treatment[
  match(ITGA3_df$Cell, rownames(Wouter_science_human_epithelial_normal@meta.data))
]

p_itga3 <- ggviolin(ITGA3_df, trim = TRUE, width = 1,
                     x = "treatment", y = "Exp", fill = "treatment") +
  geom_boxplot(width = 0.1, fill = "white") +
  facet_wrap(~CellType, scales = "free") +
  labs(title = "ITGA3 (Karthaus, et al. 2020)",
       y = "log2 Normalized Counts") +
  scale_fill_manual(values  = c("#cf385e", "#468dc1")) +
  scale_color_manual(values = c("#cf385e", "#468dc1")) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text       = element_text(size = 14, face = "bold"))

pdf("output/Wouter_science_human_ITGA3_violin_luminal_merged.pdf", width = 9, height = 5)
p_itga3
dev.off()
