# =============================================================================
# 07_stromal_mesenchymal_analysis.R
# Stromal cell analysis: Mesenchymal + Smooth Muscle subpopulations
# with Ptn/Igf1/Ar expression dynamics across castration time points
#
# Data: Karthaus et al. 2020 Science (Wouter cohort)
#
# Input:  Wouter_science_obj.RDS
# Output: Mesenchymal Seurat objects, Ptn expression plots, PDF figures
# =============================================================================

source("R/00_setup.R")

CASTRATION_TIMEPOINTS <- c("T00_Intact", "T02_Cast_Day7",
                             "T03_Cast_Day14", "T04_Cast_Day28")

# =============================================================================
# Helper: standard SCTransform normalization + dimensional reduction
# =============================================================================
normalize_and_reduce_sct <- function(obj, nfeatures = 3000, dims = 1:40,
                                      seed = 1234, vars_regress = NULL) {
  obj <- obj %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = nfeatures) %>%
    ScaleData() %>%
    SCTransform(vars.to.regress = vars_regress)
  obj <- RunPCA(obj, assay = "SCT", npcs = 50, verbose = FALSE)
  obj <- RunUMAP(obj, assay = "SCT", dims = dims, seed.use = seed, verbose = FALSE)
  obj <- RunTSNE(obj, assay = "SCT", dims = dims, seed.use = seed, verbose = FALSE)
  return(obj)
}


# =============================================================================
# 1. Mesenchymal cells – per-timepoint analysis
# =============================================================================
Wouter_science_obj          <- readRDS("Wouter_science_obj.RDS")
Wouter_science_mesenchymal  <- subset(Wouter_science_obj, IntType == "Str_Mesenchymal")
saveRDS(Wouter_science_mesenchymal, "Wouter_science_mesenchymal.RDS")
rm(Wouter_science_obj)

# Per-timepoint UMAPs
tp_labels <- c(T00_Intact      = "T00_Intact",
               T01_Cast_Day1   = "T01_Cast_Day1",
               T02_Cast_Day7   = "T02_Cast_Day7",
               T03_Cast_Day14  = "T03_Cast_Day14",
               T04_Cast_Day28  = "T04_Cast_Day28")

tp_colors <- c(Str_Mesenchymal_1 = "#12b97b", Str_Mesenchymal_2 = "#3cb1f6")

umap_plots    <- list()
feature_plots <- list()

for (tp in names(tp_labels)) {
  sub <- subset(Wouter_science_mesenchymal, timepoint == tp)
  sub <- normalize_and_reduce_sct(sub)

  p1 <- DimPlot(sub, reduction = "umap", group.by = "FullType",
                 cols = tp_colors, label = FALSE) +
    labs(title = paste("Mesenchymal cells at", tp)) +
    theme_classic() +
    theme(legend.position = "bottom")

  p2 <- FeaturePlot(sub, reduction = "umap", features = "Ptn", max.cutoff = 3) +
    theme(legend.position = "right")

  umap_plots[[tp]]    <- p1
  feature_plots[[tp]] <- p2
}

pdf("output/Wouter_Mes_umap_per_timepoint.pdf", width = 10, height = 17)
patchwork::wrap_plots(c(umap_plots, feature_plots), ncol = 2)
dev.off()


# =============================================================================
# 2. Merged mesenchymal analysis: Intact + 1w + 2w + 4w (no batch correction)
# =============================================================================
Wouter_mes_cas <- subset(Wouter_science_mesenchymal,
                          timepoint %in% CASTRATION_TIMEPOINTS)
Wouter_mes_cas <- normalize_and_reduce_sct(Wouter_mes_cas)
saveRDS(Wouter_mes_cas, "Wouter_science_mesenchymal_Cas.RDS")

# Classify cells by Ar expression
Wouter_mes_cas$Ar_positive <- GetAssayData(Wouter_mes_cas,
                                             assay = "RNA", layer = "data")["Ar", ] >= 1
Wouter_mes_cas@meta.data <- Wouter_mes_cas@meta.data %>%
  mutate(Ar_positive = case_when(Ar_positive ~ "AR_positive", TRUE ~ "AR_negative"))

# Ptn violin plots across timepoints (all / AR+ / AR-)
plot_ptn_vln <- function(obj, title_suffix = "") {
  Idents(obj) <- "timepoint"
  VlnPlot(obj, features = "Ptn", pt.size = 0, log = FALSE) +
    labs(title = paste("Ptn expression in Mesenchymal cells", title_suffix,
                        "\n(Karthaus et al., Science 2020)")) +
    theme(legend.position = "none",
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 15, hjust = 0.5))
}
ggsave("output/Ptn_all_Mes_all_TP.pdf",
       plot_ptn_vln(Wouter_mes_cas, "(all)"), width = 12, height = 6)
ggsave("output/Ptn_ARpos_Mes_all_TP.pdf",
       plot_ptn_vln(subset(Wouter_mes_cas, Ar_positive == "AR_positive"), "(AR+)"),
       width = 12, height = 6)
ggsave("output/Ptn_ARneg_Mes_all_TP.pdf",
       plot_ptn_vln(subset(Wouter_mes_cas, Ar_positive == "AR_negative"), "(AR-)"),
       width = 12, height = 6)


# =============================================================================
# 3. Mesenchymal subtype 1 – Ptn expression dynamics with statistical tests
# =============================================================================
Wouter_mes_subtype1 <- subset(Wouter_mes_cas, FullType == "Str_Mesenchymal_1")
Wouter_mes_subtype1 <- PrepSCTFindMarkers(Wouter_mes_subtype1)
Idents(Wouter_mes_subtype1) <- "timepoint"

# Compute DEG for each timepoint vs Intact
for (tp in c("T01_Cast_Day1", "T02_Cast_Day7", "T03_Cast_Day14", "T04_Cast_Day28")) {
  df <- FindMarkers(Wouter_mes_subtype1,
                    ident.1 = tp, ident.2 = "T00_Intact",
                    min.pct = 0.2, logfc.threshold = 0.2)
  if ("Ptn" %in% rownames(df)) {
    message(tp, " - Ptn adj.p = ", df["Ptn", "p_val_adj"])
  }
}

# Statistical comparison violin (Intact, 1w, 2w)
Wouter_mes_sub_filter <- subset(Wouter_mes_subtype1,
                                 !timepoint %in% c("T01_Cast_Day1", "T04_Cast_Day28"))
my_comparisons <- list(c("T00_Intact", "T02_Cast_Day7"),
                        c("T00_Intact", "T03_Cast_Day14"))

pc2 <- VlnPlot(Wouter_mes_sub_filter, features = "Ptn", pt.size = 0, log = FALSE) +
  labs(title = "Ptn expression in subtype1 Mesenchymal cells\n(Karthaus et al., Science 2020)") +
  stat_compare_means(comparisons = my_comparisons) +
  theme(legend.position = "none",
        axis.title = element_text(size = 20, hjust = 0.5),
        plot.title = element_text(size = 15, hjust = 0.5)) +
  scale_y_continuous(limits = c(-0.1, 4.5))

ggsave("output/Ptn_in_Mes1_Intact_1w_2w.pdf", pc2, width = 8, height = 6, dpi = 300)


# =============================================================================
# 4. Mesenchymal + Smooth Muscle combined: Intact vs 4w
# =============================================================================
Wouter_science_obj <- readRDS("Wouter_science_obj.RDS")

Wouter_SmoothMuscle_MES <- subset(Wouter_science_obj,
                                   FullTypeMerged %in% c("Str_SmoothMuscle", "Str_Mesenchymal") &
                                   timepoint %in% CASTRATION_TIMEPOINTS)
Wouter_SmoothMuscle_MES_Glial <- subset(Wouter_science_obj,
                                          FullTypeMerged %in% c("Str_SmoothMuscle",
                                                                  "Str_Mesenchymal",
                                                                  "Str_Glial") &
                                          timepoint %in% CASTRATION_TIMEPOINTS)
rm(Wouter_science_obj)

Wouter_SmoothMuscle_MES <- normalize_and_reduce_sct(Wouter_SmoothMuscle_MES)
Wouter_SmoothMuscle_MES <- FindNeighbors(Wouter_SmoothMuscle_MES) %>%
  FindClusters(resolution = seq(0.1, 0.5, 0.1), verbose = FALSE)

saveRDS(Wouter_SmoothMuscle_MES_Glial, "Wouter_SmoothMuscle_MES_Glial.RDS")

# Subset to Intact vs 4w for publication figure
Wouter_mes_SM_selected <- subset(Wouter_SmoothMuscle_MES,
                                  IntType %in% c("Str_Mesenchymal", "Str_SmoothMuscle") &
                                  timepoint %in% c("T00_Intact", "T04_Cast_Day28"))

Wouter_mes_SM_selected@meta.data <- Wouter_mes_SM_selected@meta.data %>%
  mutate(
    SampleType = case_when(
      IntType == "Str_Mesenchymal" & timepoint == "T00_Intact"     ~ "Mes_Intact",
      IntType == "Str_Mesenchymal" & timepoint == "T04_Cast_Day28" ~ "Mes_D28",
      IntType == "Str_SmoothMuscle" & timepoint == "T00_Intact"    ~ "SmoothMuscle_Intact",
      IntType == "Str_SmoothMuscle" & timepoint == "T04_Cast_Day28"~ "SmoothMuscle_D28"
    ),
    SampleType_v2 = case_when(
      FullType == "Str_Mesenchymal_1"   & timepoint == "T00_Intact"     ~ "Mes_1_Intact",
      FullType == "Str_Mesenchymal_1"   & timepoint == "T04_Cast_Day28" ~ "Mes_1_D28",
      FullType == "Str_Mesenchymal_2"   & timepoint == "T00_Intact"     ~ "Mes_2_Intact",
      FullType == "Str_Mesenchymal_2"   & timepoint == "T04_Cast_Day28" ~ "Mes_2_D28",
      FullType == "Str_SmoothMuscle_1"  & timepoint == "T00_Intact"     ~ "SmoothMuscle_1_Intact",
      FullType == "Str_SmoothMuscle_1"  & timepoint == "T04_Cast_Day28" ~ "SmoothMuscle_1_D28",
      FullType == "Str_SmoothMuscle_2"  & timepoint == "T00_Intact"     ~ "SmoothMuscle_2_Intact",
      FullType == "Str_SmoothMuscle_2"  & timepoint == "T04_Cast_Day28" ~ "SmoothMuscle_2_D28"
    )
  )
Wouter_mes_SM_selected@meta.data$SampleType <- factor(
  Wouter_mes_SM_selected@meta.data$SampleType,
  levels = c("Mes_Intact", "Mes_D28", "SmoothMuscle_Intact", "SmoothMuscle_D28")
)

# AR stratification
Wouter_mes_SM_selected$Ar_positive <-
  GetAssayData(Wouter_mes_SM_selected, assay = "RNA", layer = "data")["Ar", ] >= 1
Wouter_mes_SM_selected@meta.data <- Wouter_mes_SM_selected@meta.data %>%
  mutate(Ar_positive = case_when(Ar_positive ~ "AR_positive", TRUE ~ "AR_negative"))

saveRDS(Wouter_mes_SM_selected, "Wouter_mesenchymal_SmoothMuscle_selected.RDS")

# Ptn violin – all, AR+, AR-
Idents(Wouter_mes_SM_selected) <- "SampleType_v2"

vln_args <- list(features = "Ptn", pt.size = 0, log = FALSE)

p_all <- do.call(VlnPlot, c(list(object = Wouter_mes_SM_selected), vln_args)) +
  labs(title = "Ptn – Mesenchymal & Smooth Muscle (Karthaus 2020)") +
  theme(legend.position = "none")

p_ARpos <- do.call(VlnPlot, c(list(object = subset(Wouter_mes_SM_selected,
                                                     Ar_positive == "AR_positive")), vln_args)) +
  labs(title = "Ptn – AR+ Cells") + theme(legend.position = "none")

p_ARneg <- do.call(VlnPlot, c(list(object = subset(Wouter_mes_SM_selected,
                                                     Ar_positive == "AR_negative")), vln_args)) +
  labs(title = "Ptn – AR- Cells") + theme(legend.position = "none")

ggsave("output/Ptn_all_Mes_SmoothMuscle_0d_28d.pdf",  p_all,   width = 12, height = 6)
ggsave("output/Ptn_ARpos_Mes_SmoothMuscle_0d_28d.pdf", p_ARpos, width = 12, height = 6)
ggsave("output/Ptn_ARneg_Mes_SmoothMuscle_0d_28d.pdf", p_ARneg, width = 12, height = 6)

message("Stromal/mesenchymal analysis complete. Outputs saved to output/")
