# =============================================================================
# 06_gsea_cell_death.R
# Gene Set Enrichment Analysis (GSEA) for cell-death pathway signatures:
#   - Mouse luminal cells: Intact vs Castration Day7 (Wouter cohort)
#   - Human mCRPC: Enzalutamide responders vs non-responders
#   - Human CSPC vs CRPC tumor cells (GSE264573)
#
# Input:  Luminal Seurat objects, pre-ranked gene lists, GMT pathway files
# Output: GSEA result data frames and publication-quality PDF figures
# =============================================================================

source("R/00_setup.R")
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)


# =============================================================================
# Helper: rank gene list by avg_log2FC from FindMarkers output
# =============================================================================
rank_gene_list <- function(markers_df) {
  gene_list        <- markers_df$avg_log2FC
  names(gene_list) <- rownames(markers_df)
  gene_list        <- na.omit(gene_list)
  sort(gene_list, decreasing = TRUE)
}


# =============================================================================
# 1. Mouse luminal cells – Intact vs Cast Day7 (Ext Fig 2g)
# =============================================================================
Wouter_science_luminal <- readRDS("Wouter_science_luminal.filter.RDS")
Cell_death_signature   <- read.gmt("data/Pathway/Cell_death_signature_mouse_update.gmt")

# Subset and run DEG
Wouter_lum_0d_7d <- subset(Wouter_science_luminal,
                             timepoint %in% c("T00_Intact", "T02_Cast_Day7"))
DefaultAssay(Wouter_lum_0d_7d) <- "RNA"
Idents(Wouter_lum_0d_7d)       <- "timepoint"

Wouter_luminal_all_diff <- FindMarkers(Wouter_lum_0d_7d,
                                        ident.1 = "T02_Cast_Day7", ident.2 = "T00_Intact",
                                        group.by = "timepoint", min.pct = 0.2,
                                        logfc.threshold = 0.25)
write.csv(Wouter_luminal_all_diff, "output/Wouter_luminal_7d_vs_0d_DEG.csv")
saveRDS(Wouter_luminal_all_diff,   "Wouter_science_luminal_all_diff.RDS")

gene_list_luminal <- rank_gene_list(Wouter_luminal_all_diff)

Wouter_luminal_gsea <- GSEA(gene_list_luminal,
                              TERM2GENE      = Cell_death_signature,
                              minGSSize      = 5,
                              pvalueCutoff   = 1)

# Ridge plot (Ext Fig 2g)
p_ridge <- ridgeplot(Wouter_luminal_gsea) +
  scale_fill_gradientn(colours = c("red", "orange", "gray40"),
                        breaks  = c(0, 0.025, 0.05),
                        limits  = c(0, 0.05),
                        labels  = c("0", "0.025", ">0.05"),
                        oob     = scales::squish) +
  labs(x = "Enrichment distribution", fill = "Adjusted p-value") +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.title      = element_text(size = 22, hjust = 0.5),
    axis.text       = element_text(size = 20, hjust = 0.5, colour = "black"),
    axis.ticks.y    = element_blank()
  )
ggsave("output/Wouter_luminal_7d_vs_0d_cell_death_ridgeplot.pdf",
       p_ridge, width = 14, height = 4, dpi = 300)


# =============================================================================
# 2. Pairwise Ptn comparisons across time points in mesenchymal Ptn+ cells
# =============================================================================
Wouter_science_stromal   <- readRDS("Wouter_science_stromal.RDS")
Wouter_science_mes       <- subset(Wouter_science_stromal, SCT_snn_res.0.1 %in% c(0, 1, 7))

Wouter_science_mes@meta.data <- Wouter_science_mes@meta.data %>%
  mutate(Group = case_when(SCT_snn_res.0.2 == "0" ~ "PtnNeg", TRUE ~ "PtnPos"))

Wouter_science_mes_PtnPos <- subset(Wouter_science_mes, Group == "PtnPos")
Wouter_science_mes_PtnPos <- PrepSCTFindMarkers(Wouter_science_mes_PtnPos)

timepoints <- c("T00_Intact", "T02_Cast_Day7", "T03_Cast_Day14", "T04_Cast_Day28")

pairwise_results <- lapply(
  combn(timepoints, 2, simplify = FALSE),
  function(pair) {
    tp1 <- pair[1]; tp2 <- pair[2]
    sub_obj <- subset(Wouter_science_mes_PtnPos, timepoint %in% pair)
    diff_df  <- FindMarkers(sub_obj, ident.1 = tp1, ident.2 = tp2,
                             group.by = "timepoint", min.pct = 0.2,
                             logfc.threshold = 0, recorrect_umi = FALSE)
    diff_df$comparison <- paste0(tp1, "_vs_", tp2)
    diff_df$gene       <- rownames(diff_df)
    diff_df %>% filter(gene == "Ptn")
  }
)

pairwise_ptn <- do.call(rbind, pairwise_results)
pairwise_ptn$p_val_adj <- p.adjust(pairwise_ptn$p_val, method = "BH")
write.csv(pairwise_ptn, "output/Ptn_pairwise_timepoint_comparisons.csv", row.names = FALSE)


# =============================================================================
# 3. Stromal/Mesenchymal Ptn p-values for all timepoints (with AR stratification)
# =============================================================================

#' For a given Seurat object, return Ptn adj. p-value between two timepoints
get_ptn_padj <- function(seurat_obj, ident_1, ident_2) {
  Idents(seurat_obj) <- "timepoint"
  df <- FindMarkers(seurat_obj, ident.1 = ident_1, ident.2 = ident_2,
                    min.pct = 0.2, logfc.threshold = 0.2)
  if ("Ptn" %in% rownames(df)) df["Ptn", "p_val_adj"] else NA
}

#' Compute Ptn p-values across all timepoints for each cell type and AR status
get_ptn_pvalues_all_tp <- function(input_obj,
                                    timepoints         = c("T02_Cast_Day7", "T03_Cast_Day14", "T04_Cast_Day28"),
                                    reference_timepoint = "T00_Intact",
                                    cell_types          = c("Str_Mesenchymal", "Str_SmoothMuscle")) {
  results_list <- list()

  for (tp in timepoints) {
    subset_obj <- subset(input_obj,
                          IntType %in% cell_types &
                          timepoint %in% c(reference_timepoint, tp))
    Idents(subset_obj) <- "timepoint"

    # AR stratification
    subset_obj$Ar_positive <- GetAssayData(subset_obj, assay = "RNA", slot = "data")["Ar", ] >= 1
    subset_obj@meta.data   <- subset_obj@meta.data %>%
      mutate(Ar_positive = case_when(Ar_positive ~ "AR_positive", TRUE ~ "AR_negative"))

    for (ct in unique(subset_obj$FullType)) {
      sub_ct     <- subset(subset_obj, FullType == ct)
      sub_AR_pos <- subset(sub_ct, Ar_positive == "AR_positive")
      sub_AR_neg <- subset(sub_ct, Ar_positive == "AR_negative")

      results_list[[paste(tp, ct, "all", sep = "_")]] <- data.frame(
        timepoint = tp, FullType = ct, AR_status = "all",
        p_val_adj = tryCatch(get_ptn_padj(sub_ct,     tp, reference_timepoint), error = function(e) NA)
      )
      results_list[[paste(tp, ct, "ARpos", sep = "_")]] <- data.frame(
        timepoint = tp, FullType = ct, AR_status = "AR_positive",
        p_val_adj = tryCatch(get_ptn_padj(sub_AR_pos, tp, reference_timepoint), error = function(e) NA)
      )
      results_list[[paste(tp, ct, "ARneg", sep = "_")]] <- data.frame(
        timepoint = tp, FullType = ct, AR_status = "AR_negative",
        p_val_adj = tryCatch(get_ptn_padj(sub_AR_neg, tp, reference_timepoint), error = function(e) NA)
      )
    }
  }

  do.call(rbind, results_list)
}

Wouter_SmoothMuscle_MES_Glial <- readRDS("Wouter_SmoothMuscle_MES_Glial.RDS")
result <- get_ptn_pvalues_all_tp(Wouter_SmoothMuscle_MES_Glial)
write.csv(result, "output/Ptn_Seurat_adjusted_p_value_for_each_timepoint.csv", row.names = FALSE)


# =============================================================================
# 4. Human – ENZA Responder vs No-responder GSEA (Ext data Fig 6)
# =============================================================================
library(readxl); library(data.table); library(cowplot); library(grid)

Cell_death_hs <- read.gmt("data/hallmark_add_cell_death.Hs.v2503.gmt")

gene_list_enza <- read.delim(
  "data/GSEA/20Joshi_PNAS_mCRPC_hallmark_cell_death_responder_vs_noresponder.Gsea.1741380719376/ranked_gene_list_Responder_versus_No_responder_1741380719376.tsv"
)
gl_enza        <- sort(setNames(gene_list_enza$SCORE, gene_list_enza$NAME), decreasing = TRUE)

Responder_vs_Noresponder_gsea <- GSEA(gl_enza, TERM2GENE = Cell_death_hs,
                                        minGSSize = 5, pvalueCutoff = 1)
gsea_df <- Responder_vs_Noresponder_gsea@result
gsea_df$rank <- seq_len(nrow(gsea_df))

# Print key pathways
print(gsea_df %>%
  filter(Description %in% c("FERROPTOSIS", "REACTOME_PYROPTOSIS",
                             "GOBP_NECROPTOTIC_SIGNALING_PATHWAY", "HALLMARK_APOPTOSIS")))

# GSEA plots for selected gene sets (indices from gsea_df)
i_values <- c(25, 28, 36, 38)
pdf("output/cell_death_signature_GSEA_ENZA_treatment.pdf", width = 8, height = 6)
for (i in i_values) {
  anno <- Responder_vs_Noresponder_gsea[i, c("NES", "pvalue", "qvalue")]
  lab  <- paste0(names(anno), "=", round(anno, 3), collapse = "\n")

  p1 <- gseaplot2(Responder_vs_Noresponder_gsea, geneSetID = i, subplots = 1,
                   color = "black", title = gsea_df$Description[i]) +
    annotate("text", 40000, gsea_df$enrichmentScore[i] * 0.75,
             label = lab, hjust = 0, vjust = 0) +
    ylab("Enrichment Score") + xlab(NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  p2 <- gseaplot2(Responder_vs_Noresponder_gsea, geneSetID = i, subplots = 2,
                   color = "black", title = "") +
    ylab(NULL) + xlab(NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    annotation_custom(textGrob("Responder",    gp = gpar(fontsize = 13, col = "#ae272c")),
                       xmin = 3000,  xmax = 3000,  ymin = -0.1, ymax = -0.1) +
    annotation_custom(textGrob("No-responder", gp = gpar(fontsize = 13, col = "#2062a6")),
                       xmin = 50000, xmax = 50000, ymin = -0.1, ymax = -0.1) +
    coord_cartesian(clip = "off")

  print(p1 / p2 + plot_layout(heights = c(3, 1)))
}
dev.off()

# Ferroptosis only
pdf("output/ENZA_Responder_vs_NoResponder_FERROPTOSIS.pdf", width = 4, height = 4)
gseaplot2(Responder_vs_Noresponder_gsea, geneSetID = 25, color = "black",
           title = "FERROPTOSIS (p=0.014, q=0.016, NES=1.45)\nResponder vs NoResponder")
dev.off()


# =============================================================================
# 5. Human tumor cells – CSPC vs CRPC ferroptosis GSEA (Ext Data Fig 8)
# =============================================================================
library(clustree)

TumorOnly_obj <- readRDS("data/GSE264573_msk.integrated.remove.cellcycle.tumor.cells.rds")
Cell_death_hs2 <- read.gmt("data/hallmark_add_cell_death.Hs.v2503.gmt")

Tumor_CSPC_CRPC_diff <- FindMarkers(TumorOnly_obj,
                                      ident.1  = "CSPC", ident.2 = "CRPC",
                                      group.by = "subtype", min.pct = 0.2,
                                      logfc.threshold = 0.25)
gl_tumor <- rank_gene_list(Tumor_CSPC_CRPC_diff)

Tumor_gsea <- GSEA(gl_tumor, TERM2GENE = Cell_death_hs2, minGSSize = 2, pvalueCutoff = 1)
gsea_tumor_df        <- Tumor_gsea@result
gsea_tumor_df$rank   <- seq_len(nrow(gsea_tumor_df))

# Ferroptosis GSEA plot
pdf("output/CSPC_vs_CRPC_FERROPTOSIS.pdf", width = 4, height = 4)
gseaplot2(Tumor_gsea, color = "black", geneSetID = 15,
           title = "FERROPTOSIS (p=0.0028, q=0.0052, NES=1.98)\nTumor cells: CSPC vs CRPC")
dev.off()

message("GSEA analysis complete. Outputs saved to output/")
