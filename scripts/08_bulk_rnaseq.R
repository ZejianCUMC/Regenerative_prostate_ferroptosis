# =============================================================================
# 08_bulk_rnaseq.R
# Bulk RNA-seq analyses:
#   A. Mouse luminal cells – castration time course DESeq2 (Fig 5i volcano)
#   B. Human prostate – Pre vs Post ADT treatment (GSE48403)
#      MUFA/ferroptosis gene expression violin plots + DESeq2
#   C. Export expression matrices for GSEA (GCT/CLS format)
#
# Input:  Raw count tables, TPM files, metadata
# Output: DEG lists, violin plots, PDF figures, GCT/CLS files for GSEA
# =============================================================================

source("R/00_setup.R")
library(DESeq2)
library(ggrepel)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)


# =============================================================================
# A. Mouse luminal – castration time course (linear model over time)
#    Fig 5i: Volcano plot of genes downregulated over castration
# =============================================================================

# --- A1. Load and merge raw count tables ------------------------------------
# Load RSEM .genes.results files (one per sample in working directory)
load_rsem_counts <- function(pattern = "\\.genes\\.results$") {
  files <- list.files(pattern = pattern)
  raw_counts <- NULL
  tpm_counts <- NULL

  for (f in files) {
    data         <- read.delim(f, sep = "\t", header = TRUE)
    name         <- gsub(".genes.results", "", basename(f))
    count_subset <- data[, c(1, 5)]; colnames(count_subset)[2] <- name
    tpm_subset   <- data[, c(1, 6)]; colnames(tpm_subset)[2]   <- name

    if (is.null(raw_counts)) {
      raw_counts <- count_subset
      tpm_counts <- tpm_subset
    } else {
      raw_counts <- cbind(raw_counts, count_subset[2])
      tpm_counts <- cbind(tpm_counts, tpm_subset[2])
    }
  }

  # Map Ensembl IDs to gene symbols
  genes_v2           <- bitr(raw_counts$gene_id, fromType = "ENSEMBL",
                              toType = "SYMBOL", OrgDb = "org.Mm.eg.db")
  raw_counts$gene    <- genes_v2$SYMBOL[match(raw_counts$gene_id, genes_v2$ENSEMBL)]
  tpm_counts$gene    <- genes_v2$SYMBOL[match(tpm_counts$gene_id, genes_v2$ENSEMBL)]

  list(raw = raw_counts, tpm = tpm_counts)
}

# --- A2. Alternatively, load pre-merged count tables from xlsx ----------------
merged_raw_counts_batch1 <- readxl::read_xlsx("data/bulk/batch1_raw_count.xlsx")
merged_raw_counts_batch2 <- read.delim("data/bulk/batch2_raw_counts.txt", check.names = FALSE)

# Select luminal cell columns across three timepoints (0h, 48h, 72h post-castration)
merged_raw_counts_lum_tp <- as.data.frame(merged_raw_counts_batch1[, c(
  "lum1-c0h",  "lum2-c0h",  "lum3-c0h",
  "lum4-c48h", "lum5-c48h", "lum6-c48h",
  "lum7-c72h", "lum8-c72h", "lum9-c72h"
)])
rownames(merged_raw_counts_lum_tp) <- merged_raw_counts_batch1$gene

# Build colData with continuous time variable
sample_names <- colnames(merged_raw_counts_lum_tp)
time_numeric <- as.numeric(sub("h", "", sub(".*-c", "", sample_names)))
metaData <- data.frame(
  row.names  = sample_names,
  time       = time_numeric,
  replicate  = rep(paste0("rep", 1:3), times = 3),
  ScaleTime  = scale(time_numeric)
)

# --- A3. DESeq2 linear model over scaled time --------------------------------
dds_lum <- DESeqDataSetFromMatrix(
  countData = merged_raw_counts_lum_tp,
  colData   = metaData,
  design    = ~ ScaleTime
)
dds_lum <- DESeq(dds_lum)
res_lum  <- results(dds_lum, name = "ScaleTime")

# Genes downregulated over castration (Fig 5i volcano)
downreg <- as.data.frame(subset(res_lum, padj < 0.05 & log2FoldChange < 0))
downreg <- downreg[order(downreg$log2FoldChange), ]
downreg$Index <- seq_len(nrow(downreg))
write.csv(downreg, "output/Mouse_bulk_luminal_downregulated_over_castration.csv")
message("Nkx3-1 position: ", downreg["Nkx3-1", "Index"])


# =============================================================================
# B. Human prostate – Pre vs Post ADT treatment (GSE48403)
# =============================================================================

# --- B1. Load TPM and raw count data ----------------------------------------
GSE48403_tpm   <- read.delim("data/GSE48403_ADT_RNA/GSE48403_norm_counts_TPM_GRCh38.p13_NCBI.tsv")
ncbi_anno      <- read.delim("data/GSE48403_ADT_RNA/Human.GRCh38.p13.annot.tsv")
GSE48403_count <- read.delim("data/GSE48403_ADT_RNA/GSE48403_refcoding_read_counts.txt",
                              check.names = FALSE)

# Sample metadata
meta_info <- data.frame(
  GSM        = c("GSM1177208","GSM1177209","GSM1177210","GSM1177211","GSM1177212",
                 "GSM1177213","GSM1177214","GSM1177215","GSM1177216","GSM1177217",
                 "GSM1177218","GSM1177219","GSM1177220","GSM1177221"),
  SampleInfo = c("Prostate-Pre-R1","Prostate-Pre-R2","Prostate-Pre-R3",
                 "Prostate-Pre-R4","Prostate-Pre-R5","Prostate-Pre-R6",
                 "Prostate-Pre-R7","Prostate-Post-R1","Prostate-Post-R2",
                 "Prostate-Post-R3","Prostate-Post-R4","Prostate-Post-R5",
                 "Prostate-Post-R6","Prostate-Post-R7"),
  stringsAsFactors = FALSE
)
meta_info <- meta_info %>%
  mutate(
    Treatment = ifelse(grepl("Pre", SampleInfo), "Pre", "Post"),
    PatientID = stringr::str_split(SampleInfo, "-", simplify = TRUE)[, 3]
  )
rownames(meta_info) <- meta_info$SampleInfo

# --- B2. Map gene IDs to symbols and collapse duplicates --------------------
GSE48403_tpm$Gene <- ncbi_anno$Symbol[match(GSE48403_tpm$GeneID, ncbi_anno$GeneID)]
GSE48403_tpm      <- GSE48403_tpm %>%
  filter(!is.na(Gene) & Gene != "") %>%
  select(-GeneID) %>%
  group_by(Gene) %>%
  summarise(across(everything(), mean))

# --- B3. Visualize MUFA / ferroptosis-related genes -------------------------
vis_genes <- intersect(c(MUFA_gene_human, "PTN"), GSE48403_tpm$Gene)
vis_df <- GSE48403_tpm %>%
  filter(Gene %in% vis_genes) %>%
  pivot_longer(!Gene, names_to = "GSM", values_to = "TPM") %>%
  left_join(meta_info, by = c("GSM" = "GSM")) %>%
  mutate(log2TPM  = log2(TPM + 1),
         Treatment = factor(Treatment, levels = c("Pre", "Post")))

p_mufa_vln <- ggviolin(vis_df, x = "Treatment", y = "log2TPM",
                         fill = "Treatment", color = "Treatment",
                         trim = TRUE, size = 1) +
  geom_boxplot(fill = "white", width = 0.1) +
  facet_wrap(~Gene, scales = "free", nrow = 3) +
  stat_compare_means(comparisons = list(c("Pre", "Post")),
                      method = "t.test", paired = TRUE) +
  labs(title = "MUFA-related gene expression: Pre vs Post ADT (GSE48403)",
       y = "log2(TPM + 1)") +
  scale_fill_manual(values  = c("#5d5af7", "#ff0202")) +
  scale_color_manual(values = c("#5d5af7", "#ff0202")) +
  theme(legend.position = "none",
        axis.title      = element_text(size = 10),
        axis.text       = element_text(size = 10, colour = "black"),
        plot.title      = element_text(size = 10, hjust = 0.5))

ggsave("output/GSE48403_pre_post_ADT_MUFA_log2TPM_paired_ttest.pdf",
       p_mufa_vln, width = 12, height = 12)

# --- B4. DESeq2: Post vs Pre ------------------------------------------------
gene_id_map <- bitr(GSE48403_count$ID, fromType = "ENSEMBL",
                     toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
GSE48403_count$gene <- gene_id_map$SYMBOL[match(GSE48403_count$ID, gene_id_map$ENSEMBL)]

raw_counts_filter <- GSE48403_count %>%
  filter(!is.na(gene) & gene != "") %>%
  select(-ID) %>%
  data.table::as.data.table() %>%
  `[`(, lapply(.SD, median), by = gene) %>%
  as.data.frame()
rownames(raw_counts_filter) <- raw_counts_filter$gene
raw_counts_filter$gene      <- NULL

meta_info_update <- meta_info[colnames(raw_counts_filter), ]

dds_adt <- DESeqDataSetFromMatrix(
  countData = round(raw_counts_filter),
  colData   = meta_info_update,
  design    = ~ Treatment
)
dds_adt <- estimateSizeFactors(dds_adt)
dds_adt <- DESeq(dds_adt)
res_adt  <- results(dds_adt, contrast = c("Treatment", "Post", "Pre"))

res_adt_df <- as.data.frame(res_adt)
mufa_hit   <- intersect(MUFA_gene_human, rownames(res_adt_df))
print(res_adt_df[mufa_hit, ])
write.csv(res_adt_df, "output/GSE48403_DESeq2_Post_vs_Pre_ADT.csv")


# =============================================================================
# C. Export expression matrices for GSEA (GCT + CLS format)
# =============================================================================

#' Write a GCT-format file for GSEA
write_gct <- function(expr_df, filename) {
  gct           <- cbind(expr_df$gene, expr_df)
  colnames(gct)[1:2] <- c("gene_name", "Description")
  n_genes       <- nrow(gct)
  n_samples     <- ncol(gct) - 2
  header1       <- c("#1.2", rep("", n_samples + 1))
  header2       <- c(n_genes, n_samples, rep("", n_samples))
  names(header1) <- names(header2) <- colnames(gct)
  out <- rbind(header1, header2, setNames(data.frame(t(colnames(gct))), colnames(gct)), gct)
  write.table(out, filename, row.names = FALSE, col.names = FALSE,
              quote = FALSE, sep = "\t")
}

#' Write a CLS phenotype file for GSEA
write_cls <- function(sample_labels, classes, filename) {
  cls <- data.frame(
    row1 = c(length(sample_labels), length(classes), 1, rep("", length(sample_labels) - 3)),
    row2 = c("#", classes, rep("", length(sample_labels) - length(classes) - 1)),
    row3 = sample_labels
  )
  cls <- data.frame(t(cls))
  write.table(cls, filename, row.names = FALSE, col.names = FALSE,
              quote = FALSE, sep = "\t")
}

# Human cohort GCT
write_gct(as.data.frame(GSE48403_tpm), "output/GSE48403_tpm.gct")

# Mouse batch1 luminal GCT
lum_tpm <- merged_raw_counts_batch1[, c("gene",
                                          "lum1-c0h","lum2-c0h","lum3-c0h",
                                          "lum4-c48h","lum5-c48h","lum6-c48h",
                                          "lum7-c72h","lum8-c72h","lum9-c72h")]
write_gct(lum_tpm, "output/Mouse_bulk_batch1_lum.gct")

chars_lum <- c("0h", "48h", "72h")
write_cls(rep(chars_lum, each = 3), chars_lum, "output/Mouse_bulk_batch1_lum.cls")

# Mouse batch2 GCT
batch2_tpm <- merged_raw_counts_batch2
write_gct(batch2_tpm, "output/Mouse_bulk_batch2.gct")
chars_b2 <- c("AR", "CAS", "DKO", "LKO", "NKX31KO", "SKO")
write_cls(rep(chars_b2, each = 4), chars_b2, "output/Mouse_bulk_batch2.cls")

message("Bulk RNA-seq analysis complete. Outputs saved to output/")
