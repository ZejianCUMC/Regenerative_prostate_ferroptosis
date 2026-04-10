# =============================================================================
# 05_luminal_signature_visualization.R
# Luminal signature scatter plots for all figure panels (Figs 1e, 2, 3e, 4–5)
#
# Each panel plots cells in LumA-score × LumP-score space, coloured by
# classifier probability (LumA = blue, LumP = orange-red gradient).
#
# Input:  Signature data frames from 04_luminal_signature_scoring.R
# Output: PDFs and CSV exports for each figure panel
# =============================================================================

source("R/00_setup.R")
library(ggnewscale)


# =============================================================================
# Shared plot theme and gradient scales
# =============================================================================

theme_signature <- function() {
  theme_classic() +
    theme(
      legend.position  = "right",
      axis.title       = element_text(size = 20, hjust = 0.5),
      axis.text.x      = element_text(size = 10, hjust = 0.5),
      plot.title       = element_text(size = 20, hjust = 0.5),
      strip.text       = element_text(size = 14, face = "bold"),
      legend.text      = element_text(size = 15),
      legend.title     = element_text(size = 20),
      strip.background = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line        = element_line(size = 0.5, colour = "black")
    )
}

SCALE_LUMA <- scale_colour_gradientn(
  colours = c("#caf0f8", "#7989ff", "#0058b8"), limits = c(0, 1)
)
SCALE_LUMP <- scale_colour_gradientn(
  colours = c("#d8e2dc", "#e3b938", "#ca3e14"), limits = c(0, 1)
)

COORD_SIG <- coord_cartesian(xlim = c(-4.1, 4.1), ylim = c(-3, 5))

#' Core scatter plot with dual classifier-probability colour scales
#'
#' @param df         Data frame with Luminal_A, Luminal_P, LumType, PredictionProb columns
#' @param facet_var  Facetting variable name (string)
#' @param dot_size   Point size (default 0.5)
plot_signature <- function(df, facet_var, dot_size = 0.5) {
  df_lumA <- df %>% filter(LumType == "LumA")
  df_lumP <- df %>% filter(LumType == "LumP")

  ggplot(df, aes(Luminal_A, Luminal_P)) +
    geom_point(data = df_lumA, aes(colour = PredictionProb), size = dot_size) +
    SCALE_LUMA +
    labs(color = "Classifier Pr(LumA)") +
    new_scale_colour() +
    geom_point(data = df_lumP, aes(colour = PredictionProb), size = dot_size) +
    SCALE_LUMP +
    facet_wrap(as.formula(paste("~", facet_var))) +
    labs(x = "LumA signature", y = "LumP signature",
         color = "Classifier Pr(LumP)") +
    theme_signature() +
    COORD_SIG
}

save_plot <- function(p, filename, width = 8, height = 8) {
  ggsave(filename = file.path("output", filename),
         plot = p, width = width, height = height, dpi = 300)
}


# =============================================================================
# Load all signature data frames
# =============================================================================
Wouter_sig <- readRDS("Wouter_science_luminal_signature.RDS")
ML001_sig  <- readRDS("ML001_Luminal_signature_df_lab.RDS")
ML002_sig  <- readRDS("ML002_Luminal_signature_df_lab.RDS")
MW3_sig    <- readRDS("MW3_Luminal_signature_df.RDS")
MW4_sig    <- readRDS("MW4_Luminal_signature_df.RDS")
MW6_sig    <- readRDS("MW6_Luminal_signature_df.RDS")
MW7_sig    <- readRDS("MW7_Luminal_signature_df.RDS")
MW8_sig    <- readRDS("MW8_Luminal_signature_df.RDS")

# Standardise LumType in Wouter data
Wouter_sig <- Wouter_sig %>%
  mutate(LumType = case_when(FullType == "Epi_Luminal_1" ~ "LumA",
                              FullType == "Epi_Luminal_2Psca" ~ "LumP"))


# =============================================================================
# Fig 1e: Castration time course (Intact → 3d → 7d → 4w)
# =============================================================================
plot_df_fig1 <- bind_rows(
  MW3_sig %>%
    filter(CellType == "FVB WT cas 3d") %>%
    select(Luminal_A, Luminal_P, LumType, PredictionProb, LumA_pro, LumP_pro) %>%
    mutate(timepoint = "T01_Cast_Day3"),
  Wouter_sig %>%
    filter(timepoint %in% c("T00_Intact", "T02_Cast_Day7", "T04_Cast_Day28")) %>%
    mutate(LumType = case_when(FullType == "Epi_Luminal_1" ~ "LumA", TRUE ~ "LumP")) %>%
    select(Luminal_A, Luminal_P, LumType, PredictionProb, LumA_pro, LumP_pro, timepoint)
)
plot_df_fig1$timepoint <- factor(plot_df_fig1$timepoint,
                                  levels = c("T00_Intact", "T01_Cast_Day3",
                                             "T02_Cast_Day7", "T04_Cast_Day28"))

p_fig1 <- plot_signature(plot_df_fig1, "timepoint")
save_plot(p_fig1, "plot1_Cast_day0to28_smalldot_add_probability.pdf")
write.csv(plot_df_fig1, "output/plot1_Cast_day0to28_data.csv", row.names = FALSE)


# =============================================================================
# Fig 2: AR-intact, AR-Cas 4w, Stromal AR KO, Luminal AR KO
# =============================================================================
plot_df_fig2 <- bind_rows(
  MW3_sig %>%
    filter(CellType %in% c("AR-intact", "AR Castration 4w")) %>%
    select(Luminal_A, Luminal_P, LumType, CellType, PredictionProb, LumA_pro, LumP_pro),
  ML001_sig %>%
    filter(timepoint == "4w") %>%
    mutate(CellType = case_when(CellType == "Luminal AR KO" ~ "Luminal AR KO 4w",
                                 TRUE ~ "Stromal AR KO 4w")) %>%
    select(Luminal_A, Luminal_P, LumType, CellType, PredictionProb, LumA_pro, LumP_pro)
)
plot_df_fig2$CellType <- factor(plot_df_fig2$CellType,
                                 levels = c("AR-intact", "AR Castration 4w",
                                            "Stromal AR KO 4w", "Luminal AR KO 4w"))

p_fig2 <- plot_signature(plot_df_fig2, "CellType")
save_plot(p_fig2, "plot2_AR_ARcas_StromKO_LumKO_add_probability.pdf")
write.csv(plot_df_fig2, "output/plot2_AR_ARcas_StromKO_LumKO_data.csv", row.names = FALSE)


# =============================================================================
# Fig 3e: WT AR, WT AR Cas d28, Luminal KO, Stromal KO, Double KO
# =============================================================================
plot_df_fig3 <- bind_rows(
  MW3_sig %>%
    filter(CellType %in% c("AR-intact", "AR Castration 4w", "Stromal and Luminal AR KO")) %>%
    select(Luminal_A, Luminal_P, LumType, CellType, PredictionProb, LumA_pro, LumP_pro),
  ML001_sig %>%
    filter(timepoint == "4w") %>%
    mutate(CellType = case_when(CellType == "Luminal AR KO" ~ "Luminal AR KO 4w",
                                 TRUE ~ "Stromal AR KO 4w")) %>%
    select(Luminal_A, Luminal_P, LumType, CellType, PredictionProb, LumA_pro, LumP_pro)
)
p_fig3 <- plot_signature(plot_df_fig3, "CellType")
save_plot(p_fig3, "plot3_AR_ARcas_StromKO_LumKO_doubleKO_add_probability.pdf")
write.csv(plot_df_fig3, "output/plot3_AR_ARcas_StromKO_LumKO_doubleKO_data.csv",
          row.names = FALSE)


# =============================================================================
# Ext Data Fig 4: Atg7-WT, Atg7-Ko intact, Atg7-WT Cas 4w, Atg7-KO 4w
# =============================================================================
plot_df_fig4 <- bind_rows(
  MW4_sig %>%
    filter(HTO_maxID %in% c("B301", "B302", "B303")) %>%
    select(LumType, CellType, Luminal_A, Luminal_P, LumA_pro, LumP_pro, PredictionProb),
  ML002_sig %>%
    filter(HTO_maxID == "B304") %>%
    mutate(CellType = "Atg7-KO 4w") %>%
    select(LumType, CellType, Luminal_A, Luminal_P, LumA_pro, LumP_pro, PredictionProb)
)
plot_df_fig4$CellType <- factor(plot_df_fig4$CellType,
                                 levels = c("Atg7-WT intact", "Atg7-Ko intact",
                                            "Atg7-WT intact Castration 4w", "Atg7-KO 4w"))
p_fig4 <- plot_signature(plot_df_fig4, "CellType")
save_plot(p_fig4, "plot4_Luminal_Atg7_WT_KO_add_probability.pdf")
write.csv(plot_df_fig4, "output/plot4_Luminal_Atg7_WT_KO_data.csv", row.names = FALSE)


# =============================================================================
# Fig 4 extension: AR-Cast 4w, Lum_ARKO_4w, Nkx3-1 +/+ 4w, Nkx3-1 -/- 4w
# =============================================================================
plot_df_fig5 <- bind_rows(
  MW3_sig %>%
    filter(CellType == "AR Castration 4w") %>%
    select(LumType, CellType, Luminal_A, Luminal_P, LumA_pro, LumP_pro, PredictionProb),
  ML001_sig %>%
    filter(timepoint == "4w", CellType == "Luminal AR KO") %>%
    mutate(CellType = "Luminal AR KO 4w") %>%
    select(LumType, CellType, Luminal_A, Luminal_P, LumA_pro, LumP_pro, PredictionProb),
  MW4_sig %>%
    filter(CellType %in% c("NKX3.1+/+", "NKX3.1-/-")) %>%
    mutate(CellType = paste0(CellType, " 4w")) %>%
    select(LumType, CellType, Luminal_A, Luminal_P, LumA_pro, LumP_pro, PredictionProb)
)
plot_df_fig5$CellType <- factor(plot_df_fig5$CellType,
                                 levels = c("AR Castration 4w", "Luminal AR KO 4w",
                                            "NKX3.1+/+ 4w", "NKX3.1-/- 4w"))
p_fig5 <- plot_signature(plot_df_fig5, "CellType")
save_plot(p_fig5, "plot5_AR_cas_4w_luminal_AR_KO_4w_NKX3-1_4w_add_probability.pdf")
write.csv(plot_df_fig5, "output/plot5_AR_cas_4w_luminal_AR_KO_4w_NKX3-1_4w_data.csv",
          row.names = FALSE)


# =============================================================================
# Fig 4 panel 6: AR-Intact, Lum_ARKO_4w, Nkx3-1 +/+ 4w, Nkx3-1 -/- 4w
# =============================================================================
plot_df_fig6 <- bind_rows(
  plot_df_fig3 %>% filter(CellType == "AR-intact"),
  plot_df_fig5 %>% filter(!CellType %in% c("AR Castration 4w"))
)
plot_df_fig6$CellType <- factor(plot_df_fig6$CellType,
                                 levels = c("AR-intact", "Luminal AR KO 4w",
                                            "NKX3.1+/+ 4w", "NKX3.1-/- 4w"))
p_fig6 <- plot_signature(plot_df_fig6, "CellType")
save_plot(p_fig6, "plot6_AR_Intact_luminal_AR_KO_4w_NKX3-1_4w_add_probability.pdf")
write.csv(plot_df_fig6, "output/plot6_AR_Intact_luminal_AR_KO_4w_NKX3-1_4w_data.csv",
          row.names = FALSE)


# =============================================================================
# Fig (Ferroptosis intervention): Intact, Liproxstatin, Cast4w, Cast4w+Liprox
# =============================================================================
p_fig7 <- plot_signature(MW6_sig, "Batch")
save_plot(p_fig7, "plot7_Intact_Liproxstatin_Cast4w_add_probability.pdf")
write.csv(MW6_sig %>% select(cluster, Batch, Luminal_A, Luminal_P, LumA_pro, LumP_pro, PredictionProb),
          "output/plot7_Intact_Liproxstatin_Cast4w_data.csv", row.names = FALSE)


# =============================================================================
# Intact, Intact_VE_High, Cas4w, Cas4w_VE_High  (MW6 DMSO + MW7 VE high)
# =============================================================================
plot_df_fig8 <- bind_rows(
  MW6_sig %>% filter(Batch %in% c("WT_intact_DMSO", "WT_Cas4w_DMSO")),
  MW7_sig %>% filter(Batch %in% c("WT_intact_VE_high", "WT_Cas4w_VE_high"))
)
plot_df_fig8$Batch <- factor(plot_df_fig8$Batch,
                              levels = c("WT_intact_DMSO", "WT_intact_VE_high",
                                         "WT_Cas4w_DMSO",  "WT_Cas4w_VE_high"))
p_fig8 <- plot_signature(plot_df_fig8, "Batch")
save_plot(p_fig8, "plot8_Intact_VEhigh_Cas4w_add_probability.pdf")
write.csv(plot_df_fig8 %>% select(cluster, Batch, Luminal_A, Luminal_P, LumA_pro, LumP_pro, PredictionProb),
          "output/plot8_Intact_VEhigh_Cas4w_data.csv", row.names = FALSE)


# =============================================================================
# Fig 5c: AR-Intact, LumARKO_4w, Raptor KO 4w, Raptor AR DKO 4w
# =============================================================================
plot_df_fig10 <- bind_rows(
  MW8_sig %>%
    filter(Batch %in% c("RAPTOR_KO_4w", "RAPTOR_AR_DKO_4w")) %>%
    mutate(LumType = ClusterUpdate) %>%
    select(LumType, Batch, Luminal_A, Luminal_P, LumA_pro, LumP_pro, PredictionProb),
  plot_df_fig6 %>%
    filter(CellType %in% c("AR-intact", "Luminal AR KO 4w")) %>%
    mutate(Batch = CellType) %>%
    select(LumType, Batch, Luminal_A, Luminal_P, LumA_pro, LumP_pro, PredictionProb)
)
plot_df_fig10$Batch <- factor(plot_df_fig10$Batch,
                               levels = c("AR-intact", "Luminal AR KO 4w",
                                          "RAPTOR_KO_4w", "RAPTOR_AR_DKO_4w"))
p_fig10 <- plot_signature(plot_df_fig10, "Batch")
save_plot(p_fig10, "plot10_ARWT_LumARKO4w_RaptorKO4w_RaptorARDKO4w_add_probability.pdf")
write.csv(plot_df_fig10, "output/plot10_ARWT_LumARKO4w_RaptorKO4w_RaptorARDKO4w_data.csv",
          row.names = FALSE)


# =============================================================================
# Fig 5i: PTN_WT_4w, PTN_KO_4w, NKX3-1-/-, PTN_NKX_DKO_4w
# =============================================================================
plot_df_fig11 <- bind_rows(
  MW8_sig %>%
    filter(Batch %in% c("PTN_KO_4w", "PTN_NKX_DKO_4w")) %>%
    mutate(LumType = ClusterUpdate) %>%
    select(LumType, Batch, Luminal_A, Luminal_P, LumA_pro, LumP_pro, PredictionProb),
  plot_df_fig6 %>%
    filter(CellType %in% c("NKX3.1+/+ 4w", "NKX3.1-/- 4w")) %>%
    mutate(Batch = case_when(CellType == "NKX3.1+/+ 4w" ~ "PTN_WT_4w", TRUE ~ CellType)) %>%
    select(LumType, Batch, Luminal_A, Luminal_P, LumA_pro, LumP_pro, PredictionProb)
)
plot_df_fig11$Batch <- factor(plot_df_fig11$Batch,
                               levels = c("PTN_WT_4w", "PTN_KO_4w",
                                          "NKX3.1-/- 4w", "PTN_NKX_DKO_4w"))
p_fig11 <- plot_signature(plot_df_fig11, "Batch")
save_plot(p_fig11, "plot11_PTNWT4w_PTNKO4w_NKXKO4w_PTNNKXDKO4w_add_probability.pdf")
write.csv(plot_df_fig11, "output/plot11_PTNWT4w_PTNKO4w_NKXKO4w_PTNNKXDKO4w_data.csv",
          row.names = FALSE)


message("All figure panels saved to output/")
