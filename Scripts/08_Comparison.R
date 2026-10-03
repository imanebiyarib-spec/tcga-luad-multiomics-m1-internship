#' =============================================================================
#' Project:  Genomics Project - Phase 1 (LUAD)
#' Task:     Comprehensive Architecture Comparison & Visual Syntheses
#' Author:   Imane BIYAR
#' Date:     August 2026
#'
#' Objective:
#' Synthesize cross-validation and held-out test metrics across all THREE 
#' architectures (Random Forest, Elastic Net, and XGBoost). Generate publication 
#' tables (including MCC) and visualize the multi-omics driver rescue landscape.
#'
#' Inputs:
#'   - Results/Integration/master_table_pathways.rds
#'   - Results/ML/*.csv (Model summaries, coefficients)
#'   - Results/Integration/*.csv (Ranked candidate panels)
#'
#' Outputs:
#'   - Results/Summary/Table1_TriModel_Performance_Benchmark.csv
#'   - Results/Summary/Table2_Biological_Accounting_Funnel.csv
#'   - Results/Summary/Vis1_ElasticNet_Mechanics.png
#'   - Results/Summary/Vis2_MultiOmics_Rescue_Landscape_RF.png
#'   - Results/Summary/Vis3_MultiOmics_Rescue_Landscape_XGBoost.png
#' =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

dir.create("Results/Summary", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 8.1 - Load ML Synthesis Data & Integrated Backbone
# ==============================================================================
master_table    <- readRDS("Results/Integration/master_table_pathways.rds")
cv_summary      <- read.csv("Results/ML/cv_summary_all_models.csv")
test_summary    <- read.csv("Results/ML/held_out_test_summary_all_models.csv")
glmnet_coefs    <- read.csv("Results/ML/glmnet_coefficients.csv")

top_10_df       <- read.csv("Results/Integration/top_10_drivers.csv")
top_50_drivers  <- readRDS("Results/Integration/top_50_drivers.rds")
top_10_novel    <- read.csv("Results/Integration/top_10_novel_for_literature_validation.csv")
positives       <- readRDS("Results/ML/positive.rds")

top_10_df_xgb    <- read.csv("Results/Integration/top_10_drivers_XGBoost.csv")
top_10_novel_xgb <- read.csv("Results/Integration/top_10_novel_for_literature_validation_XGBoost.csv")
ranked_xgb_data  <- read.csv("Results/Integration/tier2_ranked_by_XGBoost.csv")

# ==============================================================================
# 8.2 - Table 1: Tri-Model Performance Benchmark (CV vs. Held-Out Test)
# ==============================================================================

# Format long CV summary table into publication-ready Mean ± SD strings
cv_formatted <- cv_summary %>%
  mutate(formatted_val = sprintf("%.3f ± %.3f", mean, sd)) %>%
  select(Model, Metric, formatted_val) %>%
  pivot_wider(names_from = Metric, values_from = formatted_val) %>%
  rename(
    ROC_CV         = ROC,
    Accuracy_CV    = Accuracy,
    BalancedAcc_CV = BalancedAccuracy,
    Sensitivity_CV = Sensitivity,
    Specificity_CV = Specificity,
    F1_CV          = F1,
    MCC_CV         = MCC
  )

# Format single held-out test set metrics
test_formatted <- test_summary %>%
  transmute(
    Model            = Model,
    ROC_Test         = sprintf("%.3f", ROC),
    Accuracy_Test    = sprintf("%.3f", Accuracy),
    BalancedAcc_Test = sprintf("%.3f", BalancedAccuracy),
    Sensitivity_Test = sprintf("%.3f", Sensitivity),
    Specificity_Test = sprintf("%.3f", Specificity),
    F1_Test          = sprintf("%.3f", F1),
    MCC_Test         = sprintf("%.3f", MCC)
  )

# Combine into unified master benchmark table
table1_benchmark <- left_join(cv_formatted, test_formatted, by = "Model")
write.csv(table1_benchmark, "Results/Summary/Table1_TriModel_Performance_Benchmark.csv", row.names = FALSE)

# ==============================================================================
# 8.3 - Table 2: Multi-Omics Accounting & Candidate Funnel
# ==============================================================================

n_total  <- nrow(master_table)
n_deg    <- sum(master_table$in_DEG, na.rm = TRUE)
n_cna    <- sum(master_table$in_CNA_gistic, na.rm = TRUE)
n_mut    <- sum(master_table$in_MUT, na.rm = TRUE)
n_tier1  <- sum(master_table$in_tier1, na.rm = TRUE)
n_tier2  <- sum(master_table$in_tier2, na.rm = TRUE)

table2_funnel <- data.frame(
  Pipeline_Stage = c(
    "1. Total Gene Universe (Multi-Omics Backbone Union)",
    "2. Transcriptomic Outliers (DESeq2 padj < 0.001 & |log2FC| > 2)",
    "3. Chromosomal Copy-Number Drivers (GISTIC2.0 Spatial Peaks)",
    "4. Somatic Evolutionary Selection (dN/dS q < 0.01 & w > 1)",
    "5. Static Intersection: Tier 1 (3-of-3 Strict Boolean Agreement)",
    "6. ML Candidate Pool: Tier 2 (>=2-of-3 Omics Agreement)"
  ),
  Gene_Count = c(n_total, n_deg, n_cna, n_mut, n_tier1, n_tier2)
)

write.csv(table2_funnel, "Results/Summary/Table2_Biological_Accounting_Funnel.csv", row.names = FALSE)

# ==============================================================================
# 8.4 - Visualization 1: Linear vs. Non-Linear Feature Mechanics
# ==============================================================================

en_df <- glmnet_coefs %>%
  filter(Variable != "(Intercept)") %>%
  mutate(
    Weight    = abs(Coefficient),
    Direction = ifelse(Coefficient > 0, "Positive / Oncogenic", "Negative / Suppressive")
  )

p_en <- ggplot(en_df, aes(x = reorder(Variable, Weight), y = Coefficient, fill = Direction)) +
  geom_col(width = 0.7, color = "black", alpha = 0.85) +
  coord_flip() +
  scale_fill_manual(values = c("Positive / Oncogenic" = "#D32F2F", "Negative / Suppressive" = "#1976D2")) +
  theme_bw(base_size = 12) +
  labs(title = "Elastic Net Mechanics: Standardized Log-Odds Coefficients",
       subtitle = "Evaluating linear directional contribution and L1/L2 penalty shrinkage",
       x = NULL, y = "Standardized Coefficient Magnitude") +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave("Results/Summary/Vis1_ElasticNet_Mechanics.png", plot = p_en, width = 8, height = 6, dpi = 300)

# ==============================================================================
# 8.5 - Visualization 2: The Multi-Omics Rescue Landscape (Random Forest)
# ==============================================================================

plot_df <- master_table %>%
  filter(in_tier2 == TRUE) %>%
  mutate(
    Label_Category = case_when(
      gene %in% top_10_novel$gene ~ "Top 10 Novel ML Candidate",
      gene %in% positives & gene %in% top_10_df$gene ~ "Top 10 Known Classical Driver",
      in_tier1 == TRUE ~ "Tier 1 Strict Overlap (SPRR1B)",
      TRUE ~ "Tier 2 Candidate Pool"
    ),
    Mut_Size = pmax(mut_count, 1)
  )

ranked_tier2_data <- read.csv("Results/Integration/tier2_ranked_by_RF.csv")
plot_df <- left_join(plot_df, ranked_tier2_data[, c("gene", "rf_probability")], by = "gene") %>%
  mutate(rf_probability = ifelse(is.na(rf_probability), 0, rf_probability))

p_landscape <- ggplot(plot_df, aes(x = log2FC, y = Pct_Amplified, size = Mut_Size, color = rf_probability)) +
  geom_point(alpha = 0.8) +
  scale_color_gradientn(colors = c("#E3F2FD", "#42A5F5", "#1565C0", "#B71C1C"), name = "RF Driver\nProbability") +
  scale_size_continuous(range = c(2, 9), name = "Somatic\nMutations") +
  geom_label_repel(data = filter(plot_df, Label_Category == "Top 10 Novel ML Candidate"),
                   aes(label = gene), size = 3.5, fontface = "bold", color = "#B71C1C",
                   box.padding = 0.5, point.padding = 0.3, max.overlaps = 50, show.legend = FALSE) +
  geom_label_repel(data = filter(plot_df, Label_Category == "Top 10 Known Classical Driver"),
                   aes(label = gene), size = 3.5, fontface = "bold", color = "#0D47A1",
                   box.padding = 0.5, point.padding = 0.3, max.overlaps = 50, show.legend = FALSE) +
  geom_label_repel(data = filter(plot_df, in_tier1 == TRUE),
                   aes(label = paste0(gene, " (Tier 1)")), size = 3.5, fontface = "italic",
                   color = "black", fill = "#FFF9C4", box.padding = 0.8, show.legend = FALSE) +
  theme_minimal(base_size = 12) +
  labs(title = "Multi-Omics Candidate Landscape of Lung Adenocarcinoma (TCGA-LUAD)",
       subtitle = "Random Forest rescue of classical oncogenes and novel actionable targets from Tier 2",
       x = "Transcriptional Overexpression (RNA-seq log2 Fold Change)",
       y = "Cohort Copy-Number Amplification Frequency (Pct_Amplified)") +
  theme(plot.title = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "gray80", fill = NA))

ggsave("Results/Summary/Vis2_MultiOmics_Rescue_Landscape_RF.png", plot = p_landscape, width = 12, height = 8, dpi = 300)

# ==============================================================================
# 8.6 - Visualization 3: The Multi-Omics Rescue Landscape (XGBoost)
# ==============================================================================

# 1. Data preparation and definition of dominant CNA shapes
plot_df_xgb <- master_table %>%
  filter(in_tier2 == TRUE | gene %in% positives) %>%
  left_join(ranked_xgb_data[, c("gene", "xgb_probability")], by = "gene") %>%
  mutate(
    xgb_probability = ifelse(is.na(xgb_probability), 0, xgb_probability),
    
    # Define dominant CNA state based on the 10% threshold methodology
    CNA_Status = case_when(
      Pct_Amplified >= Pct_AnyDel & Pct_Amplified >= 10 ~ "Amplified",
      Pct_AnyDel > Pct_Amplified & Pct_AnyDel >= 10 ~ "Deleted",
      TRUE ~ "Neutral"
    ),
    
    Label_Category = case_when(
      gene %in% positives & in_tier2 == TRUE ~ "Known_Tier2",              
      gene %in% positives & in_tier2 == FALSE ~ "Known_Missed",            
      !(gene %in% positives) & in_tier2 == TRUE & xgb_probability > 0.5 ~ "Novel_HighProb", 
      in_tier1 == TRUE ~ "Tier1",                                          
      TRUE ~ "Other"                                                       
    ),
    
    Mut_Size = pmax(mut_count, 1, na.rm = TRUE)
  )

plot_df_xgb$CNA_Status <- factor(plot_df_xgb$CNA_Status, levels = c("Amplified", "Neutral", "Deleted"))

# 2. Create a dummy dataframe with VALID coordinates to avoid the 'NA' axis bug
dummy_legend <- data.frame(
  Category = factor(c("Novel Candidate (Prob > 0.5)", "Known Driver (in Tier 2)", "Known Driver (Missed)"),
                    levels = c("Novel Candidate (Prob > 0.5)", "Known Driver (in Tier 2)", "Known Driver (Missed)")),
  x = 0, y = 1 
)

# 3. Plot Generation
p_landscape_xgb <- ggplot() +
  
  # Reference lines for biological interpretation
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray60", linewidth = 0.5, alpha = 0.7) +
  # --- THE HACK: Dummy layer for the Category Legend ---
  # We add alpha = 0 so the squares are 100% invisible on the plot itself
  geom_point(data = dummy_legend, aes(x = x, y = y, fill = Category), shape = 22, size = 4, color = "transparent", alpha = 0) +
  scale_fill_manual(name = "Gene Categories",
                    values = c("Novel Candidate (Prob > 0.5)" = "#D32F2F", 
                               "Known Driver (in Tier 2)" = "#2E7D32", 
                               "Known Driver (Missed)" = "#FF9800")) +
  
  # Layer A: Tier 2 Genes (Color = XGBoost Probability, Shape = CNA Status)
  geom_point(data = filter(plot_df_xgb, in_tier2 == TRUE), 
             aes(x = log2FC, y = pmax(mut_count, 1), color = xgb_probability, shape = CNA_Status), 
             size = 3, alpha = 0.8) +
  scale_color_gradientn(colors = c("#E3F2FD", "#42A5F5", "#1565C0", "#B71C1C"), 
                        name = "Probability") +
  
  # Layer B: Known Drivers missed by Tier 2
  geom_point(data = filter(plot_df_xgb, Label_Category == "Known_Missed"),
             aes(x = log2FC, y = pmax(mut_count, 1), shape = CNA_Status), 
             color = "#FF9800", size = 3, alpha = 0.8) +
  
  # Define Shapes
  scale_shape_manual(values = c("Amplified" = 16, "Neutral" = 15, "Deleted" = 17),
                     name = "Dominant CNA State") +
  
  # Logarithmic scale for Y-axis
  scale_y_log10(breaks = c(1, 5, 10, 20, 50, 100, 300)) +
  
  # LABELS
  geom_label_repel(data = filter(plot_df_xgb, Label_Category == "Novel_HighProb"),
                   aes(x = log2FC, y = pmax(mut_count, 1), label = gene), 
                   size = 3.5, fontface = "bold", color = "#D32F2F",
                   box.padding = 0.5, point.padding = 0.3, max.overlaps = 50, show.legend = FALSE) +
  
  geom_label_repel(data = filter(plot_df_xgb, Label_Category == "Known_Tier2"),
                   aes(x = log2FC, y = pmax(mut_count, 1), label = gene), 
                   size = 3.5, fontface = "bold", color = "#2E7D32",
                   box.padding = 0.4, point.padding = 0.3, max.overlaps = 50, show.legend = FALSE) +
  
  geom_label_repel(data = filter(plot_df_xgb, Label_Category == "Tier1"),
                   aes(x = log2FC, y = pmax(mut_count, 1), label = paste0(gene, " (Tier 1)")), 
                   size = 3.5, fontface = "italic", color = "black", fill = "#FFF9C4", 
                   box.padding = 0.8, show.legend = FALSE) +
  
  theme_minimal(base_size = 12) +
  labs(title = "Multi-Omics Candidate Landscape of Lung Adenocarcinoma (TCGA-LUAD)",
       subtitle = "XGBoost rescue of novel targets vs. Classical drivers missed by intersection",
       x = "Transcriptional Expression (RNA-seq Unstranded Count log2 Fold Change)",
       y = "Somatic Mutation Count (log10 scale)") +
  theme(plot.title = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "gray80", fill = NA)) +
  
  # Force the dummy legend to render properly by overriding the alpha to 1
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 5, color = NA, alpha = 1)))

# 4. Export
ggsave("Results/Summary/Vis3_MultiOmics_Rescue_Landscape_XGBoost.png", plot = p_landscape_xgb, width = 12, height = 8, dpi = 300)