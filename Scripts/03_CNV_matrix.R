#' =============================================================================
#' Project:  Genomics Project - Phase 1 (LUAD)
#' Task:     Copy Number Variation (CNV) Analysis
#' Author:   Imane BIYAR
#' Date:     July 2026
#'
#' Objective:
#' Parse and merge gene-level CNV data across the TCGA-LUAD cohort.
#' Clean missingness, resolve tibble conversion bugs, and generate 
#' descriptive per-gene features (e.g., Pct_Amplified, Pct_HomDel) 
#' for downstream Machine Learning integration.
#'
#' Inputs:
#'   - Data/CNV/*.tsv (Raw gene-level copy number estimates)
#'   - Results/CNA_segment/cna_list_gistic.rds (GISTIC2.0 significant genes)
#'
#' Outputs:
#'   - Results/CNA/CNV_matrix.rds (Cleaned numeric matrix)
#'   - Results/CNA/cnv_driver_df.rds (Feature dataframe for ML)
#'   - Results/CNA/cnv_list_threshold.rds (Threshold-based significant genes)
#' =============================================================================

library(purrr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

# ==============================================================================
# 3.1 - Build the Gene-Level CNV Matrix
# ==============================================================================

# Execute data gathering only if the matrix does not already exist
if(!file.exists('Results/CNA/CNV_matrix.rds')) {
  
  cnv_files <- list.files('Data/CNV/', pattern = '\\.tsv$', full.names = TRUE, recursive = TRUE)
  
  read_cnv <- function(f) {
    d <- read.table(f, header = TRUE, sep = '\t')
    d_clean <- d %>%
      select(gene_name, copy_number) %>%
      group_by(gene_name) %>%
      summarise(copy_number = mean(copy_number, na.rm = TRUE), .groups = 'drop')
    colnames(d_clean) <- c('gene_name', basename(f))
    return(d_clean)
  }
  
  cnv_matrix <- map(cnv_files, read_cnv) %>% reduce(left_join, by = 'gene_name')
  
  dir.create("Results/CNA", showWarnings = FALSE, recursive = TRUE)
  saveRDS(cnv_matrix, 'Results/CNA/CNV_matrix.rds')
  write.csv(as.data.frame(cnv_matrix), 'Results/CNA/CNV_df.csv', row.names = FALSE)
}

# -----------------------------------------------------------------------------
# 3.2 - Load and set rownames -- TIBBLE FIX
# -----------------------------------------------------------------------------
# BUG FIX: readRDS() returns a tibble. Tibbles do not support rownames.
# This left a character column inside what downstream code treated as a purely numeric matrix. 
# Fixed by converting away from tibble class FIRST, explicitly dropping gene_name, 
# and forcing true numeric storage mode.

cnv_matrix_final <- readRDS("Results/CNA/CNV_matrix.rds")
cnv_matrix_final <- as.data.frame(cnv_matrix_final)     # drop tibble class
rownames(cnv_matrix_final) <- cnv_matrix_final$gene_name
cnv_matrix_final$gene_name <- NULL                       # NOW safe to remove

# Verification
stopifnot(!"gene_name" %in% colnames(cnv_matrix_final))
cnv_matrix_final <- as.matrix(cnv_matrix_final)
storage.mode(cnv_matrix_final) <- "numeric"

cna_matrix <- cnv_matrix_final   

# -----------------------------------------------------------------------------
# 3.3 - Exploratory Data Analysis (EDA)
# -----------------------------------------------------------------------------

# --- Fraction of genome altered per sample ---
frac_altered <- colMeans(cna_matrix != 2, na.rm = TRUE)
ggplot(data.frame(frac = frac_altered), aes(x = frac)) +
  geom_histogram(bins = 40) +
  labs(title = "Fraction of genome CN-altered, per sample", x = "Fraction altered")
ggsave("Results/CNA/genome_altered_per_sample.png", width = 10, height = 7, dpi = 300, units = "in")

# --- Per-gene missingness ---
gene_missing_pct <- rowMeans(is.na(cna_matrix))
ggplot(data.frame(pct = gene_missing_pct), aes(x = pct)) +
  geom_histogram(bins = 50) +
  labs(title = "Per-gene missingness across the CNA cohort", x = "Fraction missing")
ggsave("Results/CNA/missingness.png", width = 10, height = 7, dpi = 300, units = "in")

# -----------------------------------------------------------------------------
# 3.4 - Missingness filter
# -----------------------------------------------------------------------------
cat("Genes with >10% missing:", sum(gene_missing_pct > 0.10), "of", nrow(cna_matrix), "\n")
cnv_matrix_clean <- cna_matrix[gene_missing_pct <= 0.10, ]

cat("Genes remaining after cleaning:", nrow(cnv_matrix_clean), "\n")
cat("Samples remaining:", ncol(cnv_matrix_clean), "\n")   # should read 503

saveRDS(cnv_matrix_clean, file = "Results/CNA/cnv_matrix_clean.rds")

# -----------------------------------------------------------------------------
# 3.5 - Descriptive per-gene CNA features (for the master table / RF features)
# -----------------------------------------------------------------------------
cnv_driver_df <- data.frame(
  Gene          = rownames(cnv_matrix_clean),
  Pct_Amplified = round(rowMeans(cnv_matrix_clean >= 5, na.rm = TRUE) * 100, 1),
  Pct_LowGain   = round(rowMeans(cnv_matrix_clean >= 3 & cnv_matrix_clean <= 4, na.rm = TRUE) * 100, 1),
  Pct_Neutral   = round(rowMeans(cnv_matrix_clean == 2, na.rm = TRUE) * 100, 1),
  Pct_HetDel    = round(rowMeans(cnv_matrix_clean == 1, na.rm = TRUE) * 100, 1),
  Pct_HomDel    = round(rowMeans(cnv_matrix_clean == 0, na.rm = TRUE) * 100, 1),
  Pct_AnyDel    = round(rowMeans(cnv_matrix_clean <= 1, na.rm = TRUE) * 100, 1)
)

# --- GISTIC significance flag (formal, segment-level test) ---
in_CNA_gistic_list <- readRDS("Results/CNA_segment/cna_list_gistic.rds")
cnv_driver_df$in_CNA_gistic <- cnv_driver_df$Gene %in% in_CNA_gistic_list

# --- Dominant event classification (descriptive, >10% recurrence heuristic) ---
cnv_driver_df$Dominant_Event <- case_when(
  cnv_driver_df$Pct_Amplified > 10 ~ "Amplification",
  cnv_driver_df$Pct_AnyDel > 10    ~ "Deletion",
  TRUE                             ~ "Neutral"
)

# --- Secondary independent CNA "significant" list (threshold-based) ---
cnv_list_threshold <- cnv_driver_df$Gene[cnv_driver_df$Dominant_Event != "Neutral"]

cat("CNV threshold-based significant genes (>10% recurrence):", length(cnv_list_threshold), "\n")
cat("GISTIC-based significant genes (formal test):", length(in_CNA_gistic_list), "\n")
cat("Overlap between the two CNA evidence sources:", length(intersect(cnv_list_threshold, in_CNA_gistic_list)), "\n")

saveRDS(cnv_driver_df, file = "Results/CNA/cnv_driver_df.rds")
saveRDS(cnv_list_threshold, file = "Results/CNA/cnv_list_threshold.rds")

# -----------------------------------------------------------------------------
# 3.6 - Tissue-type verification (CNA sample sheet)
# -----------------------------------------------------------------------------
cna_sheet <- read.delim("Data/Metadata_Clinical/gdc_sample_sheet.CNA.tsv", check.names = FALSE)