#' =============================================================================
#' Project:  Genomics Project - Phase 1 (LUAD)
#' Task:     Multi-Omics Integration (Master Matrix)
#' Author:   Imane BIYAR
#' Date:     July 2026
#'
#' Objective:
#' Integrate the three omics layers (Transcriptomics, Copy Number Alterations, 
#' and Somatic Mutations) into a single master matrix. Define Tier 1 (strict 
#' intersection) and Tier 2 (>= 2 layer union) candidate gene pools for ML.
#'
#' Inputs:
#'   - Results/DEG/res_df_clean.rds & sig.rds
#'   - Results/CNA/cnv_driver_df.rds & cnv_list_threshold.rds
#'   - Results/CNA_segment/cna_list_gistic.rds
#'   - Results/Mutation/maf_clean.rds, sel_cv.rds, mutation_list_strict.rds
#'
#' Outputs:
#'   - Results/Integration/master_table.rds (Integrated feature backbone)
#'   - Results/Integration/gene_lists_for_venn.rds
#'   - Results/Integration/venn_3way.png & upset_4way.png
#'   - Results/DEG/deg_threshold_sensitivity.csv
#' =============================================================================

library(dplyr)
library(tidyr)
library(ggVennDiagram)
library(UpSetR)
library(ggplot2)

dir.create("Results/Integration", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 5.1 - Load all layer outputs
# ==============================================================================
res_clean  <- readRDS("Results/DEG/res_df_clean.rds")      
sig        <- readRDS("Results/DEG/sig.rds")               

cnv_driver_df      <- readRDS("Results/CNA/cnv_driver_df.rds")
cnv_list_threshold <- readRDS("Results/CNA/cnv_list_threshold.rds")

cna_list_gistic    <- readRDS("Results/CNA_segment/cna_list_gistic.rds")

maf_clean <- readRDS("Results/Mutation/maf_clean.rds")
sel_cv    <- readRDS("Results/Mutation/sel_cv.rds")
mutation_list <- readRDS("Results/Mutation/mutation_list_strict.rds")

# -----------------------------------------------------------------------------
# 5.2 - Define each layer's gene universe (for the UNION backbone)
# -----------------------------------------------------------------------------
rna_genes <- rownames(res_clean)          
cna_genes <- cnv_driver_df$Gene           
mut_genes <- unique(maf_clean@data$Hugo_Symbol)

all_genes <- union(rna_genes, union(cna_genes, mut_genes))

# NOTE: This is expected to be large (dominated by cna_genes' size). This is the
# UNION of gene universes, not the intersection. Most rows will have genuine
# "neutral" values in one or two layers; that is correct, not a bug.

# -----------------------------------------------------------------------------
# 5.3 - Per-layer statistic tables
# -----------------------------------------------------------------------------
rna_df <- data.frame(gene = rownames(res_clean),
                     log2FC = res_clean$log2FoldChange,
                     padj   = res_clean$padj,
                     stat   = res_clean$stat)

cna_df <- cnv_driver_df %>% dplyr::rename(gene = Gene)   

n_samples_tested <- length(unique(maf_clean@data$Tumor_Sample_Barcode))

mut_df <- sel_cv %>%
  transmute(gene = gene_name,
            mut_count     = n_mis + n_non + n_spl,
            mut_freq_norm = mut_count / n_samples_tested,
            mut_q         = qglobal_cv)

# -----------------------------------------------------------------------------
# 5.4 - Build the master table
# -----------------------------------------------------------------------------
master_table <- data.frame(gene = all_genes) %>%
  left_join(rna_df, by = "gene") %>%
  left_join(cna_df, by = "gene") %>%
  left_join(mut_df, by = "gene") %>%
  distinct(gene, .keep_all = TRUE)

master_table <- master_table %>%
  mutate(
    log2FC          = ifelse(is.na(log2FC), 0, log2FC),
    padj            = ifelse(is.na(padj), 1, padj),
    Pct_Amplified   = ifelse(is.na(Pct_Amplified), 0, Pct_Amplified),
    Pct_HomDel      = ifelse(is.na(Pct_HomDel), 0, Pct_HomDel),
    in_CNA_gistic   = ifelse(is.na(in_CNA_gistic), FALSE, in_CNA_gistic),
    mut_count       = ifelse(is.na(mut_count), 0, mut_count),
    mut_freq_norm   = ifelse(is.na(mut_freq_norm), 0, mut_freq_norm),
    mut_q           = ifelse(is.na(mut_q), 1, mut_q)
  )

# NOTE on imputation choice: 0/1 neutral-value imputation is used deliberately
# rather than mean/KNN imputation. This is MNAR (missing not at random) data
# in a strong sense -- a gene absent from the mutation table is not a value we
# failed to measure, it IS zero recurrent mutations.

# -----------------------------------------------------------------------------
# 5.5 - DEG threshold sensitivity (transparency funnel)
# -----------------------------------------------------------------------------
thresholds_to_try <- list(
  current       = list(padj = 0.01,  fc = 1),
  stricter_fdr  = list(padj = 0.001, fc = 1),
  stricter_fc   = list(padj = 0.01,  fc = 2),
  both_stricter = list(padj = 0.001, fc = 2)
)

funnel <- lapply(names(thresholds_to_try), function(nm) {
  t <- thresholds_to_try[[nm]]
  n <- sum(res_clean$padj < t$padj & abs(res_clean$log2FoldChange) > t$fc, na.rm = TRUE)
  data.frame(threshold = nm, padj_cutoff = t$padj, fc_cutoff = t$fc, n_DEG = n)
})
funnel_table <- do.call(rbind, funnel)
write.csv(funnel_table, "Results/DEG/deg_threshold_sensitivity.csv", row.names = FALSE)

# CHOSEN threshold (principled: mirrors the q<0.01 & w>1 combined significance
# + effect-size logic already applied to the mutation layer):
rna_list <- res_clean$display_name[res_clean$padj < 0.001 & abs(res_clean$log2FoldChange) > 2]

# -----------------------------------------------------------------------------
# 5.6 - The three (plus one comparison) significant-gene lists
# -----------------------------------------------------------------------------
cna_list <- cna_list_gistic        # PRIMARY CNA evidence -- formal, segment-level test q<0.1
mut_list <- mutation_list          # tightened dN/dS list (q<0.01 & w>1)

# --- Tier 1 / Tier 2, using GISTIC as the primary CNA evidence ---
tier1_genes <- Reduce(intersect, list(rna_list, cna_list, mut_list))

all_candidate_genes <- union(rna_list, union(cna_list, mut_list))
layer_counts <- sapply(all_candidate_genes, function(g)
  sum(g %in% rna_list, g %in% cna_list, g %in% mut_list))
tier2_genes <- names(layer_counts)[layer_counts >= 2]

# --- Robustness comparison: using the CNV-threshold list instead of GISTIC ---
tier1_genes_cnv <- Reduce(intersect, list(rna_list, cnv_list_threshold, mut_list))

# -----------------------------------------------------------------------------
# 5.7 - Attach boolean columns to master_table (GISTIC-based, primary)
# -----------------------------------------------------------------------------
master_table <- master_table %>%
  mutate(
    in_DEG   = gene %in% rna_list,
    in_CNA   = gene %in% cna_list,
    in_CNA_cnv_threshold = gene %in% cnv_list_threshold,   
    in_MUT   = gene %in% mut_list,
    n_layers = in_DEG + in_CNA + in_MUT,
    in_tier1 = n_layers == 3,
    in_tier2 = n_layers >= 2
  )

saveRDS(list(rna = rna_list, cna = cna_list, mut = mut_list,
             cnv_threshold = cnv_list_threshold,
             tier1 = tier1_genes, tier2 = tier2_genes,
             tier1_cnv_threshold = tier1_genes_cnv),
        "Results/Integration/gene_lists_for_venn.rds")

saveRDS(master_table, "Results/Integration/master_table.rds")

##########################################################################################################
# 5.8 - Venn diagram (3-way, primary) and UpSet plot (4-way, robustness)
##########################################################################################################
gene_lists <- readRDS("Results/Integration/gene_lists_for_venn.rds")

venn_input <- list(
  "RNA-seq DEG"      = gene_lists$rna,
  "CNA (GISTIC)"     = gene_lists$cna,
  "Mutation (dN/dS)" = gene_lists$mut
)

p_venn <- ggVennDiagram(venn_input, label = "count") +
  scale_fill_gradient(low = "white", high = "#1976D2") +
  ggtitle("Overlap of Significant Genes Across Omics Layers") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Results/Integration/venn_3way.png", plot = p_venn, width = 8, height = 8, dpi = 300)

# 4-way UpSet plot -- UpSetR is the standard tool for 4+ set comparisons
upset_input <- list(
  "RNA-seq DEG"     = gene_lists$rna,
  "CNA (GISTIC)"    = gene_lists$cna,
  "CNV (threshold)" = gene_lists$cnv_threshold,
  "Mutation (dN/dS)"= gene_lists$mut
)

png("Results/Integration/upset_4way.png", width = 2000, height = 1200, res = 200)
upset(fromList(upset_input), order.by = "freq", nsets = 4)
dev.off()