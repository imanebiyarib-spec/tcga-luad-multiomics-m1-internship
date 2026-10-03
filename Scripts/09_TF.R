#' =============================================================================
#' Project:  Genomics Project - Phase 1 (LUAD)
#' Task:     ChIP-seq Regulatory Analysis Pipeline
#' Author:   Imane BIYAR
#' Date:     August 2026
#'
#' Objective:
#' Annotate FOXA1 physical binding sites (ChIP-seq data from A549 cells) and 
#' intersect them with the Top 10 novel machine learning candidates. 
#' Evaluate potential dual-activation mechanisms by testing for concurrent 
#' physical binding, structural amplification, and positive transcriptional 
#' correlation with the FOXA1 regulator.
#'
#' Inputs:
#'   - Data/ChIPseq/FOXA1_A549_peaks.bed
#'   - Results/Integration/master_table_pathways.rds
#'   - Results/Integration/top_10_novel_for_literature_validation_XGBoost.csv
#'   - Results/Integration/top_10_novel_for_literature_validation.csv
#'   - Results/DEG/TCGA_count_matrix.csv
#'
#' Outputs:
#'   - Results/ChIPseq/peak_annotation_pie_custom.png
#'   - Results/ChIPseq/ChIPseq_promoter_genes.csv
#'   - Results/ChIPseq/FOXA1_ML_Novel_Intersection.csv
#'   - Results/ChIPseq/Novel_ML_Dual_Activation.csv
#'   - Results/ChIPseq/venn_ml_novel_integration.png
#'   - Results/ChIPseq/FOXA1_Expression_Correlation_Validation.csv
#' =============================================================================

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggplot2)
library(ggVennDiagram)
library(AnnotationDbi) 

output_dir <- "Results/ChIPseq/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 9.1 Annotate Peaks to Genes
# ==============================================================================

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peaks <- readPeakFile('Data/ChIPseq/FOXA1_A549_peaks.bed')

anno <- annotatePeak(
  peaks,
  tssRegion = c(-3000, 3000),
  TxDb = txdb,
  annoDb = 'org.Hs.eg.db'
)

# Extract statistics and generate high-resolution plot
my_colors <- c("Promoter (<=1kb)" = "darkred", "Promoter (1-2kb)" = "blue", 
               "Promoter (2-3kb)" = "darkblue", "5' UTR" = "gold", 
               "3' UTR" = "darkgoldenrod", "1st Exon" = "green", 
               "Other Exon" = "darkgreen", "1st Intron" = "pink", 
               "Other Intron" = "hotpink", "Downstream (<=300)" = "orange", 
               "Distal Intergenic" = "darkgray")

anno_stat <- as.data.frame(anno@annoStat)
p <- ggplot(anno_stat, aes(x = "", y = Frequency, fill = Feature)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = my_colors) + 
  labs(title = "FOXA1 Binding Site Distribution", fill = "Genomic Feature") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.text = element_text(size = 10))

png(paste0(output_dir, 'peak_annotation_pie_custom.png'), width = 1000, height = 700, res = 150)
print(p)
dev.off()

# Extract genes with TF binding peaks at promoters for downstream analysis
anno_df <- as.data.frame(anno)
promoter_genes <- unique(anno_df$SYMBOL[grepl('Promoter', anno_df$annotation)])
write.csv(data.frame(gene=promoter_genes), paste0(output_dir, 'ChIPseq_promoter_genes.csv'), row.names=FALSE)

# ==============================================================================
# 9.2 Overlap ChIP Targets with Top 10 Novel ML Genes
# ==============================================================================

# Load the integrated master table for DEG/CNA context
master_table <- readRDS('Results/Integration/master_table_pathways.rds')

deg_driver_genes <- master_table$gene[master_table$in_DEG == TRUE]
cnv_amp_genes    <- master_table$gene[master_table$Pct_Amplified > 10]

# Load the Novel ML Candidates
top10_xgb <- read.csv('Results/Integration/top_10_novel_for_literature_validation_XGBoost.csv')$gene
top10_rf  <- read.csv('Results/Integration/top_10_novel_for_literature_validation.csv')$gene

combined_ml_novel <- unique(c(top10_xgb, top10_rf))

# Perform Overlaps with FOXA1 physical binding
chip_xgb_novel      <- intersect(promoter_genes, top10_xgb)
chip_rf_novel       <- intersect(promoter_genes, top10_rf)
chip_combined_novel <- intersect(promoter_genes, combined_ml_novel)

write.csv(data.frame(gene=chip_combined_novel), paste0(output_dir, 'FOXA1_ML_Novel_Intersection.csv'), row.names=FALSE)

# ==============================================================================
# 9.3: CNV + ChIP-seq Combined Interpretation for Novel Genes
# ==============================================================================

# Filter the combined ML novel pool by CNV amplification status
novel_amp_genes <- intersect(combined_ml_novel, cnv_amp_genes)

# Calculate Triple Overlap: ML Novel + CNV Amplified + FOXA1 Bound
novel_dual_activation <- intersect(novel_amp_genes, promoter_genes)

write.csv(data.frame(gene=novel_dual_activation, mechanism='Novel_ML_CNV_amp_plus_FOXA1'), 
          paste0(output_dir, 'Novel_ML_Dual_Activation.csv'), row.names=FALSE)

# ==============================================================================
# 9.4: Four-Way Venn Diagram 
# ==============================================================================

venn_data_ml <- list(
  'ChIP-seq Targets' = promoter_genes,     
  'XGBoost Novel 10' = top10_xgb,       
  'RF Novel 10'      = top10_rf,
  'DEGs'             = deg_driver_genes      
)

venn_plot_ml <- ggVennDiagram(venn_data_ml, label_alpha = 0) +
  scale_fill_gradient(low = '#D6EEF1', high = '#1D6A72') +
  ggtitle('Integration: FOXA1 Targets vs. Novel ML Candidates') +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

ggsave(paste0(output_dir, 'venn_ml_novel_integration.png'), plot = venn_plot_ml, width = 10, height = 8, dpi = 300)

# ==============================================================================
# 9.5: Mechanistic Validation (Correlation Analysis on Novel Genes)
# ==============================================================================

exp_matrix <- read.csv('Results/DEG/TCGA_count_matrix.csv', row.names=1, check.names=FALSE)

clean_ids <- gsub("\\..*", "", rownames(exp_matrix))
gene_symbols <- mapIds(org.Hs.eg.db, keys = clean_ids, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")

valid_idx <- !is.na(gene_symbols)
exp_matrix_clean <- exp_matrix[valid_idx, ]
exp_matrix_clean$Hugo_Symbol <- gene_symbols[valid_idx]

exp_matrix_final <- aggregate(exp_matrix_clean[, -ncol(exp_matrix_clean)], 
                              by = list(exp_matrix_clean$Hugo_Symbol), FUN = mean)
rownames(exp_matrix_final) <- exp_matrix_final$Group.1
exp_matrix_final$Group.1 <- NULL

if("FOXA1" %in% rownames(exp_matrix_final)) {
  foxa1_expr <- as.numeric(exp_matrix_final["FOXA1", ])
  
  cor_scores <- sapply(chip_combined_novel, function(gene) {
    if(gene %in% rownames(exp_matrix_final)) {
      cor(foxa1_expr, as.numeric(exp_matrix_final[gene, ]), method = "spearman")
    } else { NA }
  })
  
  validated_novel_03 <- chip_combined_novel[!is.na(cor_scores) & cor_scores > 0.3]
  
  if(length(validated_novel_03) > 0) {
    # Consolidate correlation results into a clean dataframe instead of raw prints
    validation_df <- data.frame(
      Gene = names(cor_scores[!is.na(cor_scores) & cor_scores > 0.3]),
      Spearman_Correlation_with_FOXA1 = as.numeric(cor_scores[!is.na(cor_scores) & cor_scores > 0.3])
    )
    validation_df <- validation_df[order(-validation_df$Spearman_Correlation_with_FOXA1), ]
    write.csv(validation_df, paste0(output_dir, 'FOXA1_Expression_Correlation_Validation.csv'), row.names = FALSE)
  }
} else {
  warning("FOXA1 not found in expression matrix.")
}