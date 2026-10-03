#' =============================================================================
#' Project:  Genomics Project - Phase 1 (LUAD)
#' Task:    Differential Expression Analysis (DEA)
#' Author:  Imane BIYAR
#' Date:    June 2026
#'
#' Objective:
#' Identify dysregulated genes (DEGs) between LUAD tumor and healthy samples.
#' Controls for covariates (Smoking, Gender) and retains the FULL
#' genome ranking to feed into Gene Set Enrichment Analysis (GSEA) and Phase 2.
#'
#' Inputs:
#'   - Results/DEG/TCGA_count_matrix.csv (Merged counts)
#'   - Metadata_DESeq2_Clean.csv (Covariate metadata)
#'
#' Outputs:
#'   - Results/DEG/DEG_results_FULL.csv (The Master list for GSEA and ML)
#'   - Results/DEG/volcano_plot.png & heatmap_top50.png
#' =============================================================================

# 3.1 - Build the TCGA Count Matrix

# Loading the tools
library(ggpubr)
library(edgeR)
library(DESeq2)
library(dplyr)
library(org.Hs.eg.db)
library(ggplot2)
library(matrixStats)
library(ggrepel)
library(EnhancedVolcano)
library(pheatmap)

# Execute data gathering and integration only if the matrix does not already exist
if(!file.exists('Results/DEG/TCGA_count_matrix.rds')) {
  
  # Gathering the files
  files <- list.files('Data/TCGA_counts/', pattern='\\.tsv$', full.names=TRUE, recursive=TRUE)
  
  # Extraction
  # The unstranded column is the GDC standard for "Harmonized Counts."
  read_one <- function(f) {
    # 1. Read the file skipping only the first metadata line (#)
    d <- read.table(f, header=TRUE, skip=1, sep='\t', check.names=FALSE)
    
    # 2. Remove the first 4 rows (N_unmapped, etc.) which are not genes
    d <- d[-(1:4), ]
    
    # 3. Extract columns and rename 'unstranded' to the filename for the matrix
    d_clean <- d[, c('gene_id', 'unstranded')]
    colnames(d_clean) <- c('gene_id', basename(f))
    
    return(d_clean)
  }
  
  # Integration
  count_list <- lapply(files, read_one) 
  count_matrix <- Reduce(function(a,b) merge(a, b, by='gene_id'), count_list)
  
  # Final formatting & serialization
  rownames(count_matrix) <- count_matrix$gene_id
  count_matrix <- count_matrix[, -1] 
  
  saveRDS(count_matrix, 'Results/DEG/TCGA_count_matrix.rds')
  write.csv(count_matrix, 'Results/DEG/TCGA_count_matrix.csv', row.names=TRUE)
}

##########################################################################################################

# 3.2 - Load Covariate Metadata

metadata <- read.csv('Results/Metadata/metadata.csv', row.names = 2)

# Extract necessary columns
metadata_short <- data.frame(
  case_id     = metadata$Case.ID,
  sample_id   = metadata$Sample.ID,
  sample_type = metadata$Tissue.Type,
  smoking_status = metadata$exposures.tobacco_smoking_status,
  gender = metadata$demographic.gender,
  stringsAsFactors = FALSE
)
rownames(metadata_short) <- rownames(metadata)

# New variable 'Condition' (Tumor vs Healthy)
metadata_short$condition <- ifelse(grepl("Normal", metadata_short$sample_type), "Healthy", "Tumor")

# Ensure factors are set for the statistical model, 'Healthy' as a reference
metadata_short$condition <- factor(metadata_short$condition, levels = c("Healthy", "Tumor"))
metadata_short$smoking_status <- factor(metadata_short$smoking_status)
metadata_short$gender <- factor(metadata_short$gender)

# Load the count matrix
count_matrix <- read.csv('Results/DEG/TCGA_count_matrix.csv', header=TRUE, check.names=FALSE, row.names = 1)

# Convert to the canonical object used throughout this document
counts_raw <- readRDS("Results/DEG/TCGA_count_matrix.rds")
counts <- as.matrix(counts_raw)
storage.mode(counts) <- "integer"

# --- Library size per sample ---
lib_sizes <- colSums(counts)

ggplot(data.frame(sample = names(lib_sizes), lib_size = lib_sizes),
       aes(x = reorder(sample, lib_size), y = lib_size)) +
  geom_col() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Library size per sample", y = "Total reads", x = "Sample (sorted)")
ggsave("Results/DEG/library_size_per_sample.png", width = 10, height = 7, dpi = 300, units = "in")

# Organize lines correctly for DESeq2
metadata_short <- metadata_short[match(colnames(count_matrix), rownames(metadata_short)), ]
rownames(metadata_short) <- colnames(count_matrix)

##########################################################################################################

# 3.3 - Filter with edgeR & Run DESeq2 (The Covariate Model)

# 1. Rigorous filtering with edgeR::filterByExpr BEFORE building the DESeq2 object
group <- metadata_short$condition   
keep <- filterByExpr(count_matrix, group = group)

cat("Genes kept by filterByExpr:", sum(keep), "of", length(keep), "\n")
count_matrix_filtered <- count_matrix[keep, ]

# 2. Build and run the DESeqDataSet object (cached to avoid redundant 20-minute execution)
if(!file.exists("Results/DEG/dds_processed_covariates.rds")) {
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix_filtered,
    colData   = metadata_short,
    design    = ~ smoking_status + gender + condition
  )
  
  dds <- DESeq(dds)
  saveRDS(dds, file = "Results/DEG/dds_processed_covariates.rds")
}

dds <- readRDS("Results/DEG/dds_processed_covariates.rds")

# 3. Extract results 
if(!file.exists("Results/DEG/res_dds.rds")) {
  res <- results(dds, 
                 contrast = c('condition', 'Tumor', 'Healthy'), 
                 independentFiltering = FALSE, # Disabled to retain all ranked genes for GSEA!
                 alpha = 0.01,
                 pAdjustMethod = 'BH')
  
  saveRDS(res, file = "Results/DEG/res_dds.rds")
}

##########################################################################################################

# 3.4 - Annotate and CLEAN the FULL Gene List for Multi-Omics Integration

res <- readRDS("Results/DEG/res_dds.rds")
res_clean <- res[res$baseMean > 0, ]
res_df <- as.data.frame(res_clean)

# 0. Save original versioned Ensembl IDs (CRITICAL bridge for linking back to dds/vsd in plots!)
res_df$ensembl_versioned <- rownames(res_df)

# 1. Clean Ensembl IDs (strip version numbers)
clean_ids <- gsub("\\..*", "", rownames(res_df))

# 2. Map Ensembl IDs to official HGNC Gene Symbols
res_df$symbol <- mapIds(org.Hs.eg.db,
                        keys = clean_ids,
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

res_df$ensembl_id <- clean_ids

# ----------------- CRITICAL CLEANING FOR THE BIG MATRIX -----------------

# A) Remove genes that lack an official symbol
res_df_clean <- res_df[!is.na(res_df$symbol), ]

# B) Deduplicate Gene Symbols: keep only the most significant one
res_df_clean <- res_df_clean[order(res_df_clean$padj, -res_df_clean$baseMean), ]
res_df_clean <- res_df_clean[!duplicated(res_df_clean$symbol), ]

# C) Set clean display names and make gene symbols the rownames
res_df_clean$display_name <- res_df_clean$symbol
rownames(res_df_clean) <- res_df_clean$symbol

cat("Initial genes in DESeq2:", nrow(res_df), "\n")
cat("Unique genes in clean Big Matrix:", nrow(res_df_clean), "\n")

# 3. SAVE THE CLEAN BIG MATRIX
write.csv(res_df_clean, 'Results/DEG/res_df_clean.csv') 
saveRDS(res_df_clean, file = "Results/DEG/res_df_clean.rds")

# -----------------------------------------------------------------------------
# 4. EXTRACT THE RNA-SEQ SPECIFIC LIST (Significant DEGs only)
sig <- subset(res_df_clean, padj < 0.01 & abs(log2FoldChange) > 1) 
sig <- sig[order(sig$padj), ] 
cat("Number of significant DEGs retained:", nrow(sig), "\n")

write.csv(sig, 'Results/DEG/sig.csv')
saveRDS(sig, file = "Results/DEG/sig.rds")

# Extract Top 20 Up and Down from the cleaned dataset
top20_up   <- head(sig[sig$log2FoldChange > 0, ], 20)
top20_down <- head(sig[sig$log2FoldChange < 0, ], 20)

top20 <- data.frame(
  Upregulated   = top20_up$display_name,
  Downregulated = top20_down$display_name
)
print(top20)

##########################################################################################################

# 3.5 Visualizations (Generated strictly from the CLEANED dataset)

### Exploratory Data Analysis
# ------------------------------Dispersion Plot----------------------------------------------------------
png("Results/DEG/dispersion_Top20_plot.png", width=2000, height=1500, res=200)
plotDispEsts(dds)

# Use our saved versioned Ensembl IDs to label the top genes on the dispersion plot
target_ids    <- c(top20_up$ensembl_versioned, top20_down$ensembl_versioned)
target_labels <- c(top20_up$display_name, top20_down$display_name)

with(mcols(dds)[target_ids, ], {
  text(baseMean, dispGeneEst,
       labels = target_labels, pos = 4, offset = 0.8, cex = 0.9, font = 2, col = "black")
})
dev.off()


### --------------------------------------PCA-----------------------------------------------------------
vsd <- vst(dds, blind = FALSE)

# FILTER VST: We restrict the PCA strictly to our cleaned, deduplicated gene set!
vsd_clean <- vsd[res_df_clean$ensembl_versioned, ]

ntop <- 500
rv <- rowVars(assay(vsd_clean))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t(assay(vsd_clean)[select, ])

pca_res <- prcomp(mat, scale. = FALSE)

df_scores <- as.data.frame(pca_res$x)
df_scores$group <- vsd_clean$condition

# Map loadings directly to our clean display names
df_loadings <- as.data.frame(pca_res$rotation)
df_loadings$display_name <- res_df_clean$display_name[match(rownames(df_loadings), res_df_clean$ensembl_versioned)]

top_genes <- df_loadings[order(abs(df_loadings$PC1), decreasing = TRUE), ][1:20, ]

scale_factor <- max(abs(df_scores$PC1)) / max(abs(top_genes$PC1)) * 0.7
top_genes$PC1_scaled <- top_genes$PC1 * scale_factor
top_genes$PC2_scaled <- top_genes$PC2 * scale_factor

ggplot() +
  geom_point(data = df_scores, aes(x = PC1, y = PC2, color = group), size = 3, alpha = 1) +
  geom_segment(data = top_genes, aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.5) +
  geom_text(data = top_genes, aes(x = PC1_scaled * 1.1, y = PC2_scaled * 1.1, label = display_name),
            color = "black", fontface = "bold", check_overlap = TRUE) +
  scale_color_manual(values = c("Tumor" = "#FF5722", "Healthy" = "#1976D2")) +
  theme_bw() +
  labs(
    title = "PCA Biplot: Sample Clusters and Driving Genes (Cleaned Cohort)",
    x = paste0("PC1: ", round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2))[1]), "% variance"),
    y = paste0("PC2: ", round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2))[2]), "% variance"),
    color = "Condition"
  )

ggsave("Results/DEG/pca_biplot.png", width = 10, height = 7, dpi = 300, units = "in")


# -----------------------------------MA-Plot ----------------------------------------------
p_ma <- ggmaplot(res_df_clean,
                 main = "MA Plot: LUAD Tumor vs Healthy (Cleaned Cohort)",
                 fdr = 0.01, fc = 2, size = 0.5,
                 palette = c("#FF5722", "#1976D2", "darkgray"),
                 genenames = as.character(res_df_clean$display_name),
                 top = 20, font.label = c(12, "bold", "black"),
                 label.rectangle = TRUE, font.main = "bold",
                 ggtheme = theme_bw())

ggsave("Results/DEG/ggma_plot.png", plot = p_ma, width = 10, height = 8, dpi = 300)


# ----------------------------------Volcano Plot -----------------------------------------------
genes_to_show <- c(top20_up$display_name, top20_down$display_name)

EnhancedVolcano(res_df_clean,
                lab = res_df_clean$display_name,
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = genes_to_show,
                pCutoff = 0.01,
                FCcutoff = 1.0,
                xlim = c(-10, 15),
                ylim = c(0, 250),
                labSize = 4.0,
                pointSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                lengthConnectors = unit(0.01, "npc"),
                typeConnectors = "closed",
                endsConnectors = "last",
                colConnectors = "grey30",
                title = 'TCGA-LUAD: Transcriptomic Landscape',
                subtitle = 'Differential Expression: Tumor vs Healthy (Covariate Adjusted)',
                caption = paste0('Total = ', nrow(res_df_clean), 'genes'),
                legendPosition = 'bottom',
                legendLabSize = 10,
                axisLabSize = 12,
                gridlines.major = FALSE,
                gridlines.minor = FALSE)

ggsave('Results/DEG/volcano_plot.png', width=12, height=10, dpi=300, bg="white")


# --------------------------------------Heatmap of Top 50 DEGs-----------------------------------
res_viz <- res_df_clean[res_df_clean$baseMean > 50, ]
top50_indices <- head(order(res_viz$padj, na.last = NA), 50)

# Pull normalized expression using the versioned Ensembl ID bridge!
mat <- assay(vsd)[res_viz$ensembl_versioned[top50_indices], ]
mat <- mat - rowMeans(mat)

# Label rows cleanly with gene symbols
rownames(mat) <- res_viz$display_name[top50_indices]

anno <- as.data.frame(colData(vsd)[, 'condition', drop = FALSE])

pheatmap(mat,
         annotation_col = anno,
         show_rownames  = TRUE,
         show_colnames  = FALSE,
         cluster_rows   = TRUE,
         cluster_cols   = TRUE,
         color          = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main           = "Top 50 Differentially Expressed Genes in LUAD (Cleaned)",
         annotation_colors = list(Condition = c(Healthy = "#1976D2", Tumor = "#FF5722")),
         filename       = 'Results/DEG/heatmap_top50_labeled.png',
         width=12, height=10)