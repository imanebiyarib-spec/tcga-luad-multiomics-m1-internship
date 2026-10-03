#' =============================================================================
#' Project:  Genomics Project - Phase 1 (LUAD)
#' Task:     MAF Cleaning & dN/dS Driver Significance Testing
#' Author:   Imane BIYAR
#' Date:     July 2026
#'
#' Objective:
#' Merge and clean somatic mutation data (MAF), filter out hypermutator samples,
#' and run the dNdScv algorithm to identify genes under positive evolutionary 
#' selection (driver genes).
#'
#' Inputs:
#'   - Data/MAF/*.maf (Raw somatic mutation files)
#'
#' Outputs:
#'   - Data/MAF/TCGA_LUAD_MAF.rds (Merged raw MAF)
#'   - Results/Mutation/maf_clean.rds (Hypermutators removed)
#'   - Results/Mutation/dndscv_output.rds (Raw dNdScv statistics)
#'   - Results/Mutation/mutation_list_strict.rds (Significant driver genes)
#' =============================================================================

library(maftools)
library(ggplot2)
library(pheatmap)
library(purrr)
library(dplyr)
library(tidyr)
library(tibble)
library(dndscv)

dir.create("Results/Mutation", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 4.1 - Build the Merged MAF Object
# ==============================================================================

# Execute data gathering and integration only if the merged file does not already exist
if(!file.exists("Data/MAF/TCGA_LUAD_MAF.rds")) {
  maf_files <- list.files("Data/MAF/", pattern = "\\.maf$", full.names = TRUE)
  maf <- merge_mafs(maf_files)
  saveRDS(maf, "Data/MAF/TCGA_LUAD_MAF.rds")
}

# -----------------------------------------------------------------------------
# 4.2 - Load and inspect
# -----------------------------------------------------------------------------
maf_obj <- readRDS("Data/MAF/TCGA_LUAD_MAF.rds")
tsb_counts <- getSampleSummary(maf_obj)

ggplot(tsb_counts, aes(x = total)) +
  geom_histogram(bins = 60) +
  scale_x_log10() +
  labs(title = "Mutation burden per sample (log10 scale)")
ggsave("Results/Mutation/Mutation_burden.png", width = 10, height = 7, dpi = 300, units = "in")

png("Results/Mutation/summary_plot.png", width = 2000, height = 1500, res = 200)
plotmafSummary(maf = maf_obj, rmOutlier = TRUE, addStat = "median", dashboard = TRUE)
dev.off()

titv_result <- titv(maf = maf_obj, plot = FALSE, useSyn = TRUE)
png("Results/Mutation/TiTv_plot.png", width = 2000, height = 1500, res = 200)
plotTiTv(res = titv_result)
dev.off()

# -----------------------------------------------------------------------------
# 4.3 - Cleaning: hypermutator exclusion
# -----------------------------------------------------------------------------
png('Results/Mutation/mutation_burden_log.png', width = 1600, height = 1200, res = 150)
hist(log10(tsb_counts$total), breaks = 60,
     main = "Mutation burden per sample (log10)", xlab = "log10(total mutations)")
dev.off()

cutoff <- 800   # set from the histogram above: 3rd quartile ~319, max ~1877 --
# 800 sits well clear of the bulk, catching only the extreme tail
hypermutators <- tsb_counts$Tumor_Sample_Barcode[tsb_counts$total > cutoff]
cat("Excluding", length(hypermutators), "hypermutator samples\n")

maf_clean <- subsetMaf(maf_obj, tsb = setdiff(maf_obj@data$Tumor_Sample_Barcode, hypermutators))
saveRDS(maf_clean, "Results/Mutation/maf_clean.rds")

gene_summary_clean <- getGeneSummary(maf_clean)
min_recurrent_genes <- gene_summary_clean[MutatedSamples >= 3, Hugo_Symbol]
cat("Genes with >=3 mutated samples:", length(min_recurrent_genes), "\n")
saveRDS(min_recurrent_genes, "Results/Mutation/min_recurrent_genes.rds")

png('Results/Mutation/maf_clean_summary.png', width = 1600, height = 1200, res = 150)
plotmafSummary(maf = maf_clean, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()

png('Results/Mutation/oncoplot_maf_clean__top20.png', width = 1400, height = 900, res = 150)
oncoplot(maf = maf_clean, top = 20, showTumorSampleBarcodes = FALSE)
dev.off()

# -----------------------------------------------------------------------------
# 4.4 - dN/dS statistical test (dndscv)
# -----------------------------------------------------------------------------

# Merge maftools data slots:
# By default, maftools separates non-synonymous variants (@data) and synonymous 
# variants (@maf.silent). We MUST combine both so dndscv receives synonymous 
# mutations to calculate its neutral background baseline.

full_mutations <- rbind(maf_clean@data, maf_clean@maf.silent)

dndscv_input <- data.frame(
  sampleID = full_mutations$Tumor_Sample_Barcode, 
  chr      = gsub("chr", "", full_mutations$Chromosome), 
  pos      = full_mutations$Start_Position, 
  ref      = full_mutations$Reference_Allele, 
  mut      = full_mutations$Tumor_Seq_Allele2,
  stringsAsFactors = FALSE
)

# Run the algorithm only if the output does not already exist
if(!file.exists("Results/Mutation/dndscv_output.rds")) {
  dndsout <- dndscv(dndscv_input, refdb = "hg38")
  saveRDS(dndsout, "Results/Mutation/dndscv_output.rds")
}

dndsout <- readRDS("Results/Mutation/dndscv_output.rds")
sel_cv <- dndsout$sel_cv
saveRDS(sel_cv, "Results/Mutation/sel_cv.rds")

# --- Raw q<0.1 list ---
mutation_list_01 <- sel_cv$gene_name[sel_cv$qglobal_cv < 0.1]
cat("Significant driver genes by dN/dS (q<0.1):", length(mutation_list_01), "\n")
saveRDS(mutation_list_01, "Results/Mutation/mutation_list_01.rds")

# --- Tightened list ---
mutation_list_strict <- sel_cv %>%
  filter(qglobal_cv < 0.01,
         wmis_cv > 1 | wnon_cv > 1 | wspl_cv > 1) %>%
  pull(gene_name)

cat("Genes at q<0.1 only:", sum(sel_cv$qglobal_cv < 0.1), "\n")
cat("Genes at q<0.01 only:", sum(sel_cv$qglobal_cv < 0.01), "\n")
cat("Genes at q<0.01 AND w>1:", length(mutation_list_strict), "\n")
saveRDS(mutation_list_strict, "Results/Mutation/mutation_list_strict.rds")