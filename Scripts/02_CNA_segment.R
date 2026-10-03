#' =============================================================================
#' Project:  Genomics Project - Phase 1 (LUAD)
#' Task:     CNA Segment Formatting & GISTIC2.0 Parsing
#' Author:   Imane BIYAR
#' Date:     July 2026
#'
#' Objective:
#' Build a Master Segment File (.seg) for GISTIC2.0 from 503 TCGA ASCAT Folders.
#' Parse the local GISTIC2.0 Docker output (-twosides 1) to identify genes 
#' undergoing significant focal amplification or deletion.
#'
#' Inputs:
#'   - Data/CNA_segment/ (Raw TCGA ASCAT segment files)
#'   - Results/CNA_segment/amp_genes.conf_99.txt (GISTIC2.0 output)
#'   - Results/CNA_segment/del_genes.conf_99.txt (GISTIC2.0 output)
#'
#' Outputs:
#'   - Results/CNA_segment/gistic_input.seg
#'   - Results/CNA_segment/segment_per_sample.png (QC histogram)
#'   - Results/CNA_segment/cna_list_gistic.csv & .rds (Significant CNA drivers)
#' =============================================================================

library(dplyr)

# -----------------------------------------------------------------------------
# 2.1 - Build the .seg file 
# -----------------------------------------------------------------------------
# Execute data gathering and integration only if the segment file does not already exist
if(!file.exists("Results/CNA_segment/gistic_input.seg")) {
  
  seg_files <- list.files(path = "Data/CNA_segment",
                          pattern = "^TCGA.*\\.(txt|tsv|seg)$",
                          full.names = TRUE, recursive = TRUE)
  seg_files <- seg_files[!grepl("annotations|logs", seg_files, ignore.case = TRUE)]
  
  read_and_format_seg <- function(file_path) {
    df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    safe_cn <- ifelse(df$Copy_Number == 0, 0.1, df$Copy_Number)   # avoid log2(0)
    seg_mean <- log2(safe_cn / 2)
    est_markers <- pmax(10, round((df$End - df$Start) / 5000))
    data.frame(
      Sample       = df$GDC_Aliquot,
      Chromosome   = gsub("chr", "", df$Chromosome),
      Start        = df$Start,
      End          = df$End,
      Num_Markers  = est_markers,
      Segment_Mean = round(seg_mean, 4),
      stringsAsFactors = FALSE
    )
  }
  
  gistic_list <- lapply(seg_files, read_and_format_seg)
  master_seg  <- do.call(rbind, gistic_list)
  master_seg <- subset(master_seg, !Chromosome %in% c("X", "Y", "chrX", "chrY"))
  
  dir.create("Results/CNA_segment", showWarnings = FALSE, recursive = TRUE)
  write.table(master_seg, "Results/CNA_segment/gistic_input.seg",
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# -----------------------------------------------------------------------------
# 2.2 - Segment-count QC (CNA equivalent of the MAF hypermutator check)
# -----------------------------------------------------------------------------
master_seg <- read.delim("Results/CNA_segment/gistic_input.seg", header = TRUE)

sample_seg_counts <- master_seg %>%
  group_by(Sample) %>%
  summarise(Num_Segments = n())

write.table(sample_seg_counts, "Results/CNA_segment/sample_seg_counts",
            sep = "\t", row.names = FALSE, quote = FALSE)

png("Results/CNA_segment/segment_per_sample.png", width = 2000, height = 1500, res = 200)
hist(sample_seg_counts$Num_Segments, breaks = 40, col = "#1976D2", border = "white",
     main = "Number of CNA segments per sample",
     xlab = "Segment count (Breakpoints)", ylab = "Number of Samples")
dev.off()

##########################################################################################################
# 2.3 - Parse Local GISTIC2.0 Output (Ubuntu Docker run, -twosides 1)
##########################################################################################################
amp_file <- "Results/CNA_segment/amp_genes.conf_99.txt"
del_file <- "Results/CNA_segment/del_genes.conf_99.txt"

# GISTIC formats these files with columns as peaks and rows as metadata/genes:
# Row 1: Cytoband | Row 2: q-value | Row 3: Residual q-value |
# Row 4: Wide peak boundaries | Row 5+: Gene symbols
extract_gistic_genes <- function(file_path, q_threshold = 0.1) {
  
  if(!file.exists(file_path)) return(character(0))
  
  df <- read.delim(file_path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  q_values <- as.numeric(df[2, -1])
  sig_cols <- which(q_values < q_threshold)
  if (length(sig_cols) == 0) {
    return(character(0))
  }
  raw_genes <- unlist(df[5:nrow(df), sig_cols + 1])   # genes start at row 5
  clean_genes <- unique(raw_genes[raw_genes != "" & !is.na(raw_genes)])
  return(clean_genes)
}

amp_genes <- extract_gistic_genes(amp_file, q_threshold = 0.1)
del_genes <- extract_gistic_genes(del_file, q_threshold = 0.1)
cna_list_gistic <- unique(c(amp_genes, del_genes))

# OPTIONAL: try a stricter q-value if Section 9 (Tier1 reduction) needs it
amp_genes_strict <- extract_gistic_genes(amp_file, q_threshold = 0.05)
del_genes_strict <- extract_gistic_genes(del_file, q_threshold = 0.05)
cna_list_gistic_strict <- unique(c(amp_genes_strict, del_genes_strict))

write.csv(data.frame(gene_symbol = cna_list_gistic),
          "Results/CNA_segment/cna_list_gistic.csv", row.names = FALSE)
saveRDS(cna_list_gistic, file = "Results/CNA_segment/cna_list_gistic.rds")

write.csv(data.frame(gene_symbol = cna_list_gistic_strict),
          "Results/CNA_segment/cna_list_gistic_strict.csv", row.names = FALSE)
saveRDS(cna_list_gistic_strict, file = "Results/CNA_segment/cna_list_gistic_strict.rds")