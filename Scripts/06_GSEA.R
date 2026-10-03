#' =============================================================================
#' Project:  Genomics Project - Phase 1 (LUAD)
#' Task:     Gene Set Enrichment Analysis (GSEA)
#' Author:   Imane BIYAR
#' Date:     July 2026
#'
#' Objective:
#' Map Entrez IDs and perform whole-transcriptome GSEA using the Wald statistic 
#' (from DESeq2) to identify system-level oncogenic reprogramming. 
#' Merge the resulting pathway flags (GO/KEGG) back onto the master_table 
#' for downstream Machine Learning integration.
#'
#' Inputs:
#'   - Results/Integration/master_table.rds
#'
#' Outputs:
#'   - Results/GSEA/gse_kegg_LUAD.rds & gse_go_LUAD.rds
#'   - Results/GSEA/*.png (Dotplots, Ridgeplots, Enrichment Maps, Cnetplots)
#'   - Results/Integration/master_table_pathways.rds (Updated backbone)
#' =============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(ggplot2)

dir.create("Results/GSEA", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 6.1 Load Master Data and ID Conversion
# ==============================================================================
master_table <- readRDS("Results/Integration/master_table.rds")

entrez_map <- bitr(master_table$gene, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

gsea_df <- left_join(master_table, entrez_map, by = c("gene" = "SYMBOL"))

gsea_df <- gsea_df %>%
  filter(!is.na(ENTREZID) & !is.na(stat))

gsea_df <- gsea_df %>%
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup()

# ==============================================================================
# 6.2 Generate the Ranked Gene List
# ==============================================================================
gsea_df <- gsea_df %>% arrange(desc(stat))

gene_list <- gsea_df$stat
names(gene_list) <- gsea_df$ENTREZID

# ==============================================================================
# 6.3 Execute GSEA (Gene Ontology & KEGG)
# ==============================================================================
set.seed(8)

if(!file.exists("Results/GSEA/gse_go_LUAD.rds")) {
  gse_go <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, ont = "BP",
                  minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05,
                  pAdjustMethod = "BH", eps = 0, seed = 8)
  
  gse_kegg <- gseKEGG(geneList = gene_list, organism = 'hsa',
                      minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05,
                      pAdjustMethod = "BH", eps = 0, seed = 8)
  
  saveRDS(gse_kegg, file = "Results/GSEA/gse_kegg_LUAD.rds")
  saveRDS(gse_go, file = "Results/GSEA/gse_go_LUAD.rds")
  
  write.csv(as.data.frame(gse_go), 'Results/GSEA/GSEA_GO_Results.csv')
  write.csv(as.data.frame(gse_kegg), 'Results/GSEA/GSEA_KEGG_Results.csv')
}

gse_go   <- readRDS("Results/GSEA/gse_go_LUAD.rds")
gse_kegg <- readRDS("Results/GSEA/gse_kegg_LUAD.rds")

# ==============================================================================
# 6.4 Visualizations
# ==============================================================================

png('Results/GSEA/GSEA_GO_Dotplot.png', width = 1200, height = 1000, res = 150)
enrichplot::dotplot(gse_go, showCategory = 10, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("Global Pathway Dysregulation (GO: Biological Process)")
dev.off()

kegg_plot <- enrichplot::dotplot(gse_kegg, showCategory = 10, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("Global Pathway Dysregulation (KEGG Signaling Pathways)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Results/GSEA/KEGG_Dysregulation_Dotplot.png", plot = kegg_plot,
       width = 14, height = 8, dpi = 300)

top_kegg_id <- gse_kegg$ID[1]
top_kegg_name <- gse_kegg$Description[1]
png('Results/GSEA/Top_KEGG_Enrichment.png', width = 1000, height = 800, res = 150)
gseaplot2(gse_kegg, geneSetID = top_kegg_id, title = top_kegg_name)
dev.off()

top_go_id <- gse_go$ID[1]
top_go_name <- gse_go$Description[1]
go_enrichment_plot <- gseaplot2(gse_go, geneSetID = top_go_id, title = paste("GO:", top_go_name))
ggsave("Results/GSEA/Top_GO_Enrichment_Plot.png", plot = go_enrichment_plot,
       width = 10, height = 8, dpi = 300)

png('Results/GSEA/GSEA_Ridgeplot_GO.png', width = 1200, height = 800, res = 150)
ridgeplot(gse_go, showCategory = 10) + ggtitle("Expression Distribution within Dysregulated GO Pathways")
dev.off()

png('Results/GSEA/GSEA_Ridgeplot_KEGG.png', width = 1200, height = 800, res = 150)
ridgeplot(gse_kegg, showCategory = 10) + ggtitle("Expression Distribution within Dysregulated KEGG Pathways")
dev.off()

gse_go_sim <- pairwise_termsim(gse_go, method = "JC")
png('Results/GSEA/GSEA_EnrichmentMap_GO.png', width = 1200, height = 1000, res = 150)
emapplot(gse_go_sim, showCategory = 20, color = "NES", layout = "kk") +
  ggtitle("Enrichment Map: Network of Dysregulated Pathways - GO")
dev.off()

gse_kegg_sim <- pairwise_termsim(gse_kegg, method = "JC")
png('Results/GSEA/GSEA_EnrichmentMap_KEGG.png', width = 1200, height = 1000, res = 150)
emapplot(gse_kegg_sim, showCategory = 20, color = "NES", layout = "kk") +
  ggtitle("Enrichment Map: Network of Dysregulated Pathways - KEGG")
dev.off()

gse_go_readable <- setReadable(gse_go, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
gene_list_symbols <- gsea_df$stat
names(gene_list_symbols) <- as.character(gsea_df$gene)
gene_list_symbols <- na.omit(gene_list_symbols)

cnet_plot <- cnetplot(gse_go_readable, showCategory = 5, foldChange = gene_list_symbols) +
  ggtitle("Gene-Concept Network: Core Drivers of Top Pathways - GO") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Results/GSEA/GO_Cnetplot_Readable.png", plot = cnet_plot, width = 14, height = 10, dpi = 300)

gse_kegg_readable <- setReadable(gse_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
cnet_plot_kegg <- cnetplot(gse_kegg_readable, showCategory = 5, foldChange = gene_list_symbols) +
  ggtitle("Gene-Concept Network: Core Drivers of Top Pathways - KEGG") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Results/GSEA/KEGG_Cnetplot_Readable.png", plot = cnet_plot_kegg, width = 14, height = 10, dpi = 300)

# ==============================================================================
# 6.5 Genes extraction: pathway_flag features
# ==============================================================================
kegg_results_df <- as.data.frame(gse_kegg)
core_kegg_drivers <- kegg_results_df$core_enrichment %>%
  strsplit(split = "/") %>% unlist() %>% unique()

go_results_df <- as.data.frame(gse_go)
core_go_drivers <- go_results_df$core_enrichment %>%
  strsplit(split = "/") %>% unlist() %>% unique()

gsea_df <- gsea_df %>%
  mutate(
    influences_kegg_pathway = ENTREZID %in% core_kegg_drivers,
    influences_go_pathway   = ENTREZID %in% core_go_drivers,
    is_any_pathway_driver   = influences_kegg_pathway | influences_go_pathway,
    is_both_pathway_driver  = influences_kegg_pathway & influences_go_pathway
  )

write.csv(gsea_df, "Results/GSEA/gsea_df.csv", row.names = FALSE)
saveRDS(gsea_df, file = "Results/GSEA/gsea_df.rds")

# ==============================================================================
# 6.6 Merge THREE explicit pathway columns back onto master_table
# ==============================================================================
pathway_flag_df <- gsea_df %>%
  transmute(gene = gene,
            in_GO_pathway          = as.numeric(influences_go_pathway),
            in_KEGG_pathway        = as.numeric(influences_kegg_pathway),
            in_GO_and_KEGG_pathway = as.numeric(is_both_pathway_driver)) %>%
  arrange(desc(in_GO_pathway), desc(in_KEGG_pathway), desc(in_GO_and_KEGG_pathway)) %>%
  distinct(gene, .keep_all = TRUE)

master_table <- readRDS("Results/Integration/master_table.rds") %>%
  dplyr::select(-any_of(c("pathway_flag", "in_GO_pathway", "in_KEGG_pathway",
                          "in_GO_and_KEGG_pathway"))) %>%   
  left_join(pathway_flag_df, by = "gene") %>%
  mutate(
    in_GO_pathway          = ifelse(is.na(in_GO_pathway), 0, in_GO_pathway),
    in_KEGG_pathway        = ifelse(is.na(in_KEGG_pathway), 0, in_KEGG_pathway),
    in_GO_and_KEGG_pathway = ifelse(is.na(in_GO_and_KEGG_pathway), 0, in_GO_and_KEGG_pathway)
  )

saveRDS(master_table, "Results/Integration/master_table_pathways.rds")

# ==============================================================================
# 6.7 Compact figures for the report appendix
# ==============================================================================
gse_kegg <- readRDS(file = "Results/GSEA/gse_kegg_LUAD.rds")
gse_go   <- readRDS(file = "Results/GSEA/gse_go_LUAD.rds")

## --- 6.7.1 Dotplots (top 5 pathways per direction) ---------------------------
go_dotplot_top5 <- enrichplot::dotplot(gse_go, showCategory = 5, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("Global Pathway Dysregulation (GO: Biological Process)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Results/GSEA/GSEA_GO_Dotplot_Top5.png", plot = go_dotplot_top5,
       width = 10, height = 6, dpi = 300)

kegg_dotplot_top5 <- enrichplot::dotplot(gse_kegg, showCategory = 5, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("Global Pathway Dysregulation (KEGG Signaling Pathways)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Results/GSEA/GSEA_KEGG_Dotplot_Top5.png", plot = kegg_dotplot_top5,
       width = 10, height = 6, dpi = 300)

## --- 6.7.2 Identify the single most dysregulated pathway per database --------
go_results_ranked <- as.data.frame(gse_go_readable) %>%
  arrange(desc(abs(NES)))
kegg_results_ranked <- as.data.frame(gse_kegg_readable) %>%
  arrange(desc(abs(NES)))

top_go_dysreg   <- go_results_ranked$Description[1]
top_kegg_dysreg <- kegg_results_ranked$Description[1]

## --- 6.7.3 Cnetplots restricted to that single pathway ------------------------
cnet_go_top1 <- cnetplot(gse_go_readable, showCategory = top_go_dysreg,
                         foldChange = gene_list_symbols) +
  ggtitle(paste0("Gene-Concept Network: ", top_go_dysreg, " (GO)")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Results/GSEA/GO_Cnetplot_TopPathway.png", plot = cnet_go_top1,
       width = 10, height = 8, dpi = 300)

cnet_kegg_top1 <- cnetplot(gse_kegg_readable, showCategory = top_kegg_dysreg,
                           foldChange = gene_list_symbols) +
  ggtitle(paste0("Gene-Concept Network: ", top_kegg_dysreg, " (KEGG)")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Results/GSEA/KEGG_Cnetplot_TopPathway.png", plot = cnet_kegg_top1,
       width = 10, height = 8, dpi = 300)