
# Created by: Lin, Wei-Zhi
# Version: 1.1.0
# Last updated: 20250411


pathfind_kegg <- function(input_csv, output_csv, logfc_thresh = 1.5, padj_thresh = 0.05) {
  # Load required libraries (suppress messages for cleaner output)
  suppressPackageStartupMessages({
    library(biomaRt)
    library(gprofiler2)
    library(pathview)
    library(pathfindR)
    library(org.Hs.eg.db)
  })

  # Step 1: Load and filter differential expression data
  input_mouse <- read.csv(input_csv)
  message("Initial number of genes: ", nrow(input_mouse))

  input_mouse <- na.omit(input_mouse)
  input_mouse <- subset(input_mouse, logFC <= -logfc_thresh | logFC >= logfc_thresh)
  input_mouse <- subset(input_mouse, adj.P.Val <= padj_thresh)
  input_mouse$adj.P.Val[input_mouse$adj.P.Val == 0] <- 1e-300
  message("After logFC ≥ ", logfc_thresh, " and adj.P.Val ≤ ", padj_thresh, ": ", nrow(input_mouse), " genes")

  if (nrow(input_mouse) == 0) {
    warning("No genes passed filtering. Try relaxing thresholds.")
    return(NULL)
  }

  # Step 2: Convert mouse genes to human orthologs using gorth()
  mouse_genes <- input_mouse$Gene.symbol
  ortholog_result <- gorth(
    query = mouse_genes,
    source_organism = "mmusculus",
    target_organism = "hsapiens",
    filter_na = TRUE
  )

  # Keep only necessary columns and rename for clarity
  ortholog_result <- ortholog_result[, colSums(is.na(ortholog_result)) < nrow(ortholog_result)]
  colnames(ortholog_result)[colnames(ortholog_result) == "input"] <- "Mouse.symbol"
  colnames(ortholog_result)[colnames(ortholog_result) == "ortholog_name"] <- "Human.symbol"
  message("Orthologs mapped: ", nrow(ortholog_result))

  # Step 3: Merge expression data with ortholog mapping
  input_mouse$Gene.symbol <- as.character(input_mouse$Gene.symbol)
  ortholog_result$Mouse.symbol <- as.character(ortholog_result$Mouse.symbol)

  merged_df <- merge(
    input_mouse,
    ortholog_result,
    by.x = "Gene.symbol",
    by.y = "Mouse.symbol",
    suffixes = c("_mouse", "_human")
  )

  input_human <- merged_df[, c("Human.symbol", "logFC", "adj.P.Val")]
  colnames(input_human) <- c("Gene.symbol", "logFC", "adj.P.Val")

  # Remove duplicates based on human gene symbol
  n_duplicates <- sum(duplicated(input_human$Gene.symbol))
  input_human <- input_human[!duplicated(input_human$Gene.symbol), ]
  message("Removed ", n_duplicates, " duplicate gene(s) from input_human.")
  message("Final gene count for enrichment: ", nrow(input_human))

  # Step 4: Check if enough genes remain
  if (nrow(input_human) < 2) {
    warning("Not enough genes for enrichment analysis. Only ", nrow(input_human), " genes found after mapping.")
    return(NULL)
  }

  # Step 5: Perform KEGG enrichment
  results_kegg <- run_pathfindR(
    input_human,
    gene_sets = "KEGG",
    p_val_threshold = 0.05
  )

  # Step 6: Save enrichment results to CSV
  write.csv(results_kegg, output_csv, row.names = FALSE)

  message("Analysis complete. Results saved to: ", output_csv)
}

# run on each DE 
pathfind_kegg("DE_B-V_vs_S-V_pydeseq2_DEtool_pathfind_input.csv",
              "pathfind_results_kegg_B-V_vs_S-V_20250411.csv",
              logfc_thresh = 1.5,
              padj_thresh = 0.05)

pathfind_kegg("DE_S-PTH_vs_S-V_pydeseq2_DEtool_pathfind_input.csv", 
              "pathfind_results_kegg_S-PTH_vs_S-V_20250411")

pathfind_kegg("DE_B-PTH_vs_B-V_pydeseq2_DEtool_pathfind_input.csv",
              "pathfind_results_kegg_B-PTH_vs_B-V_20250411")

pathfind_kegg("DE_B-PTH_vs_S-PTH_pydeseq2_adj-B-KD_DEtool_pathfind_input.csv",
              "pathfind_results_kegg_B-PTH_vs_S-PTH_adj-B-KD_20250403")
