# Created by: Lin, Wei-Zhi
# Version: 1.1.0
# Last updated: 20250412

pathfind_2_KEGG <- function(input_csv, pathway_list, output_path) {
  library(gprofiler2)
  library(pathview)

  if (!file.exists(input_csv)) {
    stop("Input CSV file does not exist: ", input_csv)
  }

  # Step 1: Data Preprocessing
  pathway_list <- pathway_list[!is.na(pathway_list) & pathway_list != ""]

  input_mouse <- read.csv(input_csv)
  input_mouse <- na.omit(input_mouse)
  input_mouse <- subset(input_mouse, logFC <= -1.5 | logFC >= 1.5)
  input_mouse <- subset(input_mouse, adj.P.Val <= 0.05)
  input_mouse$adj.P.Val[input_mouse$adj.P.Val == 0] <- 1e-300

  # Step 2: Mouse ➜ Human gene mapping
  ortholog_result <- gorth(
    query = input_mouse$Gene.symbol,
    source_organism = "mmusculus",
    target_organism = "hsapiens",
    filter_na = TRUE
  )

  ortholog_result <- ortholog_result[, colSums(is.na(ortholog_result)) < nrow(ortholog_result)]
  colnames(ortholog_result)[1:2] <- c("Mouse.symbol", "Human.symbol")

  # Step 3: Merge expression + orthologs
  merged_df <- merge(
    input_mouse,
    ortholog_result,
    by.x = "Gene.symbol",
    by.y = "Human.symbol",
    suffixes = c("_mouse", "_human")
  )

  # Step 4: Prepare input for Pathview
  input_human <- merged_df[, c("ortholog_name", "logFC", "adj.P.Val")]
  colnames(input_human) <- c("Gene.symbol", "logFC", "adj.P.Val")

  n_duplicates <- sum(duplicated(input_human$Gene.symbol))
  input_human <- input_human[!duplicated(input_human$Gene.symbol), ]
  message("Removed ", n_duplicates, " duplicate gene(s) from input_human.")

  gene.vector <- input_human$logFC
  names(gene.vector) <- input_human$Gene.symbol

  # Step 5: Run Pathview
  original_wd <- getwd()
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  setwd(output_path)

  for (pid in pathway_list) {
    pathview(
      gene.data = gene.vector,
      pathway.id = pid,
      species = "hsa",
      gene.idtype = "SYMBOL",
      out.suffix = "pathview_output",
      kegg.native = TRUE
    )
  }

  # Step 6: Count PNG outputs
  png_count <- length(list.files(path = ".", pattern = "output.*\\.png$", ignore.case = TRUE))
  message("Number of output PNG files generated: ", png_count)

  setwd(original_wd)
}

# B-V vs S-V
list_pathway_id_BV_SV <- read.csv("pathfind_results_kegg_B-V_vs_S-V_20250411_label.csv")$ID
pathfind_2_KEGG(
  input_csv = "DE_B-V_vs_S-V_pydeseq2_DEtool_pathfind_input.csv",
  pathway_list = c(list_pathway_id_BV_SV),
  output_path = "KEGG_pathview_BV_vs_SV"
)

list_pathway_id_BV_SV <- list_pathway_id_BV_SV[!is.na(list_pathway_id_BV_SV) & list_pathway_id_BV_SV != ""]
message("Number of output PNG files generated: ", length(unique(list_pathway_id_BV_SV)))

# S-PTH vs S-V
list_pathway_id_SPTH_SV <- read.csv("pathfind_results_kegg_S-PTH_vs_S-V_20250411_label.csv")$ID
pathfind_2_KEGG(
  input_csv = "DE_S-PTH_vs_S-V_pydeseq2_DEtool_pathfind_input.csv",
  pathway_list = c(list_pathway_id_SPTH_SV),
  output_path = "KEGG_pathview_SPTH_vs_SV"
)

list_pathway_id_SPTH_SV <- list_pathway_id_SPTH_SV[!is.na(list_pathway_id_SPTH_SV) & list_pathway_id_SPTH_SV != ""]
message("Number of output PNG files generated: ", length(unique(list_pathway_id_SPTH_SV)))

# B-PTH vs B-V
list_pathway_id_BPTH_BV <- read.csv("pathfind_results_kegg_B-PTH_vs_B-V_20250411_label.csv")$ID
pathfind_2_KEGG(
  input_csv = "DE_B-PTH_vs_B-V_pydeseq2_DEtool_pathfind_input.csv",
  pathway_list = c(list_pathway_id_BPTH_BV),
  output_path = "KEGG_pathview_BPTH_vs_BV"
)

list_pathway_id_BPTH_BV <- list_pathway_id_BPTH_BV[!is.na(list_pathway_id_BPTH_BV) & list_pathway_id_BPTH_BV != ""]
message("Number of output PNG files generated: ", length(unique(list_pathway_id_BPTH_BV)))

# B-PTH vs S-PTH adj B-KD
list_pathway_id_BPTH_SPTH_adjBKD <- read.csv("pathfind_results_kegg_B-PTH_vs_S-PTH_adj-B-KD_20250403_label.csv")$ID
pathfind_2_KEGG(
  input_csv = "DE_B-PTH_vs_S-PTH_pydeseq2_adj-B-KD_DEtool_pathfind_input.csv",
  pathway_list = c(list_pathway_id_BPTH_SPTH_adjBKD),
  output_path = "KEGG_pathview_BPTH_vs_SPTH_adjBKD"
)

list_pathway_id_BPTH_SPTH_adjBKD <- list_pathway_id_BPTH_SPTH_adjBKD[!is.na(list_pathway_id_BPTH_SPTH_adjBKD) & list_pathway_id_BPTH_SPTH_adjBKD != ""]
message("Number of output PNG files generated: ", length(unique(list_pathway_id_BPTH_SPTH_adjBKD)))

# reset path, if necessary
setwd('H:/我的雲端硬碟/Harvard_projects/B6mince_BambiKD_PTH/')
