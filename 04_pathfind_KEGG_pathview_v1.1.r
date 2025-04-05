
# Created by: Lin, Wei-Zhi
# Version: 1.1.0
# Last updated: 20250405
# Note: need to be further packaged into a tool.

# Setup: install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("KEGGgraph", "Rgraphviz", "graph"))
BiocManager::install("pathview", force = TRUE)
BiocManager::install(c("org.Hs.eg.db","org.Mm.eg.db"))
BiocManager::install("biomaRt", force = TRUE)
install.packages("gprofiler2")
install.packages("pathfindR")

# Step 0: Load Necessary Libraries
library(biomaRt)
library(gprofiler2)
library(pathview)
library(pathfindR)
library(org.Hs.eg.db)

# Step 1: Read and Filter Mouse Differential Expression Data
input_mouse <- read.csv("DE_B-PTH_vs_S-PTH__pydeseq2_adj-B-KD_DEtool_pathfindr_input.csv")  # Should include Gene.symbol, logFC, adj.P.Val
input_mouse <- na.omit(input_mouse)
input_mouse <- subset(input_mouse, logFC <= -1.5 | logFC >= 1.5) # filter the row that logFC <=-1.5 or >= 1.5 
input_mouse <- subset(input_mouse, adj.P.Val <= 0.05) # filter the row that adj p val < 0.05 
# Replace zero p-values with 1e-300
input_mouse$adj.P.Val[input_mouse$adj.P.Val == 0] <- 1e-300
# Check input df
head(input_mouse)

# Step 2: Retrieve Human Orthologs for Mouse Genes
mouse_genes <- input_mouse$Gene.symbol
# mouse ➜ human ortholog
ortholog_result <- gorth(
  query = mouse_genes,
  source_organism = "mmusculus",
  target_organism = "hsapiens",
  filter_na = TRUE
)

# Remove columns with all NA values from ortholog_result
ortholog_result <- ortholog_result[, colSums(is.na(ortholog_result)) < nrow(ortholog_result)]
colnames(ortholog_result)[1:2] <- c("Mouse.symbol", "Human.symbol")

# Step 3: Merge Orthologs with Expression Data
merged_df <- merge(
  input_mouse,
  ortholog_result,
  by.x = "Gene.symbol",
  by.y = "Human.symbol",
  suffixes = c("_mouse", "_human")
)

# Step 4: Prepare Data for Pathway Analysis
input_human <- merged_df[, c("ortholog_name", "logFC", "adj.P.Val")]
colnames(input_human) <- c("Gene.symbol", "logFC", "adj.P.Val")

# Step 5: Save the Processed Data
write.csv(input_human, "pathfind_input_human_20250402.csv", row.names = FALSE)

# Step 6: Perform Pathway Enrichment Analysis (KEGG)
results_kegg <- run_pathfindR(
  input_human,
  gene_sets = "KEGG",
  p_val_threshold = 0.05
)

# Create the term-gene graph
p <- term_gene_graph(
  result_df = results_kegg,
  num_terms = 15,
  use_description = TRUE,
  node_size = "p_val"
)

# Step 7: Save the Enrichment Results
write.csv(results_kegg, "pathfind_results_kegg_human_ortholog_20250403.csv", row.names = FALSE)

# Step8: Viturlization: KEGG

library(pathview)

input_human <- read.csv("pathfind_input_human_20250402.csv")

gene.vector <- input_human$logFC
names(gene.vector) <- input_human$Gene.symbol
pathways <- c('hsa04928','hsa04659','hsa04630','hsa04010','hsa04148','hsa04658','hsa04933',
              'hsa01521','hsa00250','hsa05323','hsa04610','hsa04921','hsa05418','hsa04620',
              'hsa04066','hsa04917','hsa04650','hsa04152','hsa04310'
              ) 

original_wd <- getwd()
# if (!dir.exists("KEGG_pathview")) {dir.create("KEGG_pathview")}
# setwd("KEGG_pathview")
setwd("KEGG_pathview_B_PTH_vs_S_PTH 20250405")

for (pid in pathways) {
  pathview(gene.data = gene.vector,
         pathway.id = pid,
         species = "hsa",
         gene.idtype = "SYMBOL",
         out.suffix = "B_PTH_vs_S_PTH",
         kegg.native = TRUE)
}

setwd(original_wd)
