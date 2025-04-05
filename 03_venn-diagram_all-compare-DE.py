#%% 
# import packages
import pandas as pd
from matplotlib import pyplot as plt
from venn import venn

#%% 
# load the dataset
df_BV_SV = pd.read_csv('DE_B-V_vs_S-V_pydeseq2.csv')
df_SPTH_SV = pd.read_csv('DE_S-PTH_vs_S-V_pydeseq2.csv')
df_BPTH_BV = pd.read_csv('DE_B-PTH_vs_B-V_pydeseq2.csv')
df_BPTH_SPTH = pd.read_csv('DE_B-PTH_vs_S-PTH__pydeseq2_adj-B-KD.csv')

#%% 
# pre process df
adj_threshold = 0.05
log2FC_threshold = 1.5

# Extract significant gene sets from each DataFrame
genes_BV_SV = set(df_BV_SV[(df_BV_SV['padj'] < adj_threshold) & 
                            ((df_BV_SV['log2FoldChange'] > log2FC_threshold) | 
                             (df_BV_SV['log2FoldChange'] < -log2FC_threshold))]
                            ['gene_id'])

genes_SPTH_SV = set(df_SPTH_SV[(df_SPTH_SV['padj'] < adj_threshold) & 
                            ((df_SPTH_SV['log2FoldChange'] > log2FC_threshold) | 
                             (df_SPTH_SV['log2FoldChange'] < -log2FC_threshold))]
                            ['gene_id'])

genes_BPTH_BV = set(df_BPTH_BV[(df_BPTH_BV['padj'] < adj_threshold) & 
                            ((df_BPTH_BV['log2FoldChange'] > log2FC_threshold) | 
                             (df_BPTH_BV['log2FoldChange'] < -log2FC_threshold))]
                            ['gene_id'])

genes_BPTH_SPTH = set(df_BPTH_SPTH[(df_BPTH_SPTH['padj'] < adj_threshold) & 
                            ((df_BPTH_SPTH['log2FoldChange'] > log2FC_threshold) | 
                             (df_BPTH_SPTH['log2FoldChange'] < -log2FC_threshold))]
                            ['gene_id'])

print ('# of signifiacnt gene in B-V_vs_S-V is ', len(genes_BV_SV))
print ('# of signifiacnt gene in S-PTH_vs_S-V is ', len(genes_SPTH_SV))
print ('# of signifiacnt gene in B-PTH_vs_B-V is ', len(genes_BPTH_BV))
print ('# of signifiacnt gene in B-PTH_vs_S-PTH_adj-B-KD is ', len(genes_BPTH_SPTH))
#%% 
# Venn
gene_sets = {
    f"B-V_vs_S-V (n={len(genes_BV_SV)})": genes_BV_SV,
    f"S-PTH_vs_S-V (n={len(genes_SPTH_SV)})": genes_SPTH_SV,
    f"B-PTH_vs_B-V (n={len(genes_BPTH_BV)})": genes_BPTH_BV,
    f"B-PTH_vs_S-PTH adj B-KD (n={len(genes_BPTH_SPTH)})": genes_BPTH_SPTH,
}

# Plot the Venn diagram
plt.figure(figsize=(8,8))
venn(gene_sets,  
     cmap="coolwarm", 
     fontsize=10)
plt.title("Venn Diagram of Significant Genes")
plt.show()

#%% 
# extract all the overlapping gene
# Get all unique & overlapping gene sets
only_BV_SV = genes_BV_SV - genes_SPTH_SV - genes_BPTH_BV - genes_BPTH_SPTH
only_SPTH_SV = genes_SPTH_SV - genes_BV_SV - genes_BPTH_BV - genes_BPTH_SPTH
only_BPTH_BV = genes_BPTH_BV - genes_BV_SV - genes_SPTH_SV - genes_BPTH_SPTH
only_BPTH_SPTH = genes_BPTH_SPTH - genes_BV_SV - genes_SPTH_SV - genes_BPTH_BV

# Two-way overlaps
BV_SV_SPTH_SV = genes_BV_SV & genes_SPTH_SV
BV_SV_BPTH_BV = genes_BV_SV & genes_BPTH_BV
BV_SV_BPTH_SPTH = genes_BV_SV & genes_BPTH_SPTH
SPTH_SV_BPTH_BV = genes_SPTH_SV & genes_BPTH_BV
SPTH_SV_BPTH_SPTH = genes_SPTH_SV & genes_BPTH_SPTH
BPTH_BV_BPTH_SPTH = genes_BPTH_BV & genes_BPTH_SPTH

# Three-way overlaps
BV_SV_SPTH_SV_BPTH_BV = genes_BV_SV & genes_SPTH_SV & genes_BPTH_BV
BV_SV_SPTH_SV_BPTH_SPTH = genes_BV_SV & genes_SPTH_SV & genes_BPTH_SPTH
BV_SV_BPTH_BV_BPTH_SPTH = genes_BV_SV & genes_BPTH_BV & genes_BPTH_SPTH
SPTH_SV_BPTH_BV_BPTH_SPTH = genes_SPTH_SV & genes_BPTH_BV & genes_BPTH_SPTH

# All four-way overlap
all_overlap = genes_BV_SV & genes_SPTH_SV & genes_BPTH_BV & genes_BPTH_SPTH

overlap_dict = {
    "Only B-V_vs_S-V": only_BV_SV,
    "Only S-PTH_vs_S-V": only_SPTH_SV,
    "Only B-PTH_vs_B-V": only_BPTH_BV,
    "Only B-PTH_vs_S-PTH adj B-KD": only_BPTH_SPTH,
    "B-V_vs_S-V & S-PTH_vs_S-V": BV_SV_SPTH_SV,
    "B-V_vs_S-V & B-PTH_vs_B-V": BV_SV_BPTH_BV,
    "B-V_vs_S-V & B-PTH_vs_S-PTH adj B-KD": BV_SV_BPTH_SPTH,
    "S-PTH_vs_S-V & B-PTH_vs_B-V": SPTH_SV_BPTH_BV,
    "S-PTH_vs_S-V & B-PTH_vs_S-PTH adj B-KD": SPTH_SV_BPTH_SPTH,
    "B-PTH_vs_B-V & B-PTH_vs_S-PTH adj B-KD": BPTH_BV_BPTH_SPTH,
    "B-V_vs_S-V & S-PTH_vs_S-V & B-PTH_vs_B-V": BV_SV_SPTH_SV_BPTH_BV,
    "B-V_vs_S-V & S-PTH_vs_S-V & B-PTH_vs_S-PTH adj B-KD": BV_SV_SPTH_SV_BPTH_SPTH,
    "B-V_vs_S-V & B-PTH_vs_B-V & B-PTH_vs_S-PTH adj B-KD": BV_SV_BPTH_BV_BPTH_SPTH,
    "S-PTH_vs_S-V & B-PTH_vs_B-V & B-PTH_vs_S-PTH adj B-KD": SPTH_SV_BPTH_BV_BPTH_SPTH,
    "All Four Overlap": all_overlap
}

#%%
# Convert overlap dictionary into a list of (Region, Gene) pairs
overlap_list = []
for region, genes in overlap_dict.items():
    for gene in genes:
        overlap_list.append({"Region": region, "Gene": gene})

# Convert to DataFrame
df_overlap = pd.DataFrame(overlap_list)

# Save to CSV
df_overlap.to_csv("venn-diagram_all-compare-DE-list.csv", index=False)
