"""
Created by: Lin, Wei-Zhi
Last update: 20250331
version: 1.0.0
"""
#%% import packages
import pandas as pd 
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import os
#%% prepare function
def run_deseq2(
    df,
    metadata,
    list_control_major,
    list_case_major,
    list_control_minor,
    list_case_minor,
    list_design_factors,
    list_contrasts  # list of contrasts, e.g. [["B_KD", "1", "0"], ["PTH_Treatment", "1", "0"]]
    ):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    # Step 1: Prepare sample info
    sample_names = list_control_major + list_case_major + list_control_minor + list_case_minor
    metadata = metadata.loc[sample_names]

    # Step 2: Extract count data + gene info
    df = df[['gene_id'] + sample_names + ['gene_symbol']]
    gene_info = df[['gene_id', 'gene_symbol']]
    count_df = df[sample_names].round().astype(int).transpose()

    # Step 3: Run DESeq2
    dds = DeseqDataSet(
        counts=count_df,
        metadata=metadata,
        design_factors=list_design_factors,
        min_replicates=1
    )
    dds.deseq2()

    # Step 4: Run stats for each contrast
    all_results = {}
    for contrast in list_contrasts:
        stat_res = DeseqStats(dds, contrast=contrast)
        stat_res.summary()
        result_df = stat_res.results_df.copy()
        result_df["gene_id"] = gene_info["gene_id"].values
        result_df["gene_symbol"] = gene_info["gene_symbol"].values
        result_df = result_df.sort_values(by="pvalue")
        contrast_name = f"{contrast[0]}_{contrast[1]}_vs_{contrast[2]}"
        all_results[contrast_name] = result_df 

    return all_results, dds

def export_deseq2_results_to_csv(results_dict, output_dir="DESeq2_results"):
    os.makedirs(output_dir, exist_ok=True)

    for contrast_name, df in results_dict.items():
        filename = f"{contrast_name}.csv"
        filepath = os.path.join(output_dir, filename)
        df.to_csv(filepath, index=False)
        print(f"Exported: {filepath}")

def plot_volcano(
    result_df,
    title="Volcano Plot",
    padj_threshold=0.05,
    logfc_threshold=1,
    figsize=(8, 6),
    save_path=None
    ):
    import matplotlib.pyplot as plt
    import numpy as np 
    df = result_df.copy()
    df = df.dropna(subset=["padj", "log2FoldChange"])
    df["-log10(padj)"] = -np.log10(df["padj"])

    def get_color(row):
        if row["padj"] < padj_threshold and row["log2FoldChange"] > logfc_threshold:
            return "red"  
        elif row["padj"] < padj_threshold and row["log2FoldChange"] < -logfc_threshold:
            return "blue"  #
        else:
            return "grey"
    df["color"] = df.apply(get_color, axis=1)

    plt.figure(figsize=figsize)
    plt.scatter(df["log2FoldChange"], df["-log10(padj)"], c=df["color"], s=10, alpha=0.7)
    plt.axhline(-np.log10(padj_threshold), color='black', linestyle='--', linewidth=1)
    plt.axvline(logfc_threshold, color='black', linestyle='--', linewidth=1)
    plt.axvline(-logfc_threshold, color='black', linestyle='--', linewidth=1)
    plt.title(title)
    plt.xlabel("log2 FoldChange")
    plt.ylabel("-log10 padj")
    plt.grid(True)

    plt.show()

def compare_deseq2_results(df1, df2, top_n=100):
    """
    Compare two sets of DESeq2 results (must contain gene_id, log2FoldChange, and pvalue)

    Args:
        df1: First DESeq2 result DataFrame
        df2: Second DESeq2 result DataFrame
        top_n: Number of top genes with the largest differences to display

    Returns:
        merged_df: Merged comparison result DataFrame
    """
    import pandas as pd
    import matplotlib.pyplot as plt

    for df in [df1, df2]:
        assert "gene_id" in df.columns, "Missing gene_id"
        assert "log2FoldChange" in df.columns, "Missing log2FoldChange"
        assert "pvalue" in df.columns, "Missing pvalue"

    df_merged = pd.merge(
        df1[["gene_id", "log2FoldChange", "pvalue"]].rename(columns={
            "log2FoldChange": "log2FC_1", "pvalue": "pval_1"
        }),
        df2[["gene_id", "log2FoldChange", "pvalue"]].rename(columns={
            "log2FoldChange": "log2FC_2", "pvalue": "pval_2"
        }),
        on="gene_id", how="inner"
    )

    # add col difference
    df_merged["log2FC_diff"] = df_merged["log2FC_1"] - df_merged["log2FC_2"]
    df_merged["pval_diff"] = df_merged["pval_1"] - df_merged["pval_2"]

    # show the top n of log2FC
    print(f"\n Top {top_n} genes with largest log2FC difference:")
    print(df_merged.reindex(df_merged["log2FC_diff"].abs().sort_values(ascending=False).index).head(top_n)[
        ["gene_id", "log2FC_1", "log2FC_2", "log2FC_diff"]
    ])

    # calculating correlation
    corr = df_merged[["log2FC_1", "log2FC_2"]].corr().iloc[0,1]
    print(f"Pearson correlation of log2FC: {corr:.4f}")

    # plotting: difference of log2FC 
    plt.figure(figsize=(6,6))
    plt.scatter(df_merged["log2FC_1"], df_merged["log2FC_2"], s=10, alpha=0.5)
    plt.xlabel("log2FC (Version 1)")
    plt.ylabel("log2FC (Version 2)")
    plt.title("Log2 Fold Change Comparison")
    plt.axline((0, 0), slope=1, linestyle="--", color="gray")
    plt.text(
        0.05, 0.95,
        f"Pearson R = {corr:.4f}",
        transform=plt.gca().transAxes,
        ha="left", va="top",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray", alpha=0.7))
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    return df_merged

# %% 
"""
test area
"""
# import data 
df = pd.read_csv('all_compare.csv')
df = df.rename(columns={'gene_name': 'gene_symbol'})

#  metadata
list_sample = ['SC1_count', 'SC2_count', 'SC3_count', 'SC4_count', 
               'SC1P_count','SC2P_count', 'SC3P_count', 'SC4P_count', 
               'B1_count', 'B2_count', 'B3_count', 'B4_count', 
               'B1P_count', 'B2P_count', 'B3P_count','B4P_count'
              ]
df_metadata = pd.DataFrame({
    "Sample": list_sample,
    "B_KD": [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1],  # 0 = No KD, 1 = KD (Adjust accordingly)
    "PTH_Treatment": [0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1]  # 0 = No PTH, 1 = PTH (Adjust accordingly)
    })
df_metadata.set_index("Sample", inplace=True)
df_metadata["B_KD"] = df_metadata["B_KD"].astype(str).astype("category")
df_metadata["PTH_Treatment"] = df_metadata["PTH_Treatment"].astype(str).astype("category")

del list_sample

#%% test 2 factors
# run DEseq2
# test 2 factors
list_control_major = ['SC1P_count', 'SC2P_count', 'SC3P_count', 'SC4P_count']
list_case_major = ['B1P_count', 'B2P_count','B3P_count', 'B4P_count']
list_control_minor = ['SC1_count', 'SC2_count', 'SC3_count', 'SC4_count']
list_case_minor = ['B1_count', 'B2_count', 'B3_count', 'B4_count']
list_design_factors = ["PTH_Treatment", "B_KD"]
list_contrasts = [
    ["PTH_Treatment", "1", "0"], # PTH vs No PTH
    ["B_KD", "1", "0"]           # Knockdown vs WT
]

dict_results, dds = run_deseq2(
    df,
    df_metadata,
    list_control_major,
    list_case_major,
    list_control_minor,
    list_case_minor,
    list_design_factors,
    list_contrasts
)

# Virtualizaton
for contrast_name, result_df in dict_results.items():
    plot_volcano(result_df, 
                 title=f"Volcano: {contrast_name}",
                 padj_threshold=0.05,
                 logfc_threshold=1)
# export
export_deseq2_results_to_csv(dict_results)

#%% test 1 factors
# run DEseq2
# test 1 factors
list_control_major = ['SC1_count', 'SC2_count', 'SC3_count', 'SC4_count']
list_case_major = ['B1_count', 'B2_count', 'B3_count', 'B4_count']
list_control_minor = []
list_case_minor = []
list_design_factors = ["B_KD"]
list_contrasts = [
    ["B_KD", "1", "0"]           # Knockdown vs WT
]

dict_results, dds = run_deseq2(
    df,
    df_metadata,
    list_control_major,
    list_case_major,
    list_control_minor,
    list_case_minor,
    list_design_factors,
    list_contrasts
)

# Virtualizaton
for contrast_name, result_df in dict_results.items():
    plot_volcano(result_df, 
                 title=f"Volcano: {contrast_name}",
                 padj_threshold=0.05,
                 logfc_threshold=1)
# export
export_deseq2_results_to_csv(dict_results)
# %%
df_old = pd.read_csv('DE_B-PTH_vs_B-V_DEseq2_CDIAM_compare.csv')
df_new = pd.read_csv('DE_B-PTH_vs_S-PTH__pydeseq2_adj-B-KD_DEtool.csv')
merged_df = compare_deseq2_results(df_old, df_new, top_n=100)
# %%
