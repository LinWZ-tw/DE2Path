"""
Created by Lin, Wei-Zhi
Version: v2.2.0  
Last updated: 20250403
Note: gsea; pyvis.network
"""
#%% 00 packages 
# 00 import packages
import gseapy as gp
import pandas as pd
# pd.options.display.float_format = '{:.10f}'.format
import time
from datetime import datetime, timedelta

#%% 01 functions 
# 01 prepare functions
# def preprocess_df(df):
#     # Ensure required columns exist
#     required_cols = {'gene_symbol', 'log2FoldChange'}
#     if not required_cols.issubset(df.columns):
#         raise ValueError(f"Input DataFrame must contain columns: {required_cols}")
#     df = df[['gene_symbol', 'log2FoldChange']].dropna()
#     df = df.sort_values(by='log2FoldChange', ascending=False)
#     return df

def run_gsea(rnk_file, gene_sets, output_csv, title_prefix=""):
    """
    Run GSEA prerank analysis and save results to CSV.

    Parameters:
    - rnk_file: path to the .rnk input file
    - gene_sets: gene set libraries (list or str)
    - output_csv: path to save GSEA result table (.csv)
    - title_prefix: optional prefix to prepend to plot titles
    """

    print(f"Starting GSEA for: {rnk_file}")
    start_time = time.time()
    print("GSEA start at:", datetime.fromtimestamp(start_time).strftime('%Y-%m-%d %H:%M:%S'))

    # Run GSEA prerank analysis
    gsea_result = gp.prerank(
        rnk=rnk_file,
        gene_sets=gene_sets,
        outdir=None  # prevent auto file generation
    )

    # Prepare and clean result DataFrame
    df_gsea = gsea_result.res2d.copy()
    df_gsea["FDR q-val"] = pd.to_numeric(df_gsea["FDR q-val"], errors="coerce").fillna(1)

    # Save to CSV only
    df_gsea.to_csv(output_csv, index=False)
    print(f"GSEA results saved to CSV: {output_csv}")

    # Timing
    end_time = time.time()
    print("End at:", datetime.fromtimestamp(end_time).strftime('%Y-%m-%d %H:%M:%S'))
    elapsed_time = str(timedelta(seconds=int(end_time - start_time)))
    print(f"Execution Time: {elapsed_time}")

    return gsea_result

#%% 02 loading data
# 02 loading data
df_BV_SV = pd.read_csv('DE_B-V_vs_S-V_pydeseq2_DEtool.csv')
df_SPTH_SV = pd.read_csv('DE_S-PTH_vs_S-V_pydeseq2_DEtool.csv')
df_BPTH_BV = pd.read_csv('DE_B-PTH_vs_B-V_pydeseq2_DEtool.csv')
df_BPTH_SPTH = pd.read_csv('DE_B-PTH_vs_S-PTH_pydeseq2_adj-B-KD_DEtool.csv')
df_BPTH_SPTH ["gene_symbol"] = df_BPTH_SPTH ["gene_symbol"].str.replace(r"[^a-zA-Z0-9_-]", "", regex=True)

# select gene set
''' confirm the species'''
gene_sets = [
    "KEGG_2019_Mouse",
    "GO_Biological_Process_2025",
    "Reactome_Pathways_2024"
    ]

#%% 03 pre process data
# 03 pre process data
# drop not significant data
adj_threshold = 0.05
log2FC_threshold = 1.5

"""filtered based on both adj p-val and log2FC"""
df_BV_SV = df_BV_SV[(df_BV_SV["padj"] < adj_threshold) & (abs(df_BV_SV["log2FoldChange"]) > log2FC_threshold)]

df_SPTH_SV = df_SPTH_SV[(df_SPTH_SV["padj"] < adj_threshold) & (abs(df_SPTH_SV["log2FoldChange"]) > log2FC_threshold)]

df_BPTH_BV = df_BPTH_BV[(df_BPTH_BV["padj"] < adj_threshold) & (abs(df_BPTH_BV["log2FoldChange"]) > log2FC_threshold)]

df_BPTH_SPTH = df_BPTH_SPTH[(df_BPTH_SPTH["padj"] < adj_threshold) & (abs(df_BPTH_SPTH["log2FoldChange"]) > log2FC_threshold)]

"""filtered based on adj p-val"""
# df_BV_SV = df_BV_SV[(df_BV_SV["padj"] < adj_threshold) ]

# df_SPTH_SV = df_SPTH_SV[(df_SPTH_SV["padj"] < adj_threshold) ]

# df_BPTH_BV = df_BPTH_BV[(df_BPTH_BV["padj"] < adj_threshold) ]

# df_BPTH_SPTH = df_BPTH_SPTH[(df_BPTH_SPTH["padj"] < adj_threshold) ]

"""drop not necessary columns"""
df_BV_SV[["gene_symbol", "log2FoldChange"]].dropna().to_csv("gsea_input_BV_vs_SV_pydeseq2.rnk", 
                                                            sep="\t", index=False, header=False, 
                                                            encoding='utf-8')

df_SPTH_SV[["gene_symbol", "log2FoldChange"]].dropna().to_csv("gsea_input_SPTH_vs_SV_pydeseq2.rnk", 
                                                              sep="\t", index=False, header=False, 
                                                              encoding='utf-8')

df_BPTH_BV[["gene_symbol", "log2FoldChange"]].dropna().to_csv("gsea_input_BPTH_vs_BV_pydeseq2.rnk", 
                                                              sep="\t", index=False, header=False, 
                                                              encoding='utf-8')

df_BPTH_SPTH[["gene_symbol", "log2FoldChange"]].dropna().to_csv("gsea_input_BPTH_vs_SPTH_pydeseq2.rnk", 
                                                                sep="\t", index=False, header=False, 
                                                                encoding='utf-8')

del adj_threshold, log2FC_threshold

#%% 04 run gsea
# 04 run gsea
# select gene set
''' select gene set; confirm the species'''
gene_sets = [
    "KEGG_2019_Mouse",
    "GO_Biological_Process_2025",
    "Reactome_Pathways_2024",
]

gsea_BV_SV = run_gsea(
    rnk_file='gsea_input_BV_vs_SV_pydeseq2.rnk',
    gene_sets=gene_sets,
    output_csv='GSEA_output_B-V_vs_S-V_pydeseq2.csv',
    title_prefix="B-V vs S-V "
)

gsea_SPTH_SV = run_gsea(
    rnk_file='gsea_input_SPTH_vs_SV_pydeseq2.rnk',
    gene_sets=gene_sets,
    output_csv='GSEA_output_S-PTH_vs_S-V_pydeseq2.csv',
    title_prefix="S-PTH vs S-V "
)

gsea_BPTH_BV = run_gsea(
    rnk_file='gsea_input_BPTH_vs_BV_pydeseq2.rnk',
    gene_sets=gene_sets,
    output_csv='GSEA_output_B-PTH_vs_B-V_pydeseq2.csv',
    title_prefix="B-PTH vs B-V "
)

gsea_BPTH_SPTH = run_gsea(
    rnk_file='gsea_input_B-PTH_vs_S-PTH_pydeseq2_adj-B-KD.rnk',
    gene_sets=gene_sets,
    output_csv='GSEA_output_B-PTH_vs_S-PTH_pydeseq2_adj-B-KD.csv',
    title_prefix="B-PTH vs S-PTH adj B-KD"
)

#%% 05 Plotting network
# 05 plot the network
def create_gsea_network(de_file, gsea_file, 
                        output_html="GSEA_PathwayGene_Network.html", 
                        highlight_terms=None):
    import pandas as pd
    import networkx as nx
    from pyvis.network import Network
    import webbrowser
    import os
    # Read the DE and GSEA data from CSV files
    df_de = pd.read_csv(de_file)
    df_gsea = pd.read_csv(gsea_file)
    # Filter significant pathways with FDR q-value < 0.05
    df_gsea_sig = df_gsea[df_gsea["FDR q-val"] < 0.05].sort_values("FDR q-val")
    # Filter significant genes with adjusted p-value < 0.05
    df_sig = df_de[df_de["padj"] < 0.05]
    dict_logfc = dict(zip(df_sig["gene_symbol"], df_sig["log2FoldChange"]))
    # Highlight pathway of interest
    highlight_terms = set(highlight_terms or [])

    # Initialize an undirected graph
    G = nx.Graph()

    # Loop through significant enriched terms and connect to their genes
    for _, row in df_gsea_sig.iterrows():
        term = row["Term"]
        if pd.isna(row.get("Lead_genes")):
            continue
        genes = row["Lead_genes"].split(";")

        highlight = term in highlight_terms
        node_color = "red" if highlight else "lightgray"
        G.add_node(term, type="pathway", color=node_color)

        for gene in genes:
            gene = gene.strip()
            fc = dict_logfc.get(gene, 0)
            color = "green" if fc > 0 else "red" if fc < 0 else "gray"
            G.add_node(gene, type="gene", color=color)
            G.add_edge(term, gene)

    # Initialize PyVis network
    net = Network(height="800px", width="100%", notebook=False, directed=False)
    net.force_atlas_2based()

    # Add nodes and edges from NetworkX graph to PyVis network
    for node, data in G.nodes(data=True):
        color = data['color']
        size = 15 + abs(dict_logfc.get(node, 0)) * 5
        title = f"{node}<br>Type: {data['type']}"

        net.add_node(node, label=node, color=color, size=size, title=title)

    for source, target in G.edges():
        net.add_edge(source, target)
    
    # Generate and save the network graph
    net.show(output_html)
    print(f"Network graph saved as: {output_html}")

    abs_path = os.path.abspath(output_html)
    webbrowser.open(f'file://{abs_path}')

# test run:
# B-V_vs_S-V
create_gsea_network('DE_B-V_vs_S-V_pydeseq2_DEtool.csv', 
                    'GSEA_output_B-V_vs_S-V_pydeseq2.csv',
                    output_html='GSEA_network_B-V_vs_S-V_pydeseq2.html')

# S-PTH_vs_S-V
create_gsea_network('DE_S-PTH_vs_S-V_pydeseq2_DEtool.csv', 
                    'GSEA_output_S-PTH_vs_S-V_pydeseq2.csv',
                    output_html='GSEA_network_S-PTH_vs_S-V_pydeseq2.html')

# B-PTH_vs_B-V
create_gsea_network('DE_B-PTH_vs_B-V_pydeseq2_DEtool.csv', 
                    'GSEA_output_B-PTH_vs_B-V_pydeseq2.csv',
                    output_html='GSEA_network_B-PTH_vs_B-V_pydeseq2.html')

# B-PTH_vs_S-PTH_adj-B-KD
create_gsea_network('DE_B-PTH_vs_S-PTH_pydeseq2_adj-B-KD_DEtool.csv', 
                    'GSEA_output_B-PTH_vs_S-PTH_pydeseq2_adj-B-KD.csv',
                    output_html='GSEA_network_B-PTH_vs_S-PTH_pydeseq2_adj-B-KD.html',
                    highlight_terms=['GO_Biological_Process_2025__Response to Estrogen (GO:0043627)']
                    )

#%%
