"""
Created by Lin, Wei-Zhi
Version: v1.0.0  
Last updated: 20250403
Note: input from Pathfinder.R; pyvis.network
"""
#%% 00 Function
# 00 Function
def create_pathfinder_network(df_de, df_path, 
                              output_html="Pathfinder_Network.html", 
                              highlight_terms=None):
    """
    Generate an interactive gene-pathway network using Pathfinder.R output and DE results.

    Parameters:
    - df_de: DataFrame of DE results (must contain 'gene_symbol', 'log2FoldChange', 'padj')
    - df_path: DataFrame of Pathfinder results (must contain 'Term_Description', 'lowest_p',
               'Up_regulated', and 'Down_regulated')
    - output_html: filename of the output .html file
    - highlight_terms: optional list of pathway terms to highlight
    """
    import pandas as pd
    import networkx as nx
    from pyvis.network import Network
    import webbrowser
    import os

    # Filter significant pathways
    df_path_sig = df_path[df_path["lowest_p"] < 0.05].sort_values("lowest_p")

    # Filter significant DE genes
    df_sig = df_de[df_de["padj"] < 0.05]
    dict_logfc = dict(zip(df_sig["gene_symbol"], df_sig["log2FoldChange"]))

    highlight_terms = set(highlight_terms or [])
    G = nx.Graph()

    for _, row in df_path_sig.iterrows():
        term = row["Term_Description"]
        up_genes = str(row.get("Up_regulated", "")).split(",") if pd.notna(row.get("Up_regulated")) else []
        down_genes = str(row.get("Down_regulated", "")).split(",") if pd.notna(row.get("Down_regulated")) else []
        genes = up_genes + down_genes

        if not genes:
            continue

        # Set pathway as a circle with a black border and larger size
        node_color = "red" if term in highlight_terms else "lightgray"
        G.add_node(term, type="pathway", 
                   color={"border": "black", "background": node_color},
                   shape="dot")

        for gene in genes:
            gene = gene.strip()
            fc = dict_logfc.get(gene, 0)
            color = "green" if fc > 0 else "red" if fc < 0 else "gray"
            G.add_node(gene, type="gene", color=color, shape="dot")
            G.add_edge(term, gene)

    # Visualization
    net = Network(height="800px", width="100%", notebook=False, directed=False)
    net.force_atlas_2based()

    for node, data in G.nodes(data=True):
        size = 15 + abs(dict_logfc.get(node, 0)) * 5 if data["type"] == "gene" else 30
        title = f"{node}<br>Type: {data['type']}"
        shape = data.get("shape", "dot")
        color = data.get("color", "gray")
        borderWidth = data.get("borderWidth", 1)

        net.add_node(node, label=node, color=color, size=size, title=title,
                     shape=shape, borderWidth=borderWidth)

    for source, target in G.edges():
        net.add_edge(source, target)

    net.show(output_html)
    print(f"Network graph saved as: {output_html}")
    abs_path = os.path.abspath(output_html)
    webbrowser.open(f'file://{abs_path}')

#%% 01 test run
# 01 test run
# File paths
import pandas as pd
df_de = pd.read_csv('DE_B-PTH_vs_B-V_pydeseq2_DEtool_pathfinder_input.csv')
df_pathfinder = pd.read_csv("pathfindr_results_kegg_B-PTH_vs_S-PTH_adj-B-KD_20250403.csv")

# Rename columns to match expected format
df_de = df_de.rename(columns={
    "Gene.symbol": "gene_symbol",
    "logFC": "log2FoldChange",
    "adj.P.Val": "padj"
})

# Optional: highlight specific pathways
highlight_terms = ['Parathyroid hormone synthesis, secretion and action',
                   'Th17 cell differentiation',
                   'JAK-STAT signaling pathway',
                   'MAPK signaling pathway',
                   'Efferocytosis',
                   'Th1 and Th2 cell differentiation',
                   'AGE-RAGE signaling pathway in diabetic complications',
                   'EGFR tyrosine kinase inhibitor resistance',
                   'Alanine, aspartate and glutamate metabolism',
                   'Rheumatoid arthritis',
                   'Complement and coagulation cascades',
                   'Oxytocin signaling pathway',
                   'Fluid shear stress and atherosclerosis',
                   'Toll-like receptor signaling pathway',
                   'HIF-1 signaling pathway',
                   'Prolactin signaling pathway',
                   'Natural killer cell mediated cytotoxicity',
                   'AMPK signaling pathway',
                   'Wnt signaling pathway'
                   ]

# Call the function to generate interactive HTML network
create_pathfinder_network(df_de, df_pathfinder, 
                          output_html="pathfindr_network_kegg_B-PTH_vs_S-PTH_adj-B-KD.html", 
                          highlight_terms=highlight_terms)
#%%
