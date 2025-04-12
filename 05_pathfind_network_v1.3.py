"""
Created by Lin, Wei-Zhi
Version: v1.3.0
Last updated: 20250410
Note: input from Pathfinder.R; pyvis.network
"""
#%% 00 Function
# 00 Function
def create_pathfinder_network(df_de, df_path,
                              output_html="Pathway_Network.html",
                              highlight_terms=None,
                              connect_pathways=True,
                              min_shared_genes=1):
    """
    Generate an interactive gene-pathway network using Pathfinder.R output and DE results.
    Parameters:
    - df_de: DataFrame of DE results (must contain 'gene_symbol', 'log2FoldChange', 'padj')
    - df_path: DataFrame of Pathfinder results (must contain 'Term_Description', 'lowest_p',
               'Up_regulated', and 'Down_regulated')
    - output_html: filename of the output .html file
    - highlight_terms: optional list of pathway terms to highlight
    - connect_pathways: whether to connect pathways based on shared gene similarity
    - min_shared_genes: minimum number of shared genes to connect two pathways
    """
    import pandas as pd
    import networkx as nx
    from pyvis.network import Network
    import os, webbrowser

    # 01 Load input
    if isinstance(df_de, str):
        df_de = pd.read_csv(df_de)
    if isinstance(df_path, str):
        df_path = pd.read_csv(df_path)

    # 02 Validate
    for col in ["gene_symbol", "log2FoldChange", "padj"]:
        if col not in df_de.columns:
            raise ValueError(f"Missing column '{col}' in DE dataframe.")
    for col in ["Term_Description", "lowest_p"]:
        if col not in df_path.columns:
            raise ValueError(f"Missing column '{col}' in pathfindR dataframe.")

    # 03 Build graph
    df_path_sig = df_path[df_path["lowest_p"] < 0.05].sort_values("lowest_p")
    df_sig = df_de[df_de["padj"] < 0.05]
    dict_logfc = dict(zip(df_sig["gene_symbol"], df_sig["log2FoldChange"]))

    highlight_terms = set(highlight_terms or [])
    G = nx.Graph()
    pathway_gene_sets = {}

    for _, row in df_path_sig.iterrows():
        term = row["Term_Description"]
        up_genes = str(row.get("Up_regulated", "")).split(",") if pd.notna(row.get("Up_regulated")) else []
        down_genes = str(row.get("Down_regulated", "")).split(",") if pd.notna(row.get("Down_regulated")) else []
        genes = [g.strip() for g in up_genes + down_genes if g.strip()]
        if not genes:
            continue
        pathway_gene_sets[term] = set(genes)
        node_color = "red" if term in highlight_terms else "lightgray"
        G.add_node(term, type="pathway", color={"border": "black", "background": node_color}, shape="dot")
        for gene in genes:
            fc = dict_logfc.get(gene, 0)
            color = "green" if fc > 0 else "red" if fc < 0 else "gray"
            G.add_node(gene, type="gene", color=color, shape="dot")
            G.add_edge(term, gene, color="lightgray", width=1)

    if connect_pathways:
        pathway_list = list(pathway_gene_sets.keys())
        for i in range(len(pathway_list)):
            for j in range(i + 1, len(pathway_list)):
                p1, p2 = pathway_list[i], pathway_list[j]
                genes1, genes2 = pathway_gene_sets[p1], pathway_gene_sets[p2]
                intersection = genes1 & genes2
                union = genes1 | genes2
                if len(intersection) >= min_shared_genes:
                    jaccard_sim = len(intersection) / len(union)
                    spring_len = int((1 - jaccard_sim) * 300)
                    tooltip = f"Jaccard similarity: {jaccard_sim:.2f} (Shared: {len(intersection)})"
                    G.add_edge(p1, p2, length=spring_len, title=tooltip, color="black", width=2)

    net = Network(height="800px", width="100%", notebook=False, directed=False)
    net.force_atlas_2based()

    for node, data in G.nodes(data=True):
        size = 15 + abs(dict_logfc.get(node, 0)) * 5 if data["type"] == "gene" else 30
        net.add_node(node, label=node, color=data.get("color", "gray"), size=size, 
                     title=f"{node}<br>Type: {data['type']}", shape=data.get("shape", "dot"), 
                     borderWidth=data.get("borderWidth", 1), node_type=data["type"])

    for source, target, attr in G.edges(data=True):
        net.add_edge(source, target,
                     length=attr.get("length"),
                     title=attr.get("title", ""),
                     color=attr.get("color", "gray"),
                     width=attr.get("width", 1))

    safe_filename = os.path.basename(output_html)
    net.write_html(safe_filename)
    abs_path = os.path.abspath(safe_filename)

    # 04 Inject filtering UI
    pathway_options = "\n".join([
        f"<div><label><input type='checkbox' name='pathway' value='{p}'> {p}</label></div>" 
        for p in pathway_gene_sets.keys()
    ])
    menu_html = f"""<div id="menu" style="position: absolute; top: 10px; left: 10px; z-index: 10000; background-color: white; padding: 10px; border: 1px solid lightgray; max-height: 400px; overflow-y: auto;">
    <b>Pathway List</b><br>{pathway_options}<br>
    <button onclick="filterNetwork()">Filter</button>
    <button onclick="resetNetwork()">Reset</button>
    </div>
    <script>
    function filterNetwork() {{
        const selected = Array.from(document.querySelectorAll("input[name='pathway']:checked")).map(c => c.value);
        nodes.get().forEach(node => {{
            if (node.node_type === "pathway") {{
                nodes.update({{ id: node.id, hidden: selected.indexOf(node.id) === -1 }});
            }}
        }});
        nodes.get().forEach(node => {{
            if (node.node_type === "gene") {{
                const connected = network.getConnectedEdges(node.id);
                let showGene = false;
                for (const eid of connected) {{
                    const edge = edges.get(eid);
                    const otherId = edge.from === node.id ? edge.to : edge.from;
                    const other = nodes.get(otherId);
                    if (other && other.node_type === "pathway" && !other.hidden) {{
                        showGene = true; break;
                    }}
                }}
                nodes.update({{ id: node.id, hidden: !showGene }});
            }}
        }});
    }}
    function resetNetwork() {{
        nodes.get().forEach(node => nodes.update({{ id: node.id, hidden: false }}));
        document.querySelectorAll("input[name='pathway']").forEach(cb => cb.checked = false);
    }}
    </script>"""

    with open(abs_path, "r", encoding="utf-8") as f:
        html = f.read()
    with open(abs_path, "w", encoding="utf-8") as f:
        f.write(html.replace("<body>", "<body>\n" + menu_html, 1))

    print(f"✅ Network HTML saved to: {abs_path}")

#%% 01 run function
# 01 run function
import pandas as pd

#%% B-V_vs_S-V
# File paths
df_de = pd.read_csv('DE_B-V_vs_S-V_pydeseq2_DEtool_pathfind_input.csv')
df_pathfinder = pd.read_csv("pathfind_results_kegg_B-V_vs_S-V_20250411.csv")

# Rename columns to match expected format
df_de = df_de.rename(columns={
    "Gene.symbol": "gene_symbol",
    "logFC": "log2FoldChange",
    "adj.P.Val": "padj"
})

# Optional: highlight specific pathways
highlight_terms = ['Focal adhesion',
                   'ECM-receptor interaction',
                   'cGMP-PKG signaling pathway',
                   'cAMP signaling pathway',
                   'Wnt signaling pathway',
                   'Prolactin signaling pathway'
                   ]

# Call the function to generate interactive HTML network
create_pathfinder_network(df_de, df_pathfinder, 
                          output_html="pathfindr_network_kegg_B-V_vs_S-V_0412.html", 
                          highlight_terms=highlight_terms,
                          connect_pathways=True,
                          min_shared_genes=1
                          )

#%% S-PTH_vs_S-V
df_de = pd.read_csv('DE_S-PTH_vs_S-V_pydeseq2_DEtool_pathfind_input.csv')
df_pathfinder = pd.read_csv("pathfind_results_kegg_S-PTH_vs_S-V_20250411.csv")

df_de = df_de.rename(columns={
    "Gene.symbol": "gene_symbol",
    "logFC": "log2FoldChange",
    "adj.P.Val": "padj"
})

highlight_terms = ['Parathyroid hormone synthesis, secretion and action',
                   'JAK-STAT signaling pathway',
                   'Growth hormone synthesis, secretion and action',
                   'TNF signaling pathway',
                   'Osteoclast differentiation',
                   'Th17 cell differentiation',
                   'cGMP-PKG signaling pathway',
                   'Endocrine and other factor-regulated calcium reabsorption',
                   'VEGF signaling pathway',
                   'MAPK signaling pathway'
                   ]

create_pathfinder_network(df_de, df_pathfinder, 
                          output_html="pathfindr_network_kegg_S-PTH_vs_S-V_0412.html", 
                          highlight_terms=highlight_terms,
                          connect_pathways=True,
                          min_shared_genes=1
                          )

#%% B-PTH_vs_B-V
df_de = pd.read_csv('DE_B-PTH_vs_B-V_pydeseq2_DEtool_pathfind_input.csv')
df_pathfinder = pd.read_csv("pathfind_results_kegg_B-PTH_vs_B-V_20250411.csv")

df_de = df_de.rename(columns={
    "Gene.symbol": "gene_symbol",
    "logFC": "log2FoldChange",
    "adj.P.Val": "padj"
})

highlight_terms = ['Parathyroid hormone synthesis, secretion and action',
                   'Th17 cell differentiation',
                   'JAK-STAT signaling pathway',
                   'MAPK signaling pathway',
                   'Th1 and Th2 cell differentiation',
                   'AGE-RAGE signaling pathway in diabetic complications',
                   'Rheumatoid arthritis',
                   'Toll-like receptor signaling pathway',
                   'HIF-1 signaling pathway',
                   'Prolactin signaling pathway',
                   'cGMP-PKG signaling pathway',
                   'AMPK signaling pathway',
                   'Wnt signaling pathway'
                   ]

create_pathfinder_network(df_de, df_pathfinder, 
                          output_html="pathfindr_network_kegg_B-PTH_vs_B-V_0412.html", 
                          highlight_terms=highlight_terms,
                          connect_pathways=True,
                          min_shared_genes=1
                          )
#%% B-PTH_vs_S-PTH adj-B-KD

df_de = pd.read_csv('DE_B-PTH_vs_S-PTH_pydeseq2_adj-B-KD_DEtool_pathfind_input.csv')
df_pathfinder = pd.read_csv("pathfind_results_kegg_B-PTH_vs_S-PTH_adj-B-KD_20250403.csv")

df_de = df_de.rename(columns={
    "Gene.symbol": "gene_symbol",
    "logFC": "log2FoldChange",
    "adj.P.Val": "padj"
})

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

create_pathfinder_network(df_de, df_pathfinder, 
                          output_html="pathfindr_network_kegg_B-PTH_vs_S-PTH_adj-B-KD_0411.html", 
                          highlight_terms=highlight_terms,
                          connect_pathways=True,
                          min_shared_genes=1
                          )
del df_de, df_pathfinder
# %%

