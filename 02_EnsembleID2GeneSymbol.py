#%% 
# import packages
import pandas as pd
from mygene import MyGeneInfo
import time

#%% 
# load the datasets
df_BV_SV = pd.read_csv('DE_B-V_vs_S-V_pydeseq2.csv')
df_SPTH_SV = pd.read_csv('DE_S-PTH_vs_S-V_pydeseq2.csv')
df_BPTH_BV = pd.read_csv('DE_B-PTH_vs_B-V_pydeseq2.csv')
df_BPTH_SPTH = pd.read_csv('DE_B-PTH_vs_S-PTH__pydeseq2_adj-B-KD.csv')

#%%
# define function

def ensembl2gene_symbol(df, gene_id_col='gene_id', new_col='gene_symbol', species='mouse'):
    """
    Converts Ensembl Gene IDs to Gene Symbols and adds a new column to the DataFrame.
    Parameters:
    - df (pd.DataFrame): Input DataFrame containing Ensembl Gene IDs.
    - gene_id_col (str): Column name with Ensembl Gene IDs. Default is 'gene_id'.
    - new_col (str): Column name to store the gene symbols. Default is 'gene_symbol'.
    - species (str): Species for gene mapping (default: 'human').
    Returns:
    - pd.DataFrame: DataFrame with a new column containing gene symbols.
    """
    mg = MyGeneInfo()

    # Query gene names
    results = mg.querymany(df[gene_id_col].dropna().tolist(), scopes="ensembl.gene", fields="symbol", species=species)

    # Convert results to a dictionary {Ensembl_ID: Gene_Symbol}
    id_to_gene = {entry['query']: entry.get('symbol', 'N/A') for entry in results}

    # Map the results to the DataFrame
    df[new_col] = df[gene_id_col].map(id_to_gene)

    return df
#%% 
# process. it takes about 8 minutes
start_time = time.time()

print ('processing', f'df_BV_SV')
df_BV_SV = ensembl2gene_symbol(df_BV_SV)
print ('processing', f'df_SPTH_SV')
df_SPTH_SV = ensembl2gene_symbol(df_SPTH_SV)
print ('processing', f'df_BPTH_BV')
df_BPTH_BV = ensembl2gene_symbol(df_BPTH_BV)
print ('processing', f'df_BPTH_SPTH')
df_BPTH_SPTH = ensembl2gene_symbol(df_BPTH_SPTH)

end_time = time.time()
elapsed_time = end_time - start_time
hours, remainder = divmod(elapsed_time, 3600)
minutes, seconds = divmod(remainder, 60)
formatted_time = f"{hours}h {minutes}m {seconds}s"

print(f"Execution Time: {formatted_time}")

del start_time, end_time

#%% 
# export
df_BV_SV.to_csv('DE_B-V_vs_S-V_pydeseq2.csv')
df_SPTH_SV.to_csv('DE_S-PTH_vs_S-V_pydeseq2.csv')
df_BPTH_BV.to_csv('DE_B-PTH_vs_B-V_pydeseq2.csv')
df_BPTH_SPTH.to_csv('DE_B-PTH_vs_S-PTH__pydeseq2_adj-B-KD.csv')
