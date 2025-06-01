import pandas as pd

# this is just a very quick script to extract the top representative gene for each GO and KEGG term based on the DEGs from the MICT vs M comparison

degs = pd.read_csv('ORA_results.csv')
degs = degs[(degs['comparison'] == 'MICT_vs_M') & (degs['padj'] < 0.05)]

def retrieval(file):
    df = pd.read_csv(file)
    records = []

    for _, row in df.iterrows():
        term = row['Term']
        genes = row['Genes'].split(';') 
        overlap = degs[degs['level_1'].isin(genes)].copy()
        
        if overlap.empty:
            continue
        
        overlap['absFC'] = overlap['log2FoldChange'].abs()
        top = overlap.sort_values('absFC', ascending = False).iloc[0]
        symbol = top['level_1']
        url = f'https://www.genecards.org/cgi-bin/carddisp.pl?gene={symbol}'
        records.append({'GO_Term': term, 'Representative_Gene': symbol, 'log2FC': top['log2FoldChange'], 'padj': top['padj'], 'GeneCards_URL': url})

    return pd.DataFrame(records)

# ----------------------- GO_BP -----------------------

go_representatives = retrieval('GO_BP_enrichment.csv')
go_representatives.to_csv('GO_term_representatives.csv', index = False)

# ----------------------- KEGG -----------------------

kegg_representatives = retrieval('KEGG_enrichment.csv')
kegg_representatives.to_csv('KEGG_term_representatives.csv', index = False)
