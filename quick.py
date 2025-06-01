import pandas as pd

# this is just a very quick script to extract the top representative gene for each GO and KEGG term based on the DEGs from the MICT vs M comparison

degs = pd.read_csv('ORA_results.csv')
degs = degs[(degs['comparison'] == 'MICT_vs_M') & (degs['padj'] < 0.05)]

# ----------------------- GO_BP -----------------------

go = pd.read_csv('GO_BP_enrichment.csv')
go_records = []
for _, row in go.iterrows():
    term = row['Term']
    genes = row['Genes'].split(';') 
    overlap_degs = degs[degs['level_1'].isin(genes)].copy()
    
    if overlap_degs.empty:
        continue
    
    overlap_degs['absFC'] = overlap_degs['log2FoldChange'].abs()
    top_deg = overlap_degs.sort_values('absFC', ascending = False).iloc[0]
    symbol = top_deg['level_1']
    url = f'https://www.genecards.org/cgi-bin/carddisp.pl?gene={symbol}'
    go_records.append({'GO_Term': term, 'Representative_Gene': symbol, 'log2FC': top_deg['log2FoldChange'], 'padj': top_deg['padj'], 'GeneCards_URL': url})

go_representatives = pd.DataFrame(go_records)
go_representatives.to_csv('GO_term_representatives.csv', index = False)

# ----------------------- KEGG -----------------------

kegg = pd.read_csv('KEGG_enrichment.csv')
kegg_records = []

for _, row in kegg.iterrows():
    term = row['Term']
    genes = row['Genes'].split(';') 
    overlap_degs = degs[degs['level_1'].isin(genes)].copy()
    
    if overlap_degs.empty:
        continue
    
    overlap_degs['absFC'] = overlap_degs['log2FoldChange'].abs()
    top_deg = overlap_degs.sort_values('absFC', ascending = False).iloc[0]
    symbol = top_deg['level_1']
    url = f'https://www.genecards.org/cgi-bin/carddisp.pl?gene={symbol}'
    kegg_records.append({'KEGG_Term': term, 'Representative_Gene': symbol, 'log2FC': top_deg['log2FoldChange'], 'padj': top_deg['padj'], 'GeneCards_URL': url})

kegg_representatives = pd.DataFrame(kegg_records)
kegg_representatives.to_csv('KEGG_term_representatives.csv', index = False)
