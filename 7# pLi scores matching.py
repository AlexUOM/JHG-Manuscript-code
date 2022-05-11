import pandas as pd

Proteins_file = pd.read_csv('lists of prots with High impact.csv')
pLi_file = pd.read_table('pLi data.txt')

# Putting all proteins with high impact variants (from BioGRID analysis)
# into a list
all_prots = []
for protein in Proteins_file.iterrows():
    all_prots.append(protein[1]["High_impact_proteins"])


output_dict = {'Proteins': [], 'pLi': []}
for protein in all_prots:
    output_dict['Proteins'].append(protein)
    output_dict['pLi'].append(
        pLi_file.loc[pLi_file['gene'] == protein].pLI.item())

output_df = pd.DataFrame(output_dict)
output_df.to_excel('Proteins and pLis.xlsx', index=False)
