from collections import defaultdict
import pandas as pd, pickle

# Gene information

pfal_ref_dir = '/storage/NFS/ROTATION_PROJECT/daisy/REF_DATA/pfal/'

df = pd.read_csv('%s/Pf3D7_Gene_Annotations_FULL_v14.csv' % pfal_ref_dir)

new_PF3D7_genes = set(df['Gene ID'])
new_gene_desc_dict = {}; new_gene_symbol_dict = {}
new_symbol_or_desc_gene_dict = defaultdict(list)
gene_column_dict = defaultdict(dict)

for i, row in df.iterrows():
    gene = row['Gene ID']
    for column in list(df.columns):
        gene_column_dict[gene][column] = row[column]
    new_gene_desc_dict[gene] = row['Product Description']
    new_gene_symbol_dict[gene] = '' if pd.isna(row['Gene Symbol']) else row['Gene Symbol']
    if not pd.isna(row['Gene Symbol']):
        new_symbol_or_desc_gene_dict[row['Gene Symbol']].append(gene)
    elif row['Product Description'] != '':
        new_symbol_or_desc_gene_dict[row['Product Description']].append(gene)

print("Total number of protein coding genes (release 66): %i" % len(new_PF3D7_genes))

gene_essentiality_class_dict = pickle.load(open('%s/gene_essentiality_class_dict.pkl' % pfal_ref_dir, 'rb'))
gene_binding_class_dict = pickle.load(open('%s/gene_binding_class_dict.pkl' % pfal_ref_dir, 'rb'))
