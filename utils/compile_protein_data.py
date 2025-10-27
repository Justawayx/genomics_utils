# ====================================
# compile_protein_data.py
# Constructs a table of protein
# annotations (indexed by transcript)
# ====================================

from collections import defaultdict
import pandas as pd
from file_parsing_utils import parse_GFF

REF_DATA_DIR = '/storage/NFS/ROTATION_PROJECT/daisy/REF_DATA'
REF_GENOME_DIR = f'/storage/NFS/GENOME_RESOURCES/p_fal/3D7'

# Firstly get list of transcripts in Pf3D7 genome
gene_df = pd.read_csv(f'{REF_DATA_DIR}/pfal/Pf3D7_Gene_Annotations_FULL_v14.csv')
gene_df.set_index('Transcript ID', inplace=True, verify_integrity=True)

# Initialize table of protein data to build up
desired_columns = ['Transcript Product Description', 'Protein Sequence', 'Protein Length',
                   'Uniprot IDs', 'Molecular Weight', 'Isoelectric Point', 'Num TM Domains', 'PDB IDs']

PROTEIN_INFO_DICT = gene_df[desired_columns].to_dict(orient='index')

# Cross check with reference genome data
ftype_feature_info_dict, feature_sequence_dict = parse_GFF(f'{REF_GENOME_DIR}/p_fal.gff', ['transcript'])
for feature_ID in ftype_feature_info_dict['gene']:
    chrom, start, end, attribute_dict = ftype_feature_info_dict['gene'][feature_ID]
    print(attribute_dict)

# Keys are either transcript or gene IDs
transcript_protein_ID_dict = defaultdict(set)
transcript_protein_desc_dict = {}

# UniProt
with open(f'{REF_DATA_DIR}/pfal/uniprotkb_taxonomy_id_36329_2023_11_27.fasta', 'r') as f:
    for line in f:
        if line.startswith('>'):
            if 'GN' not in line:
                continue
            _, protein, rest = line.strip('\n').split('|')
            other_id = rest.split(' ')[0]
            desc = ' '.join(rest.split('OS=')[0].split()[1:])
            transcript = rest.split('GN=')[1].split()[0]
            if not transcript.startswith('PF3D7'): # Non-canonical gene name
                continue
            transcript_protein_ID_dict[transcript].add(protein)
            transcript_protein_desc_dict[transcript] = desc # Each transcript maps to one protein description

# PlasmoDB
with open(f'{REF_DATA_DIR}/pfal/Pf_Uniprot_PlasmoDB_idmapping_2024_01_12.tsv') as f:
    header = f.readline()
    for line in f:
        uniprot_id, veupathdb_str = line.strip('\n').split('\t')
        if '3D7' in veupathdb_str:
            gene_id = veupathdb_str.split('DB:')[1].strip('\n')
            transcript_protein_ID_dict[gene_id].add(uniprot_id)

# Manually fix descriptions
transcript_protein_desc_dict['PF3D7_API04300'] = "DNA-directed RNA polymerase subunit beta"
transcript_protein_desc_dict['PF3D7_API02900'] = "Elongation factor Tu, apicoplast"
transcript_protein_desc_dict['PF3D7_1335400'] = "Reticulocyte-binding protein homolog 2a"
transcript_protein_desc_dict['PF3D7_0316300'] = "Probable inorganic pyrophosphatase"

# Get AlphaFold structure information
protein_aa_info_dict = defaultdict(dict)

for gene in gene_protein_ID_dict:
    
    if '3D7' not in gene:
        continue
    
    try:
        f = open('/storage/NFS/ROTATION_PROJECT/daisy/REF_DATA/pfal/AlphaFold_CIFs/%s.cif' % gene, 'r')
    except:
        print("%s is missing .cif file" % gene)
        continue

    pLDDT_seen = False
    for line in f:    
        if 'pLDDT' in line:
            pLDDT_seen = True
        if pLDDT_seen and 'pLDDT' not in line:
            break

    for line in f:
        if line.startswith('#') or line.startswith('loop_') or line.startswith('_ma'):
            continue
        if line.startswith('A'):
            items = line.strip('\n').split()
            if len(items) == 7:
                _, aa3, aa_idx, _, pLDDT, _, _ = items
                protein_aa_info_dict[gene][int(aa_idx)] = (aa3, float(pLDDT))
        else:
            break