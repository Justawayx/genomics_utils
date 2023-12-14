from collections import defaultdict
from urllib.parse import unquote
import pickle

SPECIES_SNPEFF_ID_DICT = {"p_fal": "Pf3D7v3"}

def get_protein_aa_info_dict(WDIR):
    protein_aa_info_dict = pickle.load(open('%s/ROTATION_PROJECT/daisy/REF_DATA/pfal/protein_aa_info_dict.pkl' % WDIR, 'rb'))
    return protein_aa_info_dict

def get_gene_info_dicts(WDIR, species="p_fal"):
    
    if species == "p_fal":
        
        gff_fpath = "%s/GENOME_RESOURCES/pf/p_fal_ref/p_fal.gff" % WDIR
        
        chromosomes = []
        for line in open(gff_fpath, 'r'):
                if line[0] != '#':
                        break
                if line[:17] == '##sequence-region':
                        chromosomes.append(line.strip().split()[1])
        
        chrom_gene_exon_interval_dict = {chrom: defaultdict(dict) for chrom in chromosomes} # chrom -> gene ID -> exon ID -> (start, end, strand_direction)
        
        chrom_gene_ids_dict = {chrom: set() for chrom in chromosomes} # All protein coding gene IDs
        gene_desc_dict = {chrom: {} for chrom in chromosomes} # gene_id -> description
        gene_interval_dict = {chrom: {} for chrom in chromosomes} # gene_id -> (start, end, strand_direction)
        
        for line in open(gff_fpath, 'r'):
                if line.strip() == '##FASTA':
                        break
                if line[0] == '#':
                        continue
                chrom, source, feature_type, start_pos, end_pos, _, sdirection, _, info = line.strip().split('\t')
                start_pos = int(start_pos); end_pos = int(end_pos)
                if feature_type == 'exon':
                        exon_id = info.split(';')[0].split('ID=')[1]
                        gene_id = exon_id.split('exon_')[1].split('-')[0]
                        chrom_gene_exon_interval_dict[chrom][gene_id][exon_id] = (start_pos, end_pos, sdirection)
                if feature_type == 'gene':
                        gene_id = info.split(';')[0].split('ID=')[1]
                        gene_desc = info.split(';')[2].split('description=')[1]
                        gene_desc_dict[gene_id] = unquote(gene_desc.strip()).replace('+', ' ')
                        gene_interval_dict[gene_id] = (start_pos, end_pos, sdirection)
                if feature_type == 'CDS': # Actually coding
                        gene_id = info.split(';')[0].split('ID=')[1].split('cds_')[1].split('-')[0]
                        if gene_id != 'PF3D7_0112400' and 'pseudogene' not in gene_desc_dict[gene_id]: # Ignore pseudogenes
                                chrom_gene_ids_dict[chrom].add(gene_id)
        
        f = open('%s/ROTATION_PROJECT/daisy/Winzeler_databases/PlasmDBv29AnnotationsforMZstudy.txt' % WDIR, 'r')
        header_items = f.readline().strip('\n').split('\t')
        
        gene_class_dict = {}
        
        possible_target_classes = set()
        likely_target_values = set()
        nontarget_classes = set()
        
        for line in f:
            items = line.strip('\n').split('\t')
            gene_id = items[0]
            nonessential = items[header_items.index('Nonessential')]
            class_field = items[header_items.index('Class')]
            category_field = items[header_items.index('Category')]
            target_class = items[header_items.index('Target Class')]
            likely_target = items[header_items.index('Field 24')]
            if likely_target == 'no' or target_class == 'Not Predicted Target':
                nontarget_classes.add(target_class)
            else:
                likely_target_values.add(likely_target)
                possible_target_classes.add(target_class)
            gene_class_dict[gene_id] = (nonessential, target_class, likely_target, class_field, category_field)
        
        return chrom_gene_ids_dict, gene_desc_dict, gene_interval_dict, gene_class_dict

aa3_to_aa1_dict = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Glu": "E", "Gln": "Q", "Gly": "G", 
                   "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", 
                   "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "Ter"}

def aa_change_to_abbr(aa_change):
    if 'p.' not in aa_change:
        return aa_change
    aa_change_only = aa_change.split('/')[0].split('p.')[1]
    aa1_1 = aa3_to_aa1_dict[aa_change_only[:3]]
    if aa_change_only[-3:] in aa3_to_aa1_dict:
        aa1_2 = aa3_to_aa1_dict[aa_change_only[-3:]]
        return aa1_1 + aa_change_only[3:-3] + aa1_2
    else:
        return aa1_1 + aa_change_only[3:]

known_targets = [lambda gene: gene == 'PF3D7_0603300', # DHODH
                 lambda gene: gene == 'mal_mito_3', # Cytochrome b
                 lambda gene: '--tRNA ligase' in gene_desc_dict[gene], # Various aaRS
                 lambda gene: 'proteasome' in gene_desc_dict[gene] and 'subunit' in gene_desc_dict[gene], # Proteasome
                 lambda gene: gene == 'PF3D7_0417200', # DHFR
                 lambda gene: gene == 'PF3D7_0627800', # AcAS
                 lambda gene: gene == 'PF3D7_1238800', # ACS11
                 lambda gene: gene == 'PF3D7_1211900', # ATP4
                 lambda gene: gene == 'PF3D7_1114700', # CLK3
                 lambda gene: gene == 'PF3D7_1438500', # CPSF3
                 lambda gene: gene == 'PF3D7_1443700', # DPCK
                 lambda gene: gene == 'PF3D7_1451100', # eEF2
                 lambda gene: gene == 'PF3D7_1147500', # Ftbeta
                 lambda gene: gene == 'PF3D7_0823300', # GCN5
                 lambda gene: gene == 'PF3D7_1128400', # GGPPS
                 lambda gene: gene == 'PF3D7_0204700', # HT1
                 lambda gene: gene == 'PF3D7_0107500', # NCR1
                 lambda gene: gene == 'PF3D7_1412800', # NMT1
                 lambda gene: gene == 'PF3D7_0915400', # PFK9
                 lambda gene: gene == 'PF3D7_0509800', # PI4K
                 lambda gene: gene == 'PF3D7_1436600', # PKG
                 lambda gene: gene == 'PF3D7_0808200', # PMX
                 lambda gene: gene == 'PF3D7_0810800', # PPPK-DHPS
                ]

known_target_abbr_dict = {'PF3D7_0603300': 'DHODH', 'mal_mito_3': 'CYTB', 'PF3D7_0417200': 'DHFR', 'PF3D7_0627800': 'AcAS', 
                          'PF3D7_1238800': 'ACS11', 'PF3D7_1211900': 'ATP4', 'PF3D7_1114700': 'CLK3', 'PF3D7_1438500': 'CPSF3', 
                          'PF3D7_1443700': 'DPCK', 'PF3D7_1451100': 'eEF2', 'PF3D7_1147500': 'Ftbeta', 'PF3D7_0823300': 'GCN5', 
                          'PF3D7_1128400': 'GGPPS', 'PF3D7_0204700': 'HT1', 'PF3D7_0107500': 'NCR1', 'PF3D7_1412800': 'NMT1', 
                          'PF3D7_0915400': 'PFK9', 'PF3D7_0509800': 'PI4K', 'PF3D7_1436600': 'PKG', 'PF3D7_0808200': 'PMX', 
                          'PF3D7_0810800': 'PPPK-DHPS'}

def is_known_target(gene):
    if gene not in gene_desc_dict:
        return False
    for target_func in known_targets:
        if target_func(gene) is True:
            return True
    return False

def aa_change_to_idx(aa_change):
    num_str = ''
    first_num_seen = False
    for char in aa_change:
        if char.isnumeric():
            first_num_seen = True
            num_str += char
        elif first_num_seen:
            break
    
    return int(num_str)

def is_multigene_family(gene):
    for keyword in ['PfEMP1', 'rifin', 'stevor']:
        if keyword in gene_desc_dict[gene]:
            return True
    return False
