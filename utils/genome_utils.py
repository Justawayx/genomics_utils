from collections import defaultdict
from urllib.parse import unquote
import pickle

from collections import defaultdict
from urllib.parse import unquote

# Helper functions

def parse_info(info_str): # Parses INFO string of format key1=val1;key2=val2;...
	info_dict = {}
	for key_value_pair in info_str.strip('\n').split(';'):
		key, value = key_value_pair.split('=')
		info_dict[key] = value
	return info_dict

def parse_GFF(gff_fpath): # Parses any GFF3 file
	chrom_feature_data_dict = defaultdict(dict) # chrom -> (feature type, feature ID) -> (start, end, strand_direction, info_dict)	
	for line in open(gff_fpath, 'r'):
		if line.strip() == '##FASTA':
			break
		if line[0] == '#':
			continue
		chrom, source, feature_type, start_pos, end_pos, _, sdirection, _, info = line.strip().split('\t')
		start_pos = int(start_pos); end_pos = int(end_pos); info_dict = parse_info(info); feature_id = info_dict['ID']			
		chrom_feature_data_dict[chrom][(feature_type, feature_id)] = (start_pos, end_pos, sdirection, info_dict)	
	return chrom_feature_data_dict

def get_feature_desc_dict(chrom_feature_data_dict, desired_feature_type):
	feature_desc_dict = {}	
	for chrom in chrom_feature_data_dict:
		for feature_type, feature_id in chrom_feature_data_dict[chrom]:
			start, end, strand_direction, info_dict = chrom_feature_data_dict[chrom][(feature_type, feature_id)]
			if feature_type == desired_feature_type:
				feature_desc_dict[feature_id] = unquote(info_dict["description"]).replace('+', ' ')
	return feature_desc_dict

SPECIES_SNPEFF_ID_DICT = {
	("p_fal", "3D7"): "Pf3D7v3",
	("t_cru", "SylvioX10"): "TcSylvioX10v67",
}

'''
usage for now
x = test.GenomeAnnotation('/tscc/projects/ps-winzelerlab', 'p_fal', '3D7')
protein_aa_info_dict = x.alphafold_data
aa_codon_usage_ab_dict = x.aa_codon_usage_ab_dict
gene_desc_dict = x.get_gene_desc_dict()
gene_GO_dict = x.gene_GO_dict
gene_abbr_dict = x.gene_abbr_dict
gene_MIS_dict = x.gene_MIS_dict
gene_MFS_dict = x.gene_MFS_dict
'''

# Strain-specific genome annotation information

class GenomeAnnotation:
	
	def __init__(self, WDIR, species_abbr, strain):
		
		# Ex: GenomeAnnotation('/projects/winzeler', 'p_fal', '3D7')		
		self.WDIR = WDIR
		self.species = species_abbr
		self.strain = strain
		
		self.GFF_PATH_MAP = {
			("p_fal", "3D7"): f"{WDIR}/GENOME_RESOURCES/p_fal/p_fal_ref/p_fal.gff",
			("t_cru", "SylvioX10"): f"{WDIR}/GENOME_RESOURCES/t_cru/TcSylvioX10-1_TriTypDB_67/TriTrypDB-67_TcruziSylvioX10-1.gff",
		}
		
		self.chrom_feature_data_dict = self.load_parsed_GFF()
	
	def load_alphafold_data(self):
		protein_aa_info_dict = {}
		if self.species == 'p_fal' and self.strain == '3D7':
			protein_aa_info_dict = pickle.load(open('%s/PROJECTS/daisy/REF_DATA/p_fal/protein_aa_info_dict.pkl' % self.WDIR, 'rb'))
		
		return protein_aa_info_dict
	
	def load_parsed_GFF(self):		
		gff_fpath = self.GFF_PATH_MAP[(self.species, self.strain)]
		return parse_GFF(gff_fpath)
	
	def load_mutagenesis_data(self):
		f = open('%s/PROJECTS/daisy/REF_DATA/p_fal/NIHMS1004827-supplement-Table_S5.txt' % self.WDIR, 'rt')
		gene_MIS_dict = {}
		gene_MFS_dict = {}
		f.readline(); header = f.readline()
		for line in f:
				items = line.strip('\n').split('\t')
				chrom, gene_id, gene_desc, gene_ID, MIS, MFS = items[:6]
				gene_MIS_dict[gene_id] = float(MIS)
				gene_MFS_dict[gene_id] = float(MFS)
		
		return gene_MIS_dict, gene_MFS_dict
	
	def load_GO_data(self):
		
		f = open('%s/PROJECTS/daisy/REF_DATA/GO/go-basic.obo' % self.WDIR, 'r')
		GO_term_info_dict = {}
		term_start = False
		for line in f:
				if line.startswith('[Term]'):
						term_start = True
						key_value_dict = defaultdict(list)
				elif term_start:        
						if line.strip() == '':
								if 'is_obsolete' not in key_value_dict or key_value_dict['is_obsolete'] != 'True':
										GO_id = key_value_dict['id'][0]
										GO_term_info_dict[GO_id] = (key_value_dict['name'][0], key_value_dict['def'][0])
										if 'alt_id' in key_value_dict:
												for alt_id in key_value_dict['alt_id']:
														GO_term_info_dict[alt_id] = (key_value_dict['name'][0], key_value_dict['def'][0])
								term_start = False
						else:
								items = line.strip('\n').split(': ')
								key = items[0]; value = ': '.join(items[1:])
								key_value_dict[key].append(value)
		
		f = open('%s/GENOME_RESOURCES/p_fal/Pf%s_PlasmoDB_66/PlasmoDB-66_Pfalciparum%s_GO.gaf' % (self.WDIR, self.strain, self.strain), 'r')
		f.readline()
		
		gene_id_GO_dict = defaultdict(dict) # gene ID -> (GO ID, evidence_code, GO name, GO def)
		gene_id_symbol_dict = defaultdict(dict) # gene ID -> gene abbr
		
		for line in f:
				items = line.strip('\n').split('\t')
				db, object_id, object_symbol, qualifier, GO_id, db_reference, evidence_code, with_from, aspect, \
						object_name, object_synonym, object_type, taxon, date, assigned_by, _, _ = items
				GO_name, GO_def = GO_term_info_dict[GO_id]
				gene_id_GO_dict[object_id] = (GO_id, evidence_code, GO_name, GO_def)
				gene_id_symbol_dict[object_id] = object_symbol
		
		return gene_id_GO_dict, gene_id_symbol_dict
	
	def load_codon_usage(self):
		
		f = open("%s/GENOME_RESOURCES/p_fal/Pf%s_PlasmoDB_66/PlasmoDB-66_Pfalciparum%s_CodonUsage.txt" % (self.WDIR, self.strain, self.strain), 'r')
		header_items = f.readline()
		aa_codon_usage_ab_dict = defaultdict(dict)
		for line in f:
				codon, aa, freq, abundance = line.strip('\n').split('\t')
				freq = float(freq); abundance = float(abundance)
				aa_codon_usage_ab_dict[aa][codon] = abundance
		
		return aa_codon_usage_ab_dict
	
	def get_gene_desc_dict(self):		
		return get_feature_desc_dict(self.chrom_feature_data_dict, 'protein_coding_gene')
	
	def get_chromosomes(self):
		return sorted(self.chrom_feature_data_dict.keys())

# General utility functions
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
