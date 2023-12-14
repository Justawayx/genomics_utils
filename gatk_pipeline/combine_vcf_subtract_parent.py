from collections import defaultdict
import os, sys
from urllib.parse import unquote
from utils import genome_utils

# =======================
# Parameters
# =======================

CONFIG_PATH = os.path.abspath('config.cfg')

def get_env_value_from_config(env_variable):
	for line in open(CONFIG_PATH, 'r'):
		if line.startswith('#'):
			continue
		if '=' in line:
			key, value = line.rstrip('\n').split('=')
			if key == env_variable:
				return value
	print("Environmental variable %s not found in config.cfg" % env_variable)
	return False

WDIR = "/projects/winzeler"
SAMPLES_PATH = get_env_value_from_config('samples_file')
MAIN_DIR = get_env_value_from_config('main_dir')
SPECIES_ABBR = get_env_value_from_config('species_abbr')
GROUP_NAME = get_env_value_from_config('group_name')

samples = [line.strip() for line in open(SAMPLES_PATH, 'r')]
print("Here are the sample names:\n\t" + '\n\t'.join(samples))
answer = input("Enter name of the parent sample: ")
while answer not in samples:
	answer = input("%s is not in samples.txt, please try again: " % answer)

PARENT_SAMPLE = answer

SNPEFF_DIR = f"{WDIR}/GENOME_RESOURCES/snpEff"

def get_format_data_dict(FORMAT, data):
	return {key: val for key, val in zip(FORMAT.split(':'), data.split(':'))}

def parse_info(INFO):
	items = INFO.split(';')
	info_dict = {}
	for item in items:
		key, value = item.split('=')
		info_dict[key] = value
	return info_dict

# =======================
# Load gene annotations
# =======================

species_refstrain_dict = {'p_fal': '3D7'}

ga = genome_utils.GenomeAnnotation(WDIR, SPECIES_ABBR, species_refstrain_dict[SPECIES_ABBR])
gene_desc_dict = ga.get_gene_desc_dict()

# =======================
# Combine VCFs
# =======================

vcf_path = f"{MAIN_DIR}/{GROUP_NAME}.raw.snps.indels.vcf"
vcf_exists_flag = os.path.isfile(vcf_path)

if vcf_exists_flag:
	answer = input(f"{GROUP_NAME}.raw.snps.indels.vcf already exists (size: {os.path.getsize(vcf_path)//1e6} MB), combine VCFs? y/n ")
	while answer.lower() not in ['y', 'n']:
		answer = input("Please enter y or n: ")
	
	if answer.lower() == 'y':
		vcf_exists_flag = False # If yes, pretend it doesn't exist

if not vcf_exists_flag:
	VCFS = [fname for fname in os.listdir(MAIN_DIR) if (fname.startswith(GROUP_NAME + '-part') and fname.endswith(".raw.snps.indels.vcf"))]
	print("The following partial VCF(s) were found in main_dir:\n" + '\n'.join(VCFS))
	
	answer = input("Proceed with combining? y/n ")
	while answer.lower() not in ['y', 'n']:
		answer = input("Please enter y or n: ")
	
	if answer.lower() == 'n':
		sys.exit("Aborting")
	
	if len(VCFS) > 1:
		vcf_partnum_tups = [(vcf, int(vcf.split('-part')[1].split('.raw.snps.indels.vcf')[0])) for vcf in VCFS]
		ordered_vcfs = [tup[0] for tup in sorted(vcf_partnum_tups, key=lambda x: x[1])]
		
		o = open(f"{MAIN_DIR}/{GROUP_NAME}.raw.snps.indels.vcf", 'w')
		first_vcf = ordered_vcfs[0]
		with open(f"{MAIN_DIR}/{first_vcf}", 'r') as f:
			for line in f:
				o.write(line)
		
		for vcf in ordered_vcfs[1:]:
			with open(f"{MAIN_DIR}/{vcf}", 'r') as f:
				for line in f:
					if line.startswith("#CHROM"):
						break
				
				for line in f:
					o.write(line)
	
	o.close()

# =======================
# Run SnpEff on VCFs
# =======================

VCF = f"{GROUP_NAME}.raw.snps.indels.vcf"

vcf_ann_path = f"{MAIN_DIR}/{VCF}.ann.txt"
ann_exists_flag = os.path.isfile(vcf_ann_path)

if ann_exists_flag:
	answer = input(f"{VCF}.ann.txt already exists (size: {os.path.getsize(vcf_ann_path)//1e6} MB), rerun SnpEff? y/n ")
	while answer.lower() not in ['y', 'n']:
		answer = input("Please enter y or n: ")
	
	if answer.lower() == 'y':
		ann_exists_flag = False # If yes, pretend it doesn't exist

if not ann_exists_flag:
	SNPEFF_SPECIES_ID = genome_utils.SPECIES_SNPEFF_ID_DICT[SPECIES_ABBR]
	os.chdir(MAIN_DIR)
	cmd = f"java -jar {SNPEFF_DIR}/snpEff.jar -formatEFF -o vcf -ud 1000 -c {SNPEFF_DIR}/snpEff.config {SNPEFF_SPECIES_ID} {MAIN_DIR}/{VCF} > {MAIN_DIR}/{VCF}.ann.txt"
	print(cmd)
	os.system(cmd)

# =======================
# Load data from VCF
# =======================

f = open(f"{MAIN_DIR}/{VCF}.ann.txt")

child_samples = []
for line in f:
	if line.startswith('#CHROM'):
		header_items = line.strip().split('\t')
		for item in header_items[9:]:
			if item != PARENT_SAMPLE:
				child_samples.append(item)
		break

if set(child_samples + [PARENT_SAMPLE]) != set(samples):
	sys.exit("Samples in VCF are inconsistent with the samples in samples.txt, aborting")

chrom_pos_sample_data_dict = defaultdict(dict) # chrom -> pos -> sample -> data
chrom_pos_anno_dict = defaultdict(dict) # chrom -> pos -> (ref, alt, gene, type, effect, impact, codon_change, aa_change)

revert_sites = set() # (chrom, pos) tuples

chrom_pos_quality_dict = defaultdict(dict)
parent_col_idx = header_items.index(PARENT_SAMPLE)

for line in f:
		items = line.strip().split('\t')
		CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = items[:9]
		POS = int(POS)
		
		if float(QUAL) < 500:
				continue
		
		parent_data = items[parent_col_idx]
		parent_data_dict = get_format_data_dict(FORMAT, parent_data)
		
		if parent_data_dict['GT'] == './.' or 'DP' not in parent_data_dict or \
			 parent_data_dict['DP'] == '.' or int(parent_data_dict['DP']) < 7:
				continue
		
		info_dict = parse_info(INFO)
		if ("ReadPosRankSum" in info_dict and (float(info_dict["ReadPosRankSum"]) > 8 or float(info_dict["ReadPosRankSum"]) < -8)) \
		 or ("QD" in info_dict and float(info_dict["QD"]) < 2) \
		 or ("MQRankSum" in info_dict and float(info_dict["MQRankSum"]) < -12.5):
				continue
		
		chrom_pos_quality_dict[CHROM][POS] = QUAL
		chrom_pos_sample_data_dict[CHROM][POS] = {}
		non_parent_GTs = []
		
		for idx in range(9, len(items)):
			sample = header_items[idx]
			data = items[idx]
			data_dict = get_format_data_dict(FORMAT, data)
			GT = data_dict['GT']
			if idx != parent_col_idx:
				non_parent_GTs.append(GT)
			if 'DP' not in data_dict or data_dict['DP'] == '.' or int(data_dict['DP']) < 7:
				chrom_pos_sample_data_dict[CHROM][POS][sample] = (data_dict['GT'], 'LowDP')
			else:
				AD = data_dict['AD']
				chrom_pos_sample_data_dict[CHROM][POS][sample] = (GT, AD)
		
		if set(non_parent_GTs) == set(['0/0']):
			revert_sites.add((CHROM, POS))
		
		first_alt_allele = ALT.split(',')[0]
		
		pre_eff = INFO.split(';EFF=')[0]
		eff = INFO.split(';EFF=')[1].split(';')[0]
		
		allele_anno_dict = {}
		for anno in eff.split(','):
				last = anno.split('|')[-1]
				if 'INFO_' in last or 'WARNING_' in last:
						alt_allele = anno.split('|')[-2].split(')')[0]
				else:
						alt_allele = anno.split('|')[-1].split(')')[0]
				effect, rest = anno.split('(')
				impact, codon_change, aa_change, _, gene = rest.split('|')[1:6] # Note that codon_change is distance for intergenics
				if alt_allele in allele_anno_dict: # Already has its own annotation
						if allele_anno_dict[alt_allele][0] == 'intergenic_region': # Replace intergenic
								allele_anno_dict[alt_allele] = (effect, impact, codon_change, aa_change, gene)
				else:
						allele_anno_dict[alt_allele] = (effect, impact, codon_change, aa_change, gene)
		
		vtypes = []; effects = []; impacts = []; codon_changes = []; aa_changes = []
		for alt_allele in ALT.split(','):
				if alt_allele == '*':
					continue
				effect, impact, codon_change, aa_change, gene = allele_anno_dict[alt_allele]
				vtypes.append('SNP' if (REF in ['A', 'C', 'G', 'T'] and alt_allele in ['A', 'C', 'G', 'T']) else 'INDEL')
				effects.append(effect); impacts.append(impact); codon_changes.append(codon_change); aa_changes.append(aa_change)
		
		chrom_pos_anno_dict[CHROM][POS] = (REF, ALT, gene, ','.join(vtypes), ','.join(effects), ','.join(impacts), ','.join(codon_changes), ','.join(aa_changes))

# =======================
# Write SNV/indel table
# =======================

def AD_to_AAF(AD):
		if AD == 'LowDP':
				return 'LowDP'
		depths = [int(item) for item in AD.split(',')]
		ref_depth = depths[0]
		alt_depths = depths[1:]
		total_depth = sum(depths)
		return sum(alt_depths)/total_depth

ordered_samples = [PARENT_SAMPLE] + sorted(child_samples)

f = open(f"{MAIN_DIR}/{GROUP_NAME}_GATK-Filters.tsv", 'w')
header_items = ['Chromosome', 'Position', 'Gene_Name', 'Gene_Descrip', 'Quality', 'Ref_Base', 'Alt_Base', 'Type', 
								 'Effect', 'Impact', 'Codon_Change', 'Amino_Acid_Change']

for sample in ordered_samples:
		header_items.append("%s_MutationCall" % sample)

for sample in ordered_samples:
		header_items.append("%s_AltFrequency" % sample)

for sample in ordered_samples:
		header_items.append("%s_AlleleDepths" % sample)

header_items.append("Has_NonHet")

f.write('\t'.join(header_items) + '\n')

for chrom in chrom_pos_sample_data_dict:
		for pos in sorted(chrom_pos_sample_data_dict[chrom]):
				
				if (chrom, pos) in revert_sites:
						pass
				elif chrom_pos_sample_data_dict[chrom][pos][PARENT_SAMPLE][0] != '0/0':
						continue
				
				quality = chrom_pos_quality_dict[chrom][pos]
				ref, alt, gene, vtype, effect, impact, codon_change, aa_change = chrom_pos_anno_dict[chrom][pos]
				gene_desc = gene_desc_dict[gene] if gene in gene_desc_dict else ''
				
				items = [chrom, pos, gene, gene_desc, quality, ref, alt, vtype, effect, impact, codon_change, aa_change]
				
				GTs = []
				for sample in ordered_samples:
						GT = chrom_pos_sample_data_dict[chrom][pos][sample][0]
						items.append(' ' + GT); GTs.append(GT)
				
				ADs = []
				for sample in ordered_samples:
						AD = chrom_pos_sample_data_dict[chrom][pos][sample][1]
						AAF = AD_to_AAF(AD)
						items.append(AAF); ADs.append(AD)
				
				for sample in ordered_samples:
						AD = chrom_pos_sample_data_dict[chrom][pos][sample][1]
						items.append(AD)
				
				has_non_het = False
				for AD, GT in zip(ADs, GTs):
						if not GT.startswith('0/') and AD != "LowDP":
								has_non_het = True
				
				items.append(has_non_het)
				f.write('\t'.join([str(item) for item in items]) + '\n')

f.close()

print(f"Done! Check {MAIN_DIR}/{GROUP_NAME}_GATK-Filters.tsv")
