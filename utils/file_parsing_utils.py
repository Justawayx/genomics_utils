from collections import defaultdict
from urllib.parse import unquote
from tqdm import tqdm
from scipy.stats import binomtest, fisher_exact, barnard_exact, boschloo_exact, chi2_contingency
import numpy as np
from matplotlib import pyplot as plt
import subprocess, gzip

# ====================
# FASTA parsing tools
# ====================

def parse_FASTA(filepath, keep_full_header=False):
    identifier_sequence_dict = defaultdict(str)
    with (gzip.open(filepath, 'r') if filepath.endswith('.gz') else open(filepath, 'r')) as f:
        first_seq = True; seq = ''
        for line in f:
            if line[0] == '>':
                if not first_seq:
                    identifier_sequence_dict[identifier] = seq
                header = line.strip('\n')[1:]
                if keep_full_header:
                    identifier = header
                else:
                    identifier = header.split()[0].split('|')[0].strip()
                first_seq = False
                seq = ''
            else:
                seq += line.strip('\n')
        identifier_sequence_dict[identifier] = seq
    return identifier_sequence_dict

def parse_FASTA_from_lines(lines):
    identifier_sequence_dict = defaultdict(str)
    first_seq = True; seq = ''
    for line in lines:
        if line[0] == '>':
            if not first_seq:
                identifier_sequence_dict[identifier] = seq
            header = line.strip('\n')[1:]
            identifier = header.split()[0].split('|')[0].strip()
            first_seq = False
            seq = ''
        else:
            seq += line.strip('\n')
    identifier_sequence_dict[identifier] = seq
    return identifier_sequence_dict


# ====================
# GenBank parsing tools
# ====================

def get_GenBank_features(filepath):
	feature_location_dict = {}
	with open(filepath, 'r') as f:
		for line in f:
			if line.startswith("FEATURES"):
				break
		for line in f:
			items = line.split()
			if line.startswith("ORIGIN"):
				break
			elif len(line.split('/')[0].split()) == 2:
				feature_type, location = items
			elif items[0].startswith('/label='):
				label = items[0].split('/label=')[1]
				feature_location_dict[label] = (location, feature_type)
	return feature_location_dict

# ====================
# GFF parsing tools
# ====================

# GFF/GTF format information:
# https://useast.ensembl.org/info/website/upload/gff.html

def parse_GFF_attribute(attribute_str, delimiter=';', key_val_sep="="):
	pairs = attribute_str.rstrip(delimiter).split(delimiter)
	attribute_dict = {}
	for pair in pairs:
		key, value = pair.strip().split(key_val_sep)
		attribute_dict[key] = value
	return attribute_dict

def parse_GFF(filepath):
	
	ftype_feature_info_dict = defaultdict(dict) # feature type -> feature ID -> (chrom, start, end, attribute_dict) 
	feature_ftype_dict = {}
	feature_sequence_dict = {} # feature ID -> sequence
	
	with (gzip.open(filepath, 'r') if filepath.endswith('.gz') else open(filepath, 'r')) as f:
		for line in f:
			if line[0] == '#':
				continue
			elif line[0] == '>':
				break
			items = line.rstrip('\n').split('\t')
			chromosome, source, feature_type, start, end, score, strand, frame, attribute = items
			start = int(start); end = int(end); attribute_dict = parse_GFF_attribute(attribute)
			feature_ID = attribute_dict['ID']
			ftype_feature_info_dict[feature_type][feature_ID] = (chromosome, start, end, strand, attribute_dict)
			feature_ftype_dict[feature_ID] = feature_type
		
		if line[0] == '>':
			identifier = line.strip('\n')[1:].split()[0].split('|')[0].strip()
		else:
			return ftype_feature_info_dict, feature_sequence_dict

		seq = ''
		for line in f:
			if line[0] == '>':
				feature_sequence_dict[identifier] = seq
				header = line.strip('\n')[1:]
				identifier = header.split()[0].split('|')[0].strip()
				seq = ''
			else:
				seq += line.strip('\n')
		feature_sequence_dict[identifier] = seq
	
	return ftype_feature_info_dict, feature_sequence_dict

# ====================
# VCF parsing tools
# ====================

def AD_to_depth_and_AAF(AD):
	"""
	Converts AD (allele depth) string to alternate allele frequency
	"""
	depths = np.array([int(d) for d in AD.split(',')])
	total_depth = sum(depths)
	AAF = sum(depths[1:])/total_depth
	return total_depth, AAF

def get_major_alt_allele(AD, ALT):
	"""
	Gets major alt allele from AD (allele depth) string
	"""
	alt_depths = [int(d) for d in AD.split(',')][1:]
	alt_alleles = ALT.split(',')
	allele_depth_tups = sorted(zip(alt_alleles, alt_depths))
	major_alt_allele, major_alt_allele_depth = allele_depth_tups[-1]
	return major_alt_allele

def get_major_allele_num(AD):
	"""
	Gets major alt allele from AD (allele depth) string
	"""
	num_depth_tups = [(num, int(d)) for num, d in enumerate(AD.split(','))]
	num, depth = sorted(num_depth_tups, key=lambda tup: tup[1])[-1]
	return num

def compare_child_to_parent(child_AD, parent_AD):
	"""
	Performs statistical test for difference between two AD (allele depth) strings
	"""
	pADs = [int(d) for d in parent_AD.split(',')]; pADs = [pADs[0], sum(pADs[1:])]
	cADs = [int(d) for d in child_AD.split(',')]; cADs = [cADs[0], sum(cADs[1:])]
	table = np.array([pADs, cADs]).T
	if len(pADs) == 2:
		res = fisher_exact(table, alternative='two-sided')
		return res.pvalue
	else:
		res = chi2_contingency(table)
		return res.pvalue

def get_format_data_dict(FORMAT, data):
	"""
	Converts VCF FORMAT and data column entries to a key-value dictionary
	"""
	return {key: val for key, val in zip(FORMAT.split(':'), data.split(':'))}

def parse_info(INFO):
	"""
	Converts VCF INFO string to key-value dictionary
	"""
	items = INFO.split(';')
	info_dict = {}
	for item in items:
		try:
			key, value = item.split('=')
			info_dict[key] = value
		except:
			info_dict[item] = True
	return info_dict

def filter_variant(QUAL, info_dict):
	"""
	Returns boolean for whether a variant should be filtered out
	"""
	if QUAL == '.':
		return True
	
	if "ReadPosRankSum" in info_dict:
		rp_rs = float(info_dict["ReadPosRankSum"])
		
		if rp_rs > 10 or rp_rs < -10:
			return True
	
	if "QD" in info_dict:
		qd = float(info_dict["QD"])
		
		if qd < 1.5:
			return True
	
	if "MQRankSum" in info_dict:
		mq_rs = float(info_dict["MQRankSum"])
		
		if mq_rs < -14:
			return True
	
	return False

def store_SnpEff_effect(snpeff_effect, allele_anno_dict, force_replace=False):
	"""
	Parses SnpEff variant effect string and updates anno. dictionary entry
	"""
	
	for anno in snpeff_effect.split(','):
		
		last = anno.split('|')[-1]
		if 'INFO_' in last or 'WARNING_' in last:
			alt_allele = anno.split('|')[-2].split(')')[0]
		else:
			alt_allele = anno.split('|')[-1].split(')')[0]
		
		effect, rest = anno.split('(')
		# Note that codon_change is distance for intergenics
		impact, codon_change, aa_change, _, gene = rest.split('|')[1:6]
		gene = unquote(gene)
		
		if (not force_replace) and (alt_allele in allele_anno_dict): # Already has its own annotation
			
			if allele_anno_dict[alt_allele][0] == 'intergenic_region': # Replace intergenic
				allele_anno_dict[alt_allele] = (effect, impact, codon_change, aa_change, gene)
			elif effect == 'intergenic_region':
				pass
			elif 'stream' in allele_anno_dict[alt_allele][0]:
				if 'stream' in effect:
					if int(codon_change) < int(allele_anno_dict[alt_allele][2]): # Replace with closer down/upstream
						allele_anno_dict[alt_allele] = (effect, impact, codon_change, aa_change, gene)
				else: # Replace stream with non-stream
					allele_anno_dict[alt_allele] = (effect, impact, codon_change, aa_change, gene)
			elif 'stream' in effect:
				pass
			elif gene.endswith('.1'):
				allele_anno_dict[alt_allele] = (effect, impact, codon_change, aa_change, gene)
			elif (effect, impact, codon_change, aa_change, gene) != allele_anno_dict[alt_allele]:
				pass
				# print(snpeff_effect); print(allele_anno_dict[alt_allele])
				# print((effect, impact, codon_change, aa_change, gene)); print()
		else:
			allele_anno_dict[alt_allele] = (effect, impact, codon_change, aa_change, gene)

def inspect_VCF_row(filepath, desired_chrom, desired_pos):
	"""
	Pretty prints row in VCF at desired genomic site
	"""
	try:
		f = open(filepath, 'r')
	except:
		print("File not found: %s" % filepath)
		return False
	
	for line in f:
		if line.startswith('#CHROM'):
			header_items = line.strip().split('\t')
			break
	
	cmd = f'grep "{desired_chrom}\t{desired_pos}" {filepath} | head -n 1'
	process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
	line = process.communicate()[0][:-1].decode('ascii')
	for column, datum in zip(header_items, line.split('\t')):
		print(column, ':', datum)


def parse_VCF_simple_keyword(filepath, site_allele_anno_dict, keyword):
	'''
	Reads a VCF, simply returning a dictionary: sample -> locus -> (ref, alt, GT, AD, depth)
	Filter for lines containing certain keyword(s)
	'''
	with open(filepath, 'r') as f:
	
		for line in f:
			if line.startswith('#CHROM'):
				header_items = line.strip().split('\t')
				samples = header_items[9:]
				break
		
		sample_locus_variant_dict = defaultdict(dict)
		
		for line in f:
			if keyword not in line:
				continue
			items = line.strip().split('\t')
			CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = items[:9]
			POS = int(POS); QUAL = float(QUAL); locus = (CHROM, POS)

			info_dict = parse_info(INFO)
			try:
				allele_anno_dict = site_allele_anno_dict[CHROM][POS]
			except:
				site_allele_anno_dict[CHROM] = defaultdict(dict)
				allele_anno_dict = site_allele_anno_dict[CHROM][POS]
			store_SnpEff_effect(info_dict['EFF'], allele_anno_dict)
			
			for sample, data in zip(samples, items[9:]):
				data_dict = get_format_data_dict(FORMAT, data)
				sample_locus_variant_dict[sample][locus] = (REF, ALT, 
															data_dict['GT'], 
															data_dict['AD'] if 'AD' in data_dict else '', 
															data_dict['DP'] if 'DP' in data_dict else '')
		return sample_locus_variant_dict


def parse_VCF_simple(filepath, site_allele_anno_dict, loci_of_interest=None, basic_filter=False):
	'''
	Reads a VCF, simply returning a dictionary: sample -> locus -> (ref, alt, GT, AD, depth)
	No filters applied
	Can specify loci_of_interest
	'''
	with open(filepath, 'r') as f:
	
		for line in f:
			if line.startswith('#CHROM'):
				header_items = line.strip().split('\t')
				samples = header_items[9:]
				break
		
		sample_locus_variant_dict = defaultdict(dict)
		
		for line in f:
			items = line.strip().split('\t')
			CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = items[:9]
			POS = int(POS); QUAL = float(QUAL); locus = (CHROM, POS)

			info_dict = parse_info(INFO)
			
			if basic_filter and filter_variant(QUAL, info_dict):
				continue
			
			try:
				allele_anno_dict = site_allele_anno_dict[CHROM][POS]
			except:
				site_allele_anno_dict[CHROM] = defaultdict(dict)
				allele_anno_dict = site_allele_anno_dict[CHROM][POS]
			store_SnpEff_effect(info_dict['EFF'], allele_anno_dict)

			if (loci_of_interest is None) or (locus in loci_of_interest):

				for sample, data in zip(samples, items[9:]):
					data_dict = get_format_data_dict(FORMAT, data)
					sample_locus_variant_dict[sample][locus] = (REF, ALT, 
																data_dict['GT'], 
																data_dict['AD'] if 'AD' in data_dict else '', 
																data_dict['DP'] if 'DP' in data_dict else '')
		return sample_locus_variant_dict

def parse_VCF_simple_sample(filepath, site_allele_anno_dict, sample):
	'''
	Reads a VCF, returning dictionary for given sample: locus -> (ref, alt, GT, AD, depth)
	No filters applied
	'''
	# Get the VCF header line
	cmd_header = ["grep", "^#CHROM", filepath]
	header_line = subprocess.run(cmd_header, capture_output=True, text=True).stdout.strip()
	columns = header_line.split("\t")

	try:
		sample_idx = columns.index(sample) + 1  # Convert 0-based index to 1-based for cut
	except ValueError:
		raise ValueError(f"Sample '{sample}' not found in VCF header")

	# Use cut to extract the column
	cmd_cut = ["cut", f"-f1,2,4,5,8,9,{sample_idx}", filepath]
	process = subprocess.Popen(cmd_cut, stdout=subprocess.PIPE, text=True)

	# Process the extracted column
	for line in process.stdout:
		if line.startswith('#CHROM'):
			header_items = line.strip().split('\t')
			samples = header_items[9:]
			break
	
	locus_variant_dict = defaultdict(dict)
	
	for line in process.stdout:
		items = line.strip().split('\t')
		CHROM, POS, REF, ALT, INFO, FORMAT, data = items
		
		POS = int(POS); locus = (CHROM, POS)

		info_dict = parse_info(INFO)
		try:
			allele_anno_dict = site_allele_anno_dict[CHROM][POS]
		except:
			site_allele_anno_dict[CHROM] = defaultdict(dict)
			allele_anno_dict = site_allele_anno_dict[CHROM][POS]
		store_SnpEff_effect(info_dict['EFF'], allele_anno_dict)
		
		data_dict = get_format_data_dict(FORMAT, data)
		locus_variant_dict[locus] = (REF, ALT, data_dict['GT'], 
							   data_dict['AD'] if 'AD' in data_dict else '',
							   data_dict['DP'] if 'DP' in data_dict else '')
	
	return locus_variant_dict

def update_annotations_from_annotated_VCF(filepath, site_allele_anno_dict, force_replace=True):
	"""
	Reads a SnpEff-annotated VCF, updating variant annotation dictionary only
	"""
	try:
		f = open(filepath, 'r')
	except:
		print("File not found: %s" % filepath)
		return False
	
	for line in f:
		if line.startswith('#CHROM'):
			header_items = line.strip().split('\t')
			break
	
	for line in f:
		items = line.strip().split('\t')
		CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = items[:9]
		if QUAL == '.':
			continue
		POS = int(POS); QUAL = float(QUAL)
		info_dict = parse_info(INFO)
		
		# Initial quality filter
		if filter_variant(QUAL, info_dict):
			continue
		
		# Store annotations
		try:
			allele_anno_dict = site_allele_anno_dict[CHROM][POS]
			store_SnpEff_effect(info_dict['EFF'], allele_anno_dict, force_replace=force_replace)
		except:
			pass
	
	return


def parse_annotated_VCF(filepath, site_allele_anno_dict, parent_sample, desired_samples=None, p_threshold=0.001, tqdm_on=False):
	"""
	Reads a SnpEff-annotated VCF, updating variant annotation dictionary and returning sample variant records
	
	If parent_sample is specified, only variants in a sample that differ from the parent sample with Fisher's exact test P-value < p_threshold are retained
	"""
	
	clone_mutation_records_dict = defaultdict(list)
	
	try:
		f = open(filepath, 'r')
	except:
		print("File not found: %s" % filepath)
		return False
	
	for line in f:
		if line.startswith('#CHROM'):
			header_items = line.strip().split('\t')
			break
	
	desired_sample_idxs = []
	for idx in range(9, len(header_items)):
		sample = header_items[idx]
		if desired_samples is None or sample in desired_samples:
			desired_sample_idxs.append(idx)
	
	parent_sample_idx = header_items.index(parent_sample)
	
	tqdm_switch_func = tqdm if tqdm_on else (lambda x: x)
	
	num_records_retained = 0
	for line in tqdm_switch_func(f):
		
		items = line.strip().split('\t')
		CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = items[:9]
		POS = int(POS); QUAL = float(QUAL)
		info_dict = parse_info(INFO)
		
		# Initial quality filter
		if filter_variant(QUAL, info_dict):
			continue
		
		# Store annotations
		try:
			allele_anno_dict = site_allele_anno_dict[CHROM][POS]
		except:
			site_allele_anno_dict[CHROM] = defaultdict(dict)
			allele_anno_dict = site_allele_anno_dict[CHROM][POS]
		store_SnpEff_effect(info_dict['EFF'], allele_anno_dict)
		
		# Store child sample data
		for idx in desired_sample_idxs:
			sample = header_items[idx]; data = items[idx]
			data_dict = get_format_data_dict(FORMAT, data)
			
			if 'DP' not in data_dict or data_dict['DP'] in ['.', '0']:
				continue
			
			GT = data_dict['GT']
			AD = data_dict['AD']
			
			parent_data = items[parent_sample_idx]
			parent_data_dict = get_format_data_dict(FORMAT, parent_data)
			
			if 'DP' not in parent_data_dict or parent_data_dict['DP'] in ['.', '0', '1']:
				continue
			
			pvalue = compare_child_to_parent(AD, parent_data_dict['AD'])
			if pvalue < p_threshold:
				record = (CHROM, POS, QUAL, GT, AD, parent_data_dict['AD'], REF, ALT)
				clone_mutation_records_dict[sample].append(record)
				num_records_retained += 1
	
	print("Number of variants retained: %i" % num_records_retained)
	
	return clone_mutation_records_dict

def annotate_record(mutation_record, site_allele_anno_dict):
	CHROM, POS, QUAL, GT, AD, parent_AD, REF, ALT = mutation_record
	major_alt_allele = get_major_alt_allele(AD, ALT)
	effect, impact, codon_change, aa_change, gene = site_allele_anno_dict[CHROM][POS][major_alt_allele]
	return effect, impact, codon_change, aa_change, gene

def filter_mutation_records(clone_mutation_records_dict, min_quality=None, min_depth=None, min_AAF=None, min_AAF_diff=None, max_AAF_diff=None, max_pvalue=None):
	"""
	Using the clone_mutation_records_dict object from parsing a VCF, filter for variants based on quality, depth, allele frequency, and/or allele frequency difference
	"""
	filtered_clone_mutation_records_dict = defaultdict(list)
	
	for sample in clone_mutation_records_dict:
		for record in clone_mutation_records_dict[sample]:
			
			CHROM, POS, QUAL, GT, AD, parent_AD, REF, ALT = record
			depth, AAF = AD_to_depth_and_AAF(AD)
			parent_depth, parent_AAF = AD_to_depth_and_AAF(parent_AD)
			AAF_diff = np.abs(parent_AAF - AAF)
			
			if min_quality and QUAL < min_quality:
				continue
			if min_depth and depth < min_depth:
				continue
			if min_AAF and AAF < min_AAF:
				continue
			if min_AAF_diff and AAF_diff < min_AAF_diff:
				continue
			if max_pvalue:
				pvalue = compare_child_to_parent(AD, parent_AD)
				if pvalue > max_pvalue:
					continue
			
			filtered_clone_mutation_records_dict[sample].append(record)
	
	return filtered_clone_mutation_records_dict

def summarize_mutation_statistics(clone_mutation_records_dict, samples):
	"""
	Summarizes quality and allele depth statistics for a set of variants/mutations
	"""
	
	fig, ax = plt.subplots(3, 1)
	
	depths = []; AAFs = []; AAF_diffs = []
	for sample in samples:
		for record in clone_mutation_records_dict[sample]:
			
			CHROM, POS, QUAL, GT, AD, parent_AD, REF, ALT = record
			depth, AAF = AD_to_depth_and_AAF(AD)
			parent_depth, parent_AAF = AD_to_depth_and_AAF(parent_AD)
			AAF_diff = parent_AAF - AAF
			
			depths.append(depth); AAFs.append(AAF); AAF_diffs.append(AAF_diff)
	
	ax[0].hist(depths, bins=30); ax[0].set_xlabel("Depth")
	ax[1].hist(AAFs, bins=30); ax[1].set_xlabel("AAF")
	ax[2].hist(AAF_diffs, bins=30); ax[2].set_xlabel("AAF diff")
	plt.tight_layout()
	plt.show()
