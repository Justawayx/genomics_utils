from collections import defaultdict
import os, sys

# ====================================
# PARAMETERS
# ====================================
CONFIG_PATH = os.path.abspath('config.cfg')

def get_env_value_from_config(env_variable):
	for line in open(CONFIG_PATH, 'r'):
		if line.startswith('#'):
			continue
		if '=' in line:
			key, value = line.rstrip('\n').split('=')
			if key == env_variable:
				return value
	sys.exit("Environmental variable %s not found in config.cfg" % env_variable)

group_name = get_env_value_from_config('group_name')
strain = get_env_value_from_config('strain')
samples_file = get_env_value_from_config('samples_file')
output_dir = get_env_value_from_config('main_dir')

ordered_samples = [line.strip() for line in open(samples_file, 'r')]
# ====================================

ordered_intervals = [] # List of (contig, start, end) tuples
interval_sample_cr_dict = defaultdict(dict) # (contig, start, end) -> sample -> CR

# Store sample data

FIRST_SAMPLE = True
new_ordered_samples = []

for sample in ordered_samples:
	try:
		f = open('%s/%s.denoisedCR.tsv' % (output_dir, sample), 'r')
		new_ordered_samples.append(sample)
	except:
		print("denoisedCR.tsv file is missing for %s" % sample)
		continue
	
	DATA_STARTED = False
	for line in f:
		
		if DATA_STARTED is True:
			contig, start, end, log2_copy_ratio = line.strip().split('\t')
			start = int(start); end = int(end); log2_cr = float(log2_copy_ratio)
			interval_sample_cr_dict[(contig, start, end)][sample] = log2_cr
			if FIRST_SAMPLE is True:
				ordered_intervals.append((contig, start, end))
		
		if line[0] == '@':
			continue
		
		if line.split('\t')[0] == 'CONTIG':
			DATA_STARTED = True
	
	FIRST_SAMPLE = False

ordered_samples = new_ordered_samples

# Gene annotations

interval_gene_id_desc_class_dict = {} # (contig, start, end) -> (ID, desc, class)

f = open('/tscc/projects/ps-winzelerlab/PROJECTS/daisy/REF_DATA/p_fal/%s_interval_annotations.tsv' % strain, 'r')
header = f.readline()
for line in f:
	items = line.split('\t')
	contig, start, end, = items[:3]
	gene_id, gene_desc, miles_genome_class = items[-3:]
	start = int(start); end = int(end); miles_genome_class = miles_genome_class.strip()
	interval_gene_id_desc_class_dict[(contig, start, end)] = (gene_id, gene_desc, miles_genome_class)

f.close()

# Combine into output file

header_items = ['CONTIG', 'START', 'END', 'Gene_ID', 'Gene_Description', 'Miles-et-al-2016-Genome-Classification'] + ordered_samples

of = open('%s/%s_CNV_analysis.tsv' % (output_dir, group_name), 'w')

of.write('\t'.join(header_items) + '\n')

for interval in ordered_intervals:
	contig, start, end = interval
	gene_id, gene_desc, miles_genome_class = interval_gene_id_desc_class_dict[interval]
	items = [contig, start, end, gene_id, gene_desc, miles_genome_class]
	for sample in ordered_samples:
		log2_cr = interval_sample_cr_dict[interval][sample]
		items.append(log2_cr)
	
	of.write('\t'.join([str(item) for item in items]) + '\n')

of.close()

