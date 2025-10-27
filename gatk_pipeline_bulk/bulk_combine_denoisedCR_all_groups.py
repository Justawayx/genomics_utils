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
sample_groups_file = get_env_value_from_config('sample_groups_file')
output_dir = get_env_value_from_config('main_dir')

group_samples_dict = defaultdict(list)

with open(sample_groups_file, 'r') as f:
    for line in f:
        group, sample = line.strip('\n').split('\t')
        group_samples_dict[group].append(sample)

# ====================================

strain_interval_gene_id_desc_class_dict = {}
strain_groups_dict = defaultdict(list)

for group in group_samples_dict:
    
    if '3D7' in group:
        strain = '3D7'
    elif 'Dd2' in group:
        strain = 'Dd2'
    else:
        strain = input("Enter strain (3D7 or Dd2) for CNV analysis of %s:" % group)

    strain_groups_dict[strain].append(group)

for strain in strain_groups_dict:

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

    strain_interval_gene_id_desc_class_dict[strain] = interval_gene_id_desc_class_dict

for strain in strain_groups_dict:
    
    interval_sample_cr_dict = defaultdict(dict) # (contig, start, end) -> sample -> CR
    ordered_samples = []
    ordered_intervals = [] # List of (contig, start, end) tuples
    
    for group in group_samples_dict:
        ordered_samples += group_samples_dict[group]
    
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
    
    # Combine into output file
    
    header_items = ['CONTIG', 'START', 'END', 'Gene_ID', 'Gene_Description', 'Miles-et-al-2016-Genome-Classification'] + ordered_samples
    
    output_filepath = '%s/%s_%s_CNV_analysis.tsv' % (output_dir, group_name, strain)
    of = open(output_filepath, 'w')
    
    of.write('\t'.join(header_items) + '\n')
    
    for interval in ordered_intervals:
        contig, start, end = interval
        gene_id, gene_desc, miles_genome_class = strain_interval_gene_id_desc_class_dict[strain][interval]
        items = [contig, start, end, gene_id, gene_desc, miles_genome_class]
        for sample in ordered_samples:
            log2_cr = interval_sample_cr_dict[interval][sample]
            items.append(log2_cr)
        
        of.write('\t'.join([str(item) for item in items]) + '\n')
    
    of.close()

    print("Wrote %s" % output_filepath)

