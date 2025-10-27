# =====================================================================
# Run this to generate qsub scripts for WGS alignment, variant calling
# Make sure to update config.cfg (in same directory)
# and run `python rename_fastq_generate_samples.py` first
# =====================================================================

from collections import defaultdict
import os, sys

# ======================================
# TSCC parameters
# ======================================

# Change these

ACCOUNT = "htl124"
LOG_DIR = "/tscc/lustre/ddn/scratch/dwc001/logs"
EMAIL = "dwc001@ucsd.edu"
EMAIL_OPTIONS = "END"

# ======================================
# More parameters (change if necessary)
# ======================================

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

GROUP_NAME = get_env_value_from_config('group_name')
JOB_NAME = GROUP_NAME

STRAIN = get_env_value_from_config('strain')
SAMPLE_GROUPS_PATH = get_env_value_from_config('sample_groups_file')
SAMPLES_PATH = get_env_value_from_config('samples_file')
FASTQ_DIR = get_env_value_from_config('fastq_dir')
SPECIES_ABBR = get_env_value_from_config('species_abbr')
REF_FASTA_PATH = f"{get_env_value_from_config('ref_dir')}/{get_env_value_from_config('ref_fasta').split('/')[-1]}"

# Load genome information
chrom_length_dict = {}
chromosome_boundaries_ordered = [0]
chromosomes_ordered = []
chrom = 'dummy'
with open(REF_FASTA_PATH, 'r') as f:
    for line in f:
        if line.startswith('>'):
            if chrom != 'dummy':
                chrom_length_dict[chrom] = seq_len
                chromosome_boundaries_ordered.append(chromosome_boundaries_ordered[-1] + seq_len)
                chromosomes_ordered.append(chrom)
            seq_len = 0
            chrom = line[1:].split('|')[0].split(' ')[0]
        else:
            seq_len += len(line.strip('\n'))

chrom_length_dict[chrom] = seq_len
chromosome_boundaries_ordered.append(chromosome_boundaries_ordered[-1] + seq_len)
chromosomes_ordered.append(chrom)
genome_length = sum(chrom_length_dict.values())
print(sorted(chrom_length_dict.items()))

# Fastq file name parameters
FASTQ_PREFIX = SPECIES_ABBR + '_'
FASTQ_SUFFIX = '_R'
FASTQ_GZIP_SUFFIX = '.gz'

samples = [line.strip() for line in open(SAMPLES_PATH, 'r') if line.strip() != '']
num_samples = len(samples)

group_samples_dict = defaultdict(list)

with open(SAMPLE_GROUPS_PATH, 'r') as f:
    for line in f:
        group, sample = line.strip('\n').split('\t')
        group_samples_dict[group].append(sample)

for group in group_samples_dict:
    scripts_dir = f"scripts_{group}"
    try:
        os.mkdir(scripts_dir)
    except:
        print("Folder already exists, skipping group %s" % group)
    
    group_samples = group_samples_dict[group]

    fastqgz_sizes = []
    for sample in group_samples:
        for readnum in [1,2]:
            fastq_path = f"{FASTQ_DIR}/{FASTQ_PREFIX}{sample}{FASTQ_SUFFIX}{readnum}.fastq{FASTQ_GZIP_SUFFIX}"
            try:
                fastqgz_sizes.append(os.path.getsize(fastq_path))
            except:
                print("Couldn't find %s, assuming ~5 GB" % fastq_path)
                fastqgz_sizes.append(5_000_000_000)
    
    total_fastqgz_size = sum(fastqgz_sizes)
    biggest_fastqgz_size = max(fastqgz_sizes)
    
    qsub_2_total_walltime_hours = int((total_fastqgz_size/1_000_000_000.0)*5) + 10
    qsub_2_walltime_hours = 60 # Should be big to be safe
    
    # Number of separate qsub_2 scripts
    num_qsub_2_scripts = (qsub_2_total_walltime_hours//qsub_2_walltime_hours) + 2
    genome_breakpoints = [0] + [int((genome_length/num_qsub_2_scripts)*i) for i in range(1,num_qsub_2_scripts)] + [genome_length]
    chromosome_interval_lists = [] # List of lists for each qsub_2 script
    cur_chromosome_interval_list = []
    
    for idx in range(len(genome_breakpoints)-1):
        interval_start = genome_breakpoints[idx]
        interval_end = genome_breakpoints[idx+1]
        
        for i in range(len(chromosomes_ordered)):
            chrom = chromosomes_ordered[i]
            chrom_start = chromosome_boundaries_ordered[i]
            chrom_end = chromosome_boundaries_ordered[i+1]
            
            interval_start_internal = interval_start - chrom_start + 1 # May not be valid
            interval_end_internal = interval_end - chrom_start # May not be valid
            
            chrom_start_internal = 1
            chrom_end_internal = chrom_length_dict[chrom]
            
            if interval_start >= chrom_start and interval_start <= chrom_end: # Chromosome straddles left side of interval
                cur_chromosome_interval_list.append('%s:%i-%i' % (chrom, interval_start_internal, chrom_end_internal))
            elif interval_start <= chrom_start and interval_end >= chrom_end: # Chromosome is inside interval
                cur_chromosome_interval_list.append('%s:%i-%i' % (chrom, chrom_start_internal, chrom_end_internal))
            elif interval_end >= chrom_start and interval_end <= chrom_end: # Chromosome straddles right side of interval
                cur_chromosome_interval_list.append('%s:%i-%i' % (chrom, chrom_start_internal, interval_end_internal))
            elif interval_start >= chrom_start and interval_end <= chrom_end: # Chromosome contains entire interval
                cur_chromosome_interval_list.append('%s:%i-%i' % (chrom, interval_start_internal, interval_end_internal))
            else:
                pass
        
        chromosome_interval_lists.append(cur_chromosome_interval_list)
        cur_chromosome_interval_list = []
    
    gatk_input_combined_samples_str = '\n'.join(["\t-I $main_dir/%s.ready.bam \\" % sample for sample in group_samples])
    
    for idx, chromosome_interval_list in enumerate(chromosome_interval_lists):
        gatk_specify_chrom_options = '' if num_qsub_2_scripts == 1 else ('\n\t-L ' + ' -L '.join(chromosome_interval_list) + ' \\')
        part_str = '' if num_qsub_2_scripts == 1 else "-part%i" % (idx + 1)
        
        qsub_2_call_variants_str = \
f'''#!/bin/sh
#SBATCH --job-name {group}{part_str}_callvar
#SBATCH --partition=hotel
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=16G
#SBATCH --ntasks-per-node=1
#SBATCH --time={qsub_2_walltime_hours}:00:00 
#SBATCH --qos=hotel
#SBATCH --account={ACCOUNT}
#SBATCH --export=ALL
#SBATCH --output {LOG_DIR}/slurm-%j.%x.out-%N
#SBATCH --output {LOG_DIR}/slurm-%j.%x.err-%N
#SBATCH --mail-type {EMAIL_OPTIONS}
#SBATCH --mail-user {EMAIL}

. {CONFIG_PATH}

module load shared
module load gatk

echo "Running GATK HaplotypeCaller..."
gatk HaplotypeCaller \\
    -R $ref_fasta \\
{gatk_input_combined_samples_str}{gatk_specify_chrom_options}
    -O $main_dir/{group}{part_str}.raw.snps.indels.vcf

echo "Done!"
'''
        f = open(f"{scripts_dir}/sbatch_2_call_variants{part_str}", 'w')
        f.write(qsub_2_call_variants_str)
        f.close()
