# =====================================================================
# Run this to generate qsub scripts for WGS alignment, variant calling
# Make sure to update config.cfg (in same directory)
# and run `python rename_fastq_generate_samples.py` first
# =====================================================================

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
SAMPLES_PATH = get_env_value_from_config('samples_file')
FASTQ_DIR = get_env_value_from_config('fastq_dir')
SPECIES_ABBR = get_env_value_from_config('species_abbr')
REF_FASTA_PATH = f"{get_env_value_from_config('ref_dir')}/{get_env_value_from_config('ref_fasta').split('/')[-1]}"

# Fastq file name parameters
FASTQ_PREFIX = SPECIES_ABBR + '_'
FASTQ_SUFFIX = '_R'
FASTQ_GZIP_SUFFIX = '.gz'

samples = [line.strip() for line in open(SAMPLES_PATH, 'r') if line.strip() != '']
num_samples = len(samples)

fastqgz_sizes = []
for sample in samples:
	for readnum in [1,2]:
		fastq_path = f"{FASTQ_DIR}/{FASTQ_PREFIX}{sample}{FASTQ_SUFFIX}{readnum}.fastq{FASTQ_GZIP_SUFFIX}"
		try:
			fastqgz_sizes.append(os.path.getsize(fastq_path))
		except:
			print("Couldn't find %s, assuming ~5 GB" % fastq_path)
			fastqgz_sizes.append(5_000_000_000)

total_fastqgz_size = sum(fastqgz_sizes)
biggest_fastqgz_size = max(fastqgz_sizes)

# Reference for this calculation:
# qsub_1: max 9 GB (one read pair) took 18 hours
# qsub_2: total 36 GB, 10 samples took around 160 hours, 63 GB, 4 samples took 92 hours

qsub_1_walltime_hours = int((biggest_fastqgz_size/1_000_000_000.0)*3) + 20
qsub_2_total_walltime_hours = int((total_fastqgz_size/1_000_000_000.0)*5) + 40
qsub_2_walltime_hours = 88 # Should be big to be safe

# Number of separate qsub_2 scripts

num_qsub_2_scripts = (qsub_2_total_walltime_hours//qsub_2_walltime_hours) + 8

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

print(sorted(chrom_length_dict.items()))

genome_length = sum(chrom_length_dict.values())
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

gatk_input_combined_samples_str = '\n'.join(["\t-I $main_dir/%s.ready.bam \\" % sample for sample in samples])

for idx, chromosome_interval_list in enumerate(chromosome_interval_lists):
	gatk_specify_chrom_options = '' if num_qsub_2_scripts == 1 else ('\n\t-L ' + ' -L '.join(chromosome_interval_list) + ' \\')
	part_str = '' if num_qsub_2_scripts == 1 else "-part%i" % (idx + 1)
	
	qsub_2_call_variants_str = \
f'''#!/bin/sh
#SBATCH --job-name {JOB_NAME}{part_str}_callvar
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
	-O $main_dir/{GROUP_NAME}{part_str}.raw.snps.indels.vcf

echo "Done!"
'''
	
	f = open(f"sbatch_2_call_variants{part_str}", 'w')
	f.write(qsub_2_call_variants_str)
	f.close()

qsub_0_bwa_str = \
f'''#!/bin/sh
#SBATCH --job-name {JOB_NAME}_bwa-mem
#SBATCH --partition=hotel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time={qsub_1_walltime_hours}:00:00 
#SBATCH --qos=hotel
#SBATCH --account={ACCOUNT}
#SBATCH --export=ALL
#SBATCH --output {LOG_DIR}/slurm-%A_%a.%x.out-%N
#SBATCH --output {LOG_DIR}/slurm-%A_%a.%x.err-%N
#SBATCH --mail-type {EMAIL_OPTIONS}
#SBATCH --mail-user {EMAIL}
#SBATCH --array=1-{num_samples}

. {CONFIG_PATH}
'''

qsub_0_bwa_str += '''
module load shared
module load samtools
module load gatk
module load picard

readarray samples < $samples_file
samples=(null ${samples[@]}) # zero to one start index
sample=${samples[$SLURM_ARRAY_TASK_ID]}
echo $sample

# ================================================================================

echo "Running bwa mem..."

/tscc/projects/ps-winzelerlab/TOOLS/bin/bwa mem -M -t 16 $ref_fasta $fastq_dir/%s${sample}%s1.fastq%s $fastq_dir/%s${sample}%s2.fastq%s > $main_dir/${sample}.sam

echo "Done!"
''' % (FASTQ_PREFIX, FASTQ_SUFFIX, FASTQ_GZIP_SUFFIX, FASTQ_PREFIX, FASTQ_SUFFIX, FASTQ_GZIP_SUFFIX)

qsub_1_align_str = \
f'''#!/bin/sh
#SBATCH --job-name {JOB_NAME}_align
#SBATCH --partition=hotel
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=16GB
#SBATCH --ntasks-per-node=1
#SBATCH --time={qsub_1_walltime_hours}:00:00 
#SBATCH --qos=hotel
#SBATCH --account={ACCOUNT}
#SBATCH --export=ALL
#SBATCH --exclude tscc-11-2
#SBATCH --output {LOG_DIR}/slurm-%A_%a.%x.out-%N
#SBATCH --output {LOG_DIR}/slurm-%A_%a.%x.err-%N
#SBATCH --mail-type {EMAIL_OPTIONS}
#SBATCH --mail-user {EMAIL}
#SBATCH --array=1-{num_samples}

. {CONFIG_PATH}
'''

qsub_1_align_str += '''
module load shared
module load samtools
module load gatk
module load picard

readarray samples < $samples_file
samples=(null ${samples[@]}) # zero to one start index
sample=${samples[$SLURM_ARRAY_TASK_ID]}
echo $sample

# ================================================================================

# echo "Running bwa mem..."

# /tscc/projects/ps-winzelerlab/TOOLS/bin/bwa mem -M -t 16 $ref_fasta $fastq_dir/${species_abbr}_${sample}_R1.fastq.gz $fastq_dir/${species_abbr}_${sample}_R2.fastq.gz > $main_dir/${sample}.sam

echo "Converting sam to bam"

samtools view -b -S -o $main_dir/${sample}.bam $main_dir/${sample}.sam

echo "Sorting bam"

samtools sort -O bam -o $main_dir/${sample}_sorted.bam $main_dir/${sample}.bam

echo "Create index"

samtools index $main_dir/${sample}_sorted.bam

echo "Running Picard AddOrReplaceReadGroups..."

picard AddOrReplaceReadGroups \\
I=$main_dir/${sample}_sorted.bam \\
O=$main_dir/${sample}_sorted.rg.bam \\
RGID=4 \\
RGLB=lib1 \\
RGPL=illumina \\
RGPU=unit1 \\
RGSM=$sample \\
CREATE_INDEX=true

echo "Running Picard MarkDuplicates..."

picard MarkDuplicates \\
I=$main_dir/${sample}_sorted.rg.bam \\
O=$main_dir/${sample}_sorted.rg.md.bam \\
M=$main_dir/${sample}_sorted.rg.md.metrics.txt \\
CREATE_INDEX=true

# Skip BaseRecalibrator

echo "Running GATK PrintReads..."

gatk PrintReads \\
        -R $ref_fasta \\
        -I $main_dir/${sample}_sorted.rg.md.bam \\
        --read-filter MappingQualityNotZeroReadFilter \\
        -O $main_dir/${sample}.ready.bam

echo "Done!"
'''

f = open("sbatch_0_bwa", 'w')
f.write(qsub_0_bwa_str)
f.close()

f = open("sbatch_1_align", 'w')
f.write(qsub_1_align_str)
f.close()

qsub_cnv_analysis_str = \
f'''#!/bin/sh
#SBATCH --job-name {JOB_NAME}_CNV_analysis
#SBATCH --partition=hotel
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks-per-node=1
#SBATCH --time=08:00:00
#SBATCH --qos=hotel
#SBATCH --account={ACCOUNT}
#SBATCH --export=ALL
#SBATCH --output {LOG_DIR}/slurm-%j.out-%N
#SBATCH --output {LOG_DIR}/slurm-%j.err-%N
#SBATCH --mail-type {EMAIL_OPTIONS}
#SBATCH --mail-user {EMAIL}
#SBATCH --array=1-{num_samples}

. {CONFIG_PATH}

module load shared
module load gatk

gatk_cnv_dir=/tscc/projects/ps-winzelerlab/GENOME_RESOURCES/p_fal/GATK_CNV
pon_path=$gatk_cnv_dir/PON/{STRAIN}cnv.pon.hdf5
interval_list_path=$gatk_cnv_dir/references/{SPECIES_ABBR}.preprocessed.interval_list
'''

qsub_cnv_analysis_str += \
'''
readarray samples < $samples_file
samples=(null ${samples[@]}) # zero to one start index
sample=${samples[$SLURM_ARRAY_TASK_ID]}
echo $sample

echo "Working on sample $sample"

bam_path=$main_dir/${sample}.ready.bam

echo "Collecting read counts..."
gatk CollectReadCounts \
        -I $bam_path \
        -L $interval_list_path \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O $main_dir/${sample}.counts.hdf5

echo "Denoising read counts..."
gatk DenoiseReadCounts \
        -I $main_dir/${sample}.counts.hdf5 \
        --count-panel-of-normals $pon_path \
        --standardized-copy-ratios $main_dir/${sample}.standardizedCR.tsv \
        --denoised-copy-ratios $main_dir/${sample}.denoisedCR.tsv

echo "Done"
'''

f = open("sbatch_cnv_analysis", 'w')
f.write(qsub_cnv_analysis_str)
f.close()
