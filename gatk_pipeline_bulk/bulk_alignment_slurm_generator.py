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

qsub_1_walltime_hours = int((biggest_fastqgz_size/1_000_000_000.0)*3) + 20

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
