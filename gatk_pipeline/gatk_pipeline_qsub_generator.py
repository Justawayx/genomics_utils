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

LOG_DIR = "/home/dwc001/scratch/logs"
EMAIL = "dwc001@ucsd.edu"
EMAIL_OPTIONS = "ea"	

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

SAMPLES_PATH = get_env_value_from_config('samples_file')
FASTQ_DIR = get_env_value_from_config('fastq_dir')
SPECIES_ABBR = get_env_value_from_config('species_abbr')
REF_FASTA_PATH = f"{get_env_value_from_config('ref_dir')}/{get_env_value_from_config('ref_fasta').split('/')[-1]}"

samples = [line.strip() for line in open(SAMPLES_PATH, 'r') if line.strip() != '']
num_samples = len(samples)

fastqgz_sizes = []
for sample in samples:
	for readnum in [1,2]:
		fastq_path = f"{FASTQ_DIR}/{SPECIES_ABBR}_{sample}_R{readnum}.fastq.gz"
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
qsub_2_total_walltime_hours = int((total_fastqgz_size/1_000_000_000.0)*5) + 10
qsub_2_walltime_hours = 150 # Should be big to be safe

# Number of separate qsub_2 scripts

num_qsub_2_scripts = (qsub_2_total_walltime_hours//qsub_2_walltime_hours) + 1

chrom_length_dict = {}
chromosome_boundaries_ordered = [0]
chromosomes_ordered = []
with open(REF_FASTA_PATH, 'r') as f:
	for line in f:
		if line.startswith('>'):
			items = line[1:].split('|')
			chrom = items[0].strip()
			for item in items:
				if 'length=' in item:
					length = int(item.split('length=')[1].strip())
			chrom_length_dict[chrom] = length
			chromosome_boundaries_ordered.append(chromosome_boundaries_ordered[-1] + length)
			chromosomes_ordered.append(chrom)

genome_length = sum(chrom_length_dict.values())
genome_breakpoints = [((genome_length/num_qsub_2_scripts)*i) for i in range(1,num_qsub_2_scripts)] + [genome_length]
chromosome_lists = [] # List of lists for each qsub_2 script

cur_breakpoint_idx = 0
cur_chromosome_list = []
for i in range(len(chromosomes_ordered)):
	chrom = chromosomes_ordered[i]
	end_pos = chromosome_boundaries_ordered[i+1]
	cur_chromosome_list.append(chrom)
	if end_pos >= genome_breakpoints[cur_breakpoint_idx]:
		cur_breakpoint_idx += 1
		chromosome_lists.append(cur_chromosome_list)
		cur_chromosome_list = []

gatk_input_combined_samples_str = '\n'.join(["\t-I $main_dir/%s.ready.bam \\" % sample for sample in samples])

for idx, chromosome_list in enumerate(chromosome_lists):
	gatk_specify_chrom_options = '' if num_qsub_2_scripts == 1 else ('\n\t-L ' + ' -L '.join(chromosome_list) + ' \\')
	part_str = '' if num_qsub_2_scripts == 1 else "-part%i" % (idx + 1)
	
	qsub_2_call_variants_str = \
f'''#!/bin/sh
#PBS -A winzeler-group
#PBS -N {JOB_NAME}_callvar
#PBS -l nodes=1:ppn=16,walltime={qsub_2_walltime_hours}:00:00
#PBS -o {LOG_DIR}
#PBS -e {LOG_DIR}
#PBS -M {EMAIL}
#PBS -m {EMAIL_OPTIONS}

. {CONFIG_PATH}

echo "Running GATK HaplotypeCaller..."
java -jar $gatk_dir/GenomeAnalysisTK.jar \\
	-T HaplotypeCaller \\
	-R $ref_fasta \\
{gatk_input_combined_samples_str}{gatk_specify_chrom_options}
	-o $main_dir/{GROUP_NAME}{part_str}.raw.snps.indels.vcf

echo "Done!"
'''
	
	f = open(f"qsub_2_call_variants{part_str}", 'w')
	f.write(qsub_2_call_variants_str)
	f.close()

qsub_1_align_str = \
f'''#!/bin/sh
#PBS -A winzeler-group
#PBS -N {JOB_NAME}_align
#PBS -l nodes=1:ppn=16,walltime={qsub_1_walltime_hours}:00:00
#PBS -o {LOG_DIR}
#PBS -e {LOG_DIR}
#PBS -M {EMAIL}
#PBS -m {EMAIL_OPTIONS}
#PBS -t 1-{num_samples}

. {CONFIG_PATH}
'''

qsub_1_align_str += '''
module load bwa

readarray samples < $samples_file
samples=(null ${samples[@]}) # zero to one start index
sample=${samples[$PBS_ARRAYID]}
echo $sample

# ================================================================================

echo "Running bwa mem..."

bwa mem -M -t 16 $ref_fasta $fastq_dir/${species_abbr}_${sample}_R1.fastq.gz $fastq_dir/${species_abbr}_${sample}_R2.fastq.gz > $main_dir/${sample}.sam

echo "Converting sam to bam"

samtools view -b -S -o $main_dir/${sample}.bam $main_dir/${sample}.sam

echo "Sorting bam"

samtools sort -O bam -o $main_dir/${sample}_sorted.bam $main_dir/${sample}.bam

echo "Create index"

samtools index $main_dir/${sample}_sorted.bam

echo "Running Picard AddOrReplaceReadGroups..."

java -jar /opt/biotools/picard/picard.jar AddOrReplaceReadGroups \\
I=$main_dir/${sample}_sorted.bam \\
O=$main_dir/${sample}_sorted.rg.bam \\
RGID=4 \\
RGLB=lib1 \\
RGPL=illumina \\
RGPU=unit1 \\
RGSM=$sample \\
CREATE_INDEX=true

echo "Running Picard MarkDuplicates..."

java -Xmx6g -jar /opt/biotools/picard/picard.jar MarkDuplicates \\
I=$main_dir/${sample}_sorted.rg.bam \\
O=$main_dir/${sample}_sorted.rg.md.bam \\
M=$main_dir/${sample}_sorted.rg.md.metrics.txt \\
CREATE_INDEX=true

echo "Running GATK RelaignerTargetCreator..."

java -Xmx12g -jar $gatk_dir/GenomeAnalysisTK.jar \\
        -T RealignerTargetCreator \\
        -R $ref_fasta \\
        -I $main_dir/${sample}_sorted.rg.md.bam \\
        -o $main_dir/${sample}_sorted.rg.md.intervals

echo "Running GATK IndelRealigner..."

java -Xmx12g -jar $gatk_dir/GenomeAnalysisTK.jar \\
        -T IndelRealigner \\
        -R $ref_fasta \\
        -I $main_dir/${sample}_sorted.rg.md.bam \\
        -targetIntervals $main_dir/${sample}_sorted.rg.md.intervals \\
        -o $main_dir/${sample}_sorted.rg.md.ir.bam

# Skipping BaseRecalibrator

echo "Running GATK PrintReads..."

java -Xmx12g -jar $gatk_dir/GenomeAnalysisTK.jar \\
        -T PrintReads \\
        -R $ref_fasta \\
        -I $main_dir/${sample}_sorted.rg.md.ir.bam \\
        --read_filter MappingQualityZero \\
        -o $main_dir/${sample}.ready.bam

echo "Done!"
'''

f = open("qsub_1_align", 'w')
f.write(qsub_1_align_str)
f.close()