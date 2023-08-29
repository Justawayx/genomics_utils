# =====================================================================
# Run this to generate qsub scripts for WGS alignment, variant calling
# Make sure to update config.cfg (in same directory)
# and run `python rename_fastq_generate_samples.py` first
# =====================================================================

import os

# ======================================
# TSCC parameters (CHANGE)
# ======================================

JOB_NAME = "Sarah_Yliu4"
RESOURCES_STR = "nodes=1:ppn=16,walltime=100:00:00"
LOG_DIR = "/home/dwc001/scratch/logs"
EMAIL = "dwc001@ucsd.edu"
EMAIL_OPTIONS = "ea"
GROUP_NAME = "Yliu3-combined"

# ======================================
# More parameters (change if necessary)
# ======================================

CONFIG_PATH = os.path.abspath('config.cfg')

for line in open(CONFIG_PATH, 'r'):
	if line.startswith("samples_file="):
		SAMPLES_PATH = line.split("samples_file=")[1].strip()
		break

samples = [line.strip() for line in open(SAMPLES_PATH, 'r') if line.strip() != '']
num_samples = len(samples)

gatk_input_combined_samples_str = '\n'.join(["\t-I $main_dir/%s.ready.bam \\" % sample for sample in samples])

qsub_2_call_variants_str = \
f'''#!/bin/sh
#PBS -A winzeler-group
#PBS -N {JOB_NAME}_callvar
#PBS -l {RESOURCES_STR}
#PBS -o {LOG_DIR}
#PBS -e {LOG_DIR}
#PBS -M {EMAIL}
#PBS -m {EMAIL_OPTIONS}

. {CONFIG_PATH}

echo "Running GATK HaplotypeCaller..."
java -jar $gatk_dir/GenomeAnalysisTK.jar \\
	-T HaplotypeCaller \\
	-R $ref_fasta \\
{gatk_input_combined_samples_str}
	-o $main_dir/{GROUP_NAME}.raw.snps.indels.vcf

echo "Done!"
'''

f = open("qsub_2_call_variants", 'w')
f.write(qsub_2_call_variants_str)
f.close()

qsub_1_align_str = \
f'''#!/bin/sh
#PBS -A winzeler-group
#PBS -N {JOB_NAME}_align
#PBS -l {RESOURCES_STR}
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

bwa mem -M -t 16 $ref_fasta $fastq_dir/${sample}_R1.fastq.gz $fastq_dir/${sample}_R2.fastq.gz > $main_dir/${sample}.sam

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