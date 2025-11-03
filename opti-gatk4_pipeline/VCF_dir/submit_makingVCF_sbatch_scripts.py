import subprocess

for i in range(5, 5+1):
    with open('../Ref_data/core_chr%i.list' % i, 'r') as f:
        for line in f:
            region = line.strip('\n')
            if region != 'Pf3D7_05_v3:400001-700000' and region != 'Pf3D7_05_v3:700001-1000000':
                continue
            fname = 'sbatch_genotype_chr%i_%s' % (i, region)
            
            script = f'''#!/bin/bash
#SBATCH --job-name VCF_chr{i}_{region}
#SBATCH --partition=hotel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10GB
#SBATCH --time=26:00:00
#SBATCH --qos=hotel
#SBATCH --account=htl124
#SBATCH --export=ALL
#SBATCH --exclude tscc-11-11
#SBATCH --output /tscc/lustre/ddn/scratch/dwc001/logs/slurm-%A_%a.%x.out-%N
#SBATCH --output /tscc/lustre/ddn/scratch/dwc001/logs/slurm-%A_%a.%x.err-%N

module load shared
module load gatk
module load samtools
module load bcftools

ref_dir=/tscc/projects/ps-winzelerlab/SEQUENCING/Winzeler/Jason/LNA27_Selection_Dd2_Optimized_GATK4_pipeline_test/Ref_data
vcf_dir=/tscc/projects/ps-winzelerlab/SEQUENCING/Winzeler/Jason/LNA27_Selection_Dd2_Optimized_GATK4_pipeline_test/VCF_dir
temp_dir=/tscc/projects/ps-winzelerlab/SEQUENCING/Winzeler/Jason/LNA27_Selection_Dd2_Optimized_GATK4_pipeline_test/VCF_dir/temp

gatk --java-options "-Xmx80g -Xms80g"  GenotypeGVCFs --genomicsdb-use-bcf-codec true -R $ref_dir/Pf3D7.fasta -V gendb://$vcf_dir/chr{i}_database --max-genotype-count 1024 -O $vcf_dir/chr{i}_part{region}.vcf.gz --tmp-dir $temp_dir -stand-call-conf 30 -L {region}
'''
            
            result = subprocess.run(["sbatch"], input=script, text=True, capture_output=True)
            print("Submitted %s" % fname)
