1. Download .fastq files, make an output directory

	To download from IGM FTP server:
	wget --user="winzeler" --password="azvWp164Rh" -r ftp://igm-storage.ucsd.edu/240814_LH00444_0175_B22N3YNLT3/

2. Update relevant paths in config.cfg and source it (`. config.cfg`)

3. Run `python rename_fastq_generate_samples.py` (renames .fastq files, creates samples.txt), or manually rename and make samples.txt (helpful tool: genomics_utils/filesystem_utils.py)

	Expected fastq format: {species_abbr}_{sample}_R{1,2}.fastq.gz
	Example: p_fal_GHDDI-DSM265R-Dd2-E9_R1.fastq.gz

4. Run `python gatk_pipeline_slurm_generator.py` (creates sbatch_0_bwa, sbatch_1_align, and multiple sbatch_2_call_variants-partXX)

5. Check and submit sbatch_0_bwa (creates .sam)

5. Check and submit sbatch_1_align (creates .ready.bam)

6. Check and submit sbatch_2_call_variants (creates .raw.snps.indels.vcf)

	At the same time, if strain is 3D7 or Dd2, you can submit sbatch_cnv_analysis

7. Run `python combine_vcf_subtract_parent.py` (creates SNV-INDELs.tsv)

	This step combines the partial VCFs into one VCF, runs SnpEff on that full VCF, and computes P values for difference between parent and child variant calls.
	You will most likely need to request an interactive job before running this step, which can be done on TSCC as follows:

	srun --account=htl124 --partition=hotel --qos=hotel --pty --nodes=1 --ntasks-per-node=1 -t 00:10:00 --export=ALL /bin/bash

	If sbatch_cnv_analysis was run successfully, run `python combined_denoisedCR.py` (creates CNV_analysis.tsv)

Make sure to backup .fastq, .ready.bam, .raw.snps.indels.vcf and SNV-INDELs.tsv / CNV_analysis.tsv on the NAS
