1. Download .fastq files, make an output directory

2. Update relevant paths in config.cfg and source it (`. config.cfg`)

3. Run `python rename_fastq_generate_samples.py` (renames .fastq files, creates samples.txt), or manually rename and make samples.txt (helpful tool: genomics_utils/filesystem_utils.py)
   Expected fastq format: {species_abbr}_{sample}_R{1,2}.fastq.gz
	 Example: p_fal_GHDDI-DSM265R-Dd2-E9_R1.fastq.gz

4. Run `python gatk_pipeline_qsub_generator.py` (creates qsub_1_align and qsub_2_call_variants)

5. Check and submit qsub_1_align (creates .ready.bam)

6. Check and submit qsub_2_call_variants (creates .raw.snps.indels.vcf)

7. Run subtraction step (creates GATK-Filters.txt)

Make sure to backup .fastq, .ready.bam, .raw.snps.indels.vcf and GATK-Filters on the NAS