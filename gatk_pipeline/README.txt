1. Download .fastq files, make an output directory
2. Update relevant paths in config.cfg
3. Run `python rename_fastq_generate_samples.py` (renames .fastq files, creates samples.txt)
4. Run `python gatk_pipeline_qsub_generator.py` (creates qsub_1_align and qsub_2_call_variants)
5. Check and submit qsub_1_align (creates .ready.bam)
6. Check and submit qsub_2_call_variants (creates .raw.snps.indels.vcf)
7. Run subtraction step [TODO]
7. Backup .fastq, .ready.bam, .raw.snps.indels.vcf and GATK-Filters on the NAS