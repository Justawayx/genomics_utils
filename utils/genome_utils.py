from collections import defaultdict
from urllib.parse import unquote

SPECIES_SNPEFF_ID_DICT = {"p_fal": "Pf3D7v3"}

def get_gene_info_dicts(WDIR, species="p_fal"):
	
	if species == "p_fal":
		
		gff_fpath = "%s/GENOME_RESOURCES/pf/p_fal_ref/p_fal.gff" % WDIR
		
		chromosomes = []
		for line in open(gff_fpath, 'r'):
				if line[0] != '#':
						break
				if line[:17] == '##sequence-region':
						chromosomes.append(line.strip().split()[1])
		
		chrom_gene_exon_interval_dict = {chrom: defaultdict(dict) for chrom in chromosomes} # chrom -> gene ID -> exon ID -> (start, end, strand_direction)

		chrom_gene_ids_dict = {chrom: set() for chrom in chromosomes} # All protein coding gene IDs
		gene_desc_dict = {chrom: {} for chrom in chromosomes} # gene_id -> description
		gene_interval_dict = {chrom: {} for chrom in chromosomes} # gene_id -> (start, end, strand_direction)

		for line in open(gff_fpath, 'r'):
				if line.strip() == '##FASTA':
						break
				if line[0] == '#':
						continue
				chrom, source, feature_type, start_pos, end_pos, _, sdirection, _, info = line.strip().split('\t')
				start_pos = int(start_pos); end_pos = int(end_pos)
				if feature_type == 'exon':
						exon_id = info.split(';')[0].split('ID=')[1]
						gene_id = exon_id.split('exon_')[1].split('-')[0]
						chrom_gene_exon_interval_dict[chrom][gene_id][exon_id] = (start_pos, end_pos, sdirection)
				if feature_type == 'gene':
						gene_id = info.split(';')[0].split('ID=')[1]
						gene_desc = info.split(';')[2].split('description=')[1]
						gene_desc_dict[gene_id] = unquote(gene_desc.strip()).replace('+', ' ')
						gene_interval_dict[gene_id] = (start_pos, end_pos, sdirection)
				if feature_type == 'CDS': # Actually coding
						gene_id = info.split(';')[0].split('ID=')[1].split('cds_')[1].split('-')[0]
						if gene_id != 'PF3D7_0112400' and 'pseudogene' not in gene_desc_dict[gene_id]: # Ignore pseudogenes
								chrom_gene_ids_dict[chrom].add(gene_id)
		
		return chrom_gene_ids_dict, gene_desc_dict, gene_interval_dict