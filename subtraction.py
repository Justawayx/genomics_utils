from collections import defaultdict
import os
from urllib.parse import unquote

# =======================
# Parameters (CHANGE)
# =======================

WDIR = "/storage/NFS" # /projects/winzeler if running on TSCC, /storage/NFS if on winzelerserver
VCF_DIR = f"{WDIR}/ROTATION_PROJECT/daisy/sarah_yliu4/data" # Where VCF is located
OUTPUT_DIR = f"{WDIR}/ROTATION_PROJECT/daisy/sarah_yliu4/data" # Where output table should be written
VCF = "Yliu4-combined.raw.snps.indels.vcf" # Name of VCF file (output of qsub_2_call_variants)
OUTPUT_NAME = "Yliu4-combined" # Name of output (selection group name)
PARENT_SAMPLE = "ParentDd2" # Name of parent sample

SNPEFF_DIR = f"{WDIR}/GENOME_RESOURCES/snpEff"

def get_format_data_dict(FORMAT, data):
    return {key: val for key, val in zip(FORMAT.split(':'), data.split(':'))}

# =======================
# Load gene annotations
# =======================

# TODO: integrate with utils

gff_fpath = f"{WDIR}/GENOME_RESOURCES/pf/p_fal_ref/p_fal.gff"

chromosomes = []
for line in open(gff_fpath, 'r'):
    if line[0] != '#':
        break
    if line[:17] == '##sequence-region':
        chromosomes.append(line.strip().split()[1])

nuclear_chromosomes = sorted(chromosomes[2:])
chrom_gene_exon_interval_dict = {chrom: defaultdict(dict) for chrom in chromosomes} # chrom -> gene ID -> exon ID -> (start, end, strand_direction)

chrom_gene_ids_dict = {chrom: set() for chrom in chromosomes} # All protein coding gene IDs
chrom_gene_desc_dict = {chrom: {} for chrom in chromosomes} # gene_id -> description
chrom_gene_interval_dict = {chrom: {} for chrom in chromosomes} # gene_id -> (start, end, strand_direction)

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
        chrom_gene_desc_dict[chrom][gene_id] = unquote(gene_desc.strip()).replace('+', ' ')
        chrom_gene_interval_dict[chrom][gene_id] = (start_pos, end_pos, sdirection)
    if feature_type == 'CDS': # Actually coding
        gene_id = info.split(';')[0].split('ID=')[1].split('cds_')[1].split('-')[0]
        if gene_id != 'PF3D7_0112400' and 'pseudogene' not in chrom_gene_desc_dict[chrom][gene_id]: # Ignore pseudogenes
            chrom_gene_ids_dict[chrom].add(gene_id)

# =======================
# Run SnpEff on VCF
# =======================

os.chdir(OUTPUT_DIR)
os.system(f"java -jar {SNPEFF_DIR}/snpEff.jar -formatEFF -o vcf -ud 0 -c {SNPEFF_DIR}/snpEff.config Pf3D7v3 {VCF_DIR}/{VCF} > {OUTPUT_DIR}/{VCF}.ann.txt")

# =======================
# Load data from VCF
# =======================

f = open(f"{OUTPUT_DIR}/{VCF}.ann.txt")

child_samples = []
for line in f:
    if line.startswith('#CHROM'):
        header_items = line.strip().split('\t')
        print('\t'.join(header_items))
        for item in header_items[9:]:
            if item != PARENT_SAMPLE:
                child_samples.append(item)
        break

chrom_pos_sample_data_dict = defaultdict(dict) # chrom -> pos -> sample -> data
chrom_pos_anno_dict = defaultdict(dict) # chrom -> pos -> (ref, alt, gene, type, effect, impact, codon_change, aa_change)

chrom_pos_quality_dict = defaultdict(dict)
parent_col_idx = header_items.index(PARENT_SAMPLE)

for line in f:
    items = line.strip().split('\t')
    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = items[:9]
    POS = int(POS)
    
    if float(QUAL) < 500:
        continue
    
    parent_data = items[parent_col_idx]
    parent_data_dict = get_format_data_dict(FORMAT, parent_data)
    
    if parent_data_dict['GT'] == './.' or 'DP' not in parent_data_dict or \
       parent_data_dict['DP'] == '.' or int(parent_data_dict['DP']) < 7:
        continue
    
    chrom_pos_quality_dict[CHROM][POS] = QUAL
    chrom_pos_sample_data_dict[CHROM][POS] = {}
    non_parent_GTs = []
    for idx in range(9, len(items)):
        sample = header_items[idx]
        data = items[idx]
        data_dict = get_format_data_dict(FORMAT, data)
        GT = data_dict['GT']
        if GT not in ['0/0', '0/1', './.'] and parent_data_dict['GT'] == '0/0':
            pass # print(items)
        if idx != parent_col_idx:
            non_parent_GTs.append(GT)
        
        if 'DP' not in data_dict or data_dict['DP'] == '.' or int(data_dict['DP']) < 7:
            chrom_pos_sample_data_dict[CHROM][POS][sample] = (data_dict['GT'], 'LowDP')
            continue
        AD = data_dict['AD']
        chrom_pos_sample_data_dict[CHROM][POS][sample] = (GT, AD)
    
    if set(non_parent_GTs) == set(['0/0']):
        # print(non_parent_GTs)
        print('\t'.join(items))
        pass
    
    first_alt_allele = ALT.split(',')[0]
    
    vtype = 'SNP' if (REF in ['A', 'C', 'G', 'T'] and first_alt_allele in ['A', 'C', 'G', 'T']) else 'INDEL'        
    pre_eff = INFO.split(';EFF=')[0]
    eff = INFO.split(';EFF=')[1].split(';')[0]
    
    allele_anno_dict = {}
    for anno in eff.split(','):
        last = anno.split('|')[-1]
        if 'INFO_' in last or 'WARNING_' in last:
            alt_allele = anno.split('|')[-2].split(')')[0]
        else:
            alt_allele = anno.split('|')[-1].split(')')[0]
        effect, rest = anno.split('(')
        impact, codon_change, aa_change, _, gene = rest.split('|')[1:6]
        allele_anno_dict[alt_allele] = (effect, impact, codon_change, aa_change, gene)
    
    if first_alt_allele in allele_anno_dict:
        effect, impact, codon_change, aa_change, gene = allele_anno_dict[first_alt_allele]
        chrom_pos_anno_dict[CHROM][POS] = (REF, first_alt_allele, gene, vtype, effect, impact, codon_change, aa_change)
    else:
        print('\t'.join(items))

# =======================
# Write SNV/indel table
# =======================

ordered_samples = [PARENT_SAMPLE] + sorted(child_samples)
print('\n'.join(ordered_samples))

f = open(f"{OUTPUT_DIR}/{OUTPUT_NAME}_GATK-Filters.tsv", 'w')
header_items = ['Chromosome', 'Position', 'Gene_Name', 'Gene_Descrip', 'Quality', 'Ref_Base', 'Alt_Base', 'Type', 
                 'Effect', 'Impact', 'Codon_Change', 'Amino_Acid_Change']

for sample in ordered_samples:
    header_items.append("%s_MutationCall" % sample)

for sample in ordered_samples:
    header_items.append("%s_AlleleDepths" % sample)

f.write('\t'.join(header_items) + '\n')

new_chrom_pos_sets = defaultdict(set)
for chrom in chrom_pos_sample_data_dict:
    for pos in sorted(chrom_pos_sample_data_dict[chrom]):
        if chrom_pos_sample_data_dict[chrom][pos][PARENT_SAMPLE][0] != '0/0':
            continue
        quality = chrom_pos_quality_dict[chrom][pos]
        ref, alt, gene, vtype, effect, impact, codon_change, aa_change = chrom_pos_anno_dict[chrom][pos]
        gene_desc = chrom_gene_desc_dict[chrom][gene] if gene in chrom_gene_desc_dict[chrom] else ''
        
        items = [chrom, pos, gene, gene_desc, quality, ref, alt, vtype, effect, impact, codon_change, aa_change]
        
        for sample in ordered_samples:
            GT = chrom_pos_sample_data_dict[chrom][pos][sample][0]
            items.append(GT)
        
        for sample in ordered_samples:
            AD = chrom_pos_sample_data_dict[chrom][pos][sample][1]
            items.append(AD)
        
        f.write('\t'.join([str(item) for item in items]) + '\n')
        
        new_chrom_pos_sets[chrom].add(pos)

f.close()