from file_parsing_utils import parse_VCF_simple

# TODO: replace with genome_utils

CHROMOSOME_SEQ_DICT = {}

with open('/storage/NFS/GENOME_RESOURCES/pf/p_fal_ref/p_fal.fasta', 'r') as f:
    first = False
    for line in f:
        if line[0] == '>':
            
            if first is True:
                CHROMOSOME_SEQ_DICT[chrom] = seq
            
            first = True
            chrom = line[1:].split(' | ')[0]
            seq = ''
        else:
            seq += line.strip('\n')
    
    CHROMOSOME_SEQ_DICT[chrom] = seq

# ============================================
# Tools for identifying and comparing strains
# ============================================

SNP_ID_CHROM_POS_DICT = {} # Barcoding SNP ID -> (chrom, pos)
BARCODE_SNP_IDs = []
BARCODE_LOCI = []

with open('data/Barcoding_SNPs_Pf3D7v3_translation.tsv', 'r') as f:
    header_items = f.readline().strip('\n').split('\t')
    for line in f:
        SNP_ID, old_gene, old_chrom, old_pos, gene, chrom, pos, desc = line.strip('\n').split('\t')
        pos = int(pos)
        SNP_ID_CHROM_POS_DICT[SNP_ID] = (chrom, pos)
        BARCODE_LOCI.append((chrom, pos))
        BARCODE_SNP_IDs.append(SNP_ID)

SNP_3D7_ALLELE_DICT = {}

for SNP_ID in SNP_ID_CHROM_POS_DICT:
    chrom, pos = SNP_ID_CHROM_POS_DICT[SNP_ID]
    allele_3D7 = CHROMOSOME_SEQ_DICT[chrom][pos-1]
    SNP_3D7_ALLELE_DICT[SNP_ID] = allele_3D7

STRAIN_SNP_BARCODE_DICT = defaultdict(dict)

with open('data/common_strain_SNP_barcodes.tsv', 'r') as f:
    header_items = f.readline().strip('\n').split('\t')
    SNP_IDs = header_items[1:]
    for line in f:
        items = line.strip('\n').split('\t')
        strain = items[0]
        for SNP_ID, allele in zip(SNP_IDs, items[1:]):
            STRAIN_SNP_BARCODE_DICT[strain][SNP_ID] = allele

def find_closest_strain_by_barcode(sample, vcf_filepath):
    
    sample_locus_variant_dict = parse_VCF_simple(vcf_filepath, BARCODE_LOCI)
    
    sample_SNP_barcode_dict = defaultdict(dict)
    
    for locus, SNP_ID in zip(BARCODE_LOCI, BARCODE_SNP_IDs):
        if locus in sample_locus_variant_dict[sample]:
            REF, ALT, GT, AD, DP = sample_locus_variant_dict[sample][locus]
            alleles = [REF] + ALT.split(',')
            if GT == './.':
                major_allele = '.' # Undetermined
            else:
                allele_nums = set([int(num) for num in GT.split('/')])
                if len(allele_nums) == 1:
                    major_allele = alleles[list(allele_nums)[0]]
                else:
                    allele_depths = [int(ad) for ad in AD.split(',')]
                    major_allele, major_allele_depth = max(zip(alleles, allele_depths), key=lambda x: x[1])
        else:
            major_allele = SNP_3D7_ALLELE_DICT[SNP_ID] # Assume matches 3D7 reference
        
        sample_SNP_barcode_dict[sample][SNP_ID] = major_allele
    
    # Simply count number of mismatches
    strain_mismatches_dict = {}
    
    for strain in STRAIN_SNP_BARCODE_DICT:
        strain_mismatches_dict[strain] = 0
        for SNP_ID in BARCODE_SNP_IDs:
            ref_allele = STRAIN_SNP_BARCODE_DICT[strain][SNP_ID]
            if sample_SNP_barcode_dict[sample][SNP_ID] != ref_allele:
                strain_mismatches_dict[strain] += 1
    
    return sorted(strain_mismatches_dict.items(), key=lambda x: x[1])