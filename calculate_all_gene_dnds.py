import math, random
from collections import defaultdict
import pickle, gzip, csv
from urllib.parse import unquote

from utils.dnds_utils import (classify_snv, nuc_sequence_to_expected_n_s, nuc_sequence_to_aa_sequence)

comp_dict = {'G': 'C', 'C': 'G', 'T': 'A', 'A': 'T', 'g': 'c', 'c': 'g', 't': 'a', 'a': 't'}

# Load gene information

gff_fpath = "/storage/NFS/GENOME_RESOURCES/pf/p_fal_ref/p_fal.gff"

chromosomes = []
for line in open(gff_fpath, 'r'):
	if line[0] != '#':
		break
	if line[:17] == '##sequence-region':
		chromosomes.append(line.strip().split()[1])

chrom_gene_cds_interval_dict = {chrom: defaultdict(dict) for chrom in chromosomes} # chrom -> gene ID -> CDS ID -> (start, end, strand_direction)

all_coding_gene_ids = set() # All protein coding gene IDs

chrom_gene_desc_dict = {chrom: {} for chrom in chromosomes} # gene_id -> description
chrom_gene_interval_dict = {chrom: {} for chrom in chromosomes} # gene_id -> (start, end, strand_direction)

for line in open(gff_fpath, 'r'):
	if line.strip() == '##FASTA':
		break
	if line[0] == '#':
		continue
	chrom, source, feature_type, start_pos, end_pos, _, sdirection, _, info = line.strip().split('\t')
	start_pos = int(start_pos); end_pos = int(end_pos)
	if feature_type == 'CDS':
		cds_id = info.split(';')[0].split('ID=')[1]
		gene_id = cds_id.split('cds_')[1].split('-')[0]
		all_coding_gene_ids.add(gene_id)
		chrom_gene_cds_interval_dict[chrom][gene_id][cds_id] = (start_pos, end_pos, sdirection)
	if feature_type == 'gene':
		gene_id = info.split(';')[0].split('ID=')[1]
		gene_desc = info.split(';')[2].split('description=')[1]
		chrom_gene_desc_dict[chrom][gene_id] = unquote(gene_desc.strip()).replace('+', ' ')
		chrom_gene_interval_dict[chrom][gene_id] = (start_pos, end_pos, sdirection)
	

def reverse_complement(dna):
    complement = ''
    for nuc in dna[::-1]:
        complement += comp_dict[nuc]
    return complement

# Gene of interest

gene_name = 'PF3D7_1133400'
gene_coords = ('Pf3D7_11_v3', 1292966, 1296696)
gene_chrom, gene_start, gene_end = gene_coords

orig_full_nuc_sequence = "aataatatatatatatatatatttaaatttatggattaataacttcaatgattttccaatttttcatttcgttttatttattttatatttttatatgaattgatttaattttgttctattttatttttttgttttcgtagtggtaatttaatcccttgttttgtatatatagcatgcagtaagggaaaaaaaaaaaaaaaaaagttgttagagaaaacagatcattttgtttaaataaataaatatatacatatatatatatatatatatatattaatatttgtgtttttcttttttatattcgatggaatttattgattttattcgataattataactaataaatgtattaaaaatatattaaattaaaaaaaaaaaattataaataaatatatatatatatatataaatatatatataacaccttgttttgaaacctttacacaaacgttatacgtacacaggtttttcttttatttttatgtgcttcttttttatttcacttttgttagagtctcttattaaacgttaaaaaaaaaaaaaaaaaaaaaaattactaaataataaaagatctaataacggatatatctattttttttaatgcacaagaaaaaaaaaaaaaagaacaagggaaaaaagaagaaaaatggaaattaaaattttacacttattaataaatatttaatatatttatatattataaaaaaaagaaaaaaaaaaaaaaaaaaaaaaaattaaccatatttattggttttatttttttatttttttaaaaaaaaaaattgagaatataaattgtaatatattttttattttaatataattttaaaaacctaataatttatttgataatttttcaaattaatgtacttgttataaattgtacaaaaATGAGAAAATTATACTGCGTATTATTATTGAGCGCCTTTGAGTTTACATATATGATAAACTTTGGAAGAGGACAGAATTATTGGGAACATCCATATCAAAATAGTGATGTGTATCGTCCAATCAACGAACATAGGGAACATCCAAAAGAATACGAATATCCATTACACCAGGAACATACATACCAACAAGAAGATTCAGGAGAAGACGAAAATACATTACAACACGCATATCCAATAGACCACGAAGGTGCCGAACCCGCACCACAAGAACAAAATTTATTTTCAAGCATTGAAATAGTAGAAAGAAGTAATTATATGGGTAATCCATGGACGGAATATATGGCAAAATATGATATTGAAGAAGTTCATGGTTCAGGTATAAGAGTAGATTTAGGAGAAGATGCTGAAGTAGCTGGAACTCAATATAGACTTCCATCAGGGAAATGTCCAGTATTTGGTAAAGGTATAATTATTGAGAATTCAAATACTACTTTTTTAACACCGGTAGCTACGGGAAATCAATATTTAAAAGATGGAGGTTTTGCTTTTCCTCCAACAGAACCTCTTATGTCACCAATGACATTAGATGAAATGAGACATTTTTATAAAGATAATAAATATGTAAAAAATTTAGATGAATTGACTTTATGTTCAAGACATGCAGGAAATATGATTCCAGATAATGATAAAAATTCAAATTATAAATATCCAGCTGTTTATGATGACAAAGATAAAAAGTGTCATATATTATATATTGCAGCTCAAGAAAATAATGGTCCTAGATATTGTAATAAAGACGAAAGTAAAAGAAACAGCATGTTTTGTTTTAGACCAGCAAAAGATATATCATTTCAAAACTATACATATTTAAGTAAGAATGTAGTTGATAACTGGGAAAAAGTTTGCCCTAGAAAGAATTTACAGAATGCAAAATTCGGATTATGGGTCGATGGAAATTGTGAAGATATACCACATGTAAATGAATTTCCAGCAATTGATCTTTTTGAATGTAATAAATTAGTTTTTGAATTGAGTGCTTCGGATCAACCTAAACAATATGAACAACATTTAACAGATTATGAAAAAATTAAAGAAGGTTTCAAAAATAAGAACGCTAGTATGATCAAAAGTGCTTTTCTTCCCACTGGTGCTTTTAAAGCAGATAGATATAAAAGTCATGGTAAGGGTTATAATTGGGGAAATTATAACACAGAAACACAAAAATGTGAAATTTTTAATGTCAAACCAACATGTTTAATTAACAATTCATCATACATTGCTACTACTGCTTTGTCCCATCCCATCGAAGTTGAAAACAATTTTCCATGTTCATTATATAAAGATGAAATAATGAAAGAAATCGAAAGAGAATCAAAACGAATTAAATTAAATGATAATGATGATGAAGGGAATAAAAAAATTATAGCTCCAAGAATTTTTATTTCAGATGATAAAGACAGTTTAAAATGCCCATGTGACCCTGAAATGGTAAGTAATAGTACATGTCGTTTCTTTGTATGTAAATGTGTAGAAAGAAGGGCAGAAGTAACATCAAATAATGAAGTTGTAGTTAAAGAAGAATATAAAGATGAATATGCAGATATTCCTGAACATAAACCAACTTATGATAAAATGAAAATTATAATTGCATCATCAGCTGCTGTCGCTGTATTAGCAACTATTTTAATGGTTTATCTTTATAAAAGAAAAGGAAATGCTGAAAAATATGATAAAATGGATGAACCACAAGATTATGGGAAATCAAATTCAAGAAATGATGAAATGTTAGATCCTGAGGCATCTTTTTGGGGGGAAGAAAAAAGAGCATCACATACAACACCAGTTCTGATGGAAAAACCATACTATTAAaatgtgaactataataatttcaacgtctgatataatcagcttctcttttatgctaaaaaaaaaaaaaatatatatatatttataaatatatttatatatatttatatttatatttctatgatttcttaatatttattttttattttaaaaacacaaaaaaataattcacaaaaaatccattttatgttgtttctttactatattttttatgtattacacgaattaaaaaatgaatatacatatggatatataaatataaatatatatatatatatatattagttgttaatttatttatttttttattttttccttaaacaatttgatgttgtaactttaaagttaaagcatcataatgtaaattcttcttttacgccgtaaaaattcttatatatatatatatatatatatatatatataacatatttttattcttttaaaatttttaatatgttttcttttctttaaaataatattatttcaaacatttaatatatatatatatatatatatatatttataaatgtatataaatatatatattttttttttttactttaatataataattcataataactttaaatttataaatactgatattttctattttttttgttttccttgccatttataaaacaggcacacatataataaattttatatatgataagatatatatattatataactgacccttaaataaatatatatgattagaataatataaatttggtttattaaaaaataaaataaaaatgaagggttcataaaaatgttatgggaaaatgtgtacaacaatagatgtccattataaggaaaatgaaaggatatataaaatataaataaataaacatatatatatatatatatatatatatatattatattatgttattattcttttaaccatatttataaaattattatgacgtatgtgacattatatatataattaaaaaggaaaagaaa"
full_nuc_sequence = reverse_complement(orig_full_nuc_sequence)

orig_nuc_sequence = "ATGAGAAAATTATACTGCGTATTATTATTGAGCGCCTTTGAGTTTACATATATGATAAACTTTGGAAGAGGACAGAATTATTGGGAACATCCATATCAAAATAGTGATGTGTATCGTCCAATCAACGAACATAGGGAACATCCAAAAGAATACGAATATCCATTACACCAGGAACATACATACCAACAAGAAGATTCAGGAGAAGACGAAAATACATTACAACACGCATATCCAATAGACCACGAAGGTGCCGAACCCGCACCACAAGAACAAAATTTATTTTCAAGCATTGAAATAGTAGAAAGAAGTAATTATATGGGTAATCCATGGACGGAATATATGGCAAAATATGATATTGAAGAAGTTCATGGTTCAGGTATAAGAGTAGATTTAGGAGAAGATGCTGAAGTAGCTGGAACTCAATATAGACTTCCATCAGGGAAATGTCCAGTATTTGGTAAAGGTATAATTATTGAGAATTCAAATACTACTTTTTTAACACCGGTAGCTACGGGAAATCAATATTTAAAAGATGGAGGTTTTGCTTTTCCTCCAACAGAACCTCTTATGTCACCAATGACATTAGATGAAATGAGACATTTTTATAAAGATAATAAATATGTAAAAAATTTAGATGAATTGACTTTATGTTCAAGACATGCAGGAAATATGATTCCAGATAATGATAAAAATTCAAATTATAAATATCCAGCTGTTTATGATGACAAAGATAAAAAGTGTCATATATTATATATTGCAGCTCAAGAAAATAATGGTCCTAGATATTGTAATAAAGACGAAAGTAAAAGAAACAGCATGTTTTGTTTTAGACCAGCAAAAGATATATCATTTCAAAACTATACATATTTAAGTAAGAATGTAGTTGATAACTGGGAAAAAGTTTGCCCTAGAAAGAATTTACAGAATGCAAAATTCGGATTATGGGTCGATGGAAATTGTGAAGATATACCACATGTAAATGAATTTCCAGCAATTGATCTTTTTGAATGTAATAAATTAGTTTTTGAATTGAGTGCTTCGGATCAACCTAAACAATATGAACAACATTTAACAGATTATGAAAAAATTAAAGAAGGTTTCAAAAATAAGAACGCTAGTATGATCAAAAGTGCTTTTCTTCCCACTGGTGCTTTTAAAGCAGATAGATATAAAAGTCATGGTAAGGGTTATAATTGGGGAAATTATAACACAGAAACACAAAAATGTGAAATTTTTAATGTCAAACCAACATGTTTAATTAACAATTCATCATACATTGCTACTACTGCTTTGTCCCATCCCATCGAAGTTGAAAACAATTTTCCATGTTCATTATATAAAGATGAAATAATGAAAGAAATCGAAAGAGAATCAAAACGAATTAAATTAAATGATAATGATGATGAAGGGAATAAAAAAATTATAGCTCCAAGAATTTTTATTTCAGATGATAAAGACAGTTTAAAATGCCCATGTGACCCTGAAATGGTAAGTAATAGTACATGTCGTTTCTTTGTATGTAAATGTGTAGAAAGAAGGGCAGAAGTAACATCAAATAATGAAGTTGTAGTTAAAGAAGAATATAAAGATGAATATGCAGATATTCCTGAACATAAACCAACTTATGATAAAATGAAAATTATAATTGCATCATCAGCTGCTGTCGCTGTATTAGCAACTATTTTAATGGTTTATCTTTATAAAAGAAAAGGAAATGCTGAAAAATATGATAAAATGGATGAACCACAAGATTATGGGAAATCAAATTCAAGAAATGATGAAATGTTAGATCCTGAGGCATCTTTTTGGGGGGAAGAAAAAAGAGCATCACATACAACACCAGTTCTGATGGAAAAACCATACTATTAA"
orig_sequence_start_within_full = orig_full_nuc_sequence.index(orig_nuc_sequence[0])

nuc_sequence = reverse_complement(orig_nuc_sequence)
nuc_sequence_start_within_full = full_nuc_sequence.index(nuc_sequence[0])

gene_trans_start = gene_start+nuc_sequence_start_within_full
gene_trans_end = gene_trans_start+len(nuc_sequence)-1

aa_sequence = "MRKLYCVLLLSAFEFTYMINFGRGQNYWEHPYQNSDVYRPINEHREHPKEYEYPLHQEHTYQQEDSGEDENTLQHAYPIDHEGAEPAPQEQNLFSSIEIVERSNYMGNPWTEYMAKYDIEEVHGSGIRVDLGEDAEVAGTQYRLPSGKCPVFGKGIIIENSNTTFLTPVATGNQYLKDGGFAFPPTEPLMSPMTLDEMRHFYKDNKYVKNLDELTLCSRHAGNMIPDNDKNSNYKYPAVYDDKDKKCHILYIAAQENNGPRYCNKDESKRNSMFCFRPAKDISFQNYTYLSKNVVDNWEKVCPRKNLQNAKFGLWVDGNCEDIPHVNEFPAIDLFECNKLVFELSASDQPKQYEQHLTDYEKIKEGFKNKNASMIKSAFLPTGAFKADRYKSHGKGYNWGNYNTETQKCEIFNVKPTCLINNSSYIATTALSHPIEVENNFPCSLYKDEIMKEIERESKRIKLNDNDDEGNKKIIAPRIFISDDKDSLKCPCDPEMVSNSTCRFFVCKCVERRAEVTSNNEVVVKEEYKDEYADIPEHKPTYDKMKIIIASSAAVAVLATILMVYLYKRKGNAEKYDKMDEPQDYGKSNSRNDEMLDPEASFWGEEKRASHTTPVLMEKPYY"

# In[4]:


n_denom, s_denom = nuc_sequence_to_expected_n_s(orig_nuc_sequence)
print(n_denom, s_denom)

# In[5]:


def chrom_pos_to_pos_within_gene(pos):
    return pos - gene_start

def chrom_pos_to_pos_within_nuc_sequence(pos):
    return pos - (gene_start+nuc_sequence_start_within_full)

def chrom_pos_to_pos_within_orig_nuc_sequence(pos):
    pos_within_nuc_sequence = chrom_pos_to_pos_within_nuc_sequence(pos)    
    return ((len(nuc_sequence)-1) - pos_within_nuc_sequence)

# # Metadata

# In[6]:


# Mozambique metadata

moz_vcf_dir = '/home/dwc001/moz_vcf' # vcf files: name.g.vcf
moz_data_dir = '/storage/NFS/ROTATION_PROJECT/daisy/mozambique'
metadata_path = '%s/Moz_metadata_subset.tsv' % moz_data_dir
metadata_df = pd.read_csv(metadata_path, sep='\t', header=0)

moz_samples = ['Moz-%i' % i for i in range(1, 121)] # All sample names of the format Moz-n

study_id_sample_dict = {}
f = open('%s/Moz_sample_mapping.tsv' % moz_data_dir, 'r')
for line in f:
    sample_num, study_id = line.strip().split('\t')
    study_id_sample_dict[study_id] = 'Moz-%s' % sample_num

# Create metadata dicts

sample_locale_dict = {}
sample_hiv_dict = {}
sample_date_dict = {}

for study_id in study_id_sample_dict:
    sample = study_id_sample_dict[study_id]
    locale_datum = metadata_df[metadata_df['Study ID'] == study_id]['Local of residence']
    if not locale_datum.empty:
        locale = locale_datum.item()
        sample_locale_dict[sample] = locale
    
    hiv_datum = metadata_df[metadata_df['Study ID'] == study_id]['HIV ELISA results']
    if not hiv_datum.empty:
        sample_hiv_dict[sample] = (hiv_datum.item() == 1)
    
    admit_date_datum = metadata_df[metadata_df['Study ID'] == study_id]['Date of Admission']
    if not admit_date_datum.empty:
        sample_date_dict[sample] = admit_date_datum.item()
    
    metadata_df['Date of Admission']

# In[7]:


# Pf6 metadata

# Sample metadata

f = open('/storage/NFS/ROTATION_PROJECT/daisy/pf6/File2_Pf_6_samples.txt', 'r')
header = f.readline().split('\t')

pf6_samples = []
sample_country_dict = {}
country_samples_dict = defaultdict(list)
sample_year_dict = {}
sample_qc_pass_dict = {}

for line in f:
    items = line.strip().split('\t')
    sample, study, site, country, lat, long, year, ena, sameind, pop, pctcall, qc, _, _ = items
    year = int(year) if items[6] != 'Lab' else 'Lab'
    
    pf6_samples.append(sample)
    sample_country_dict[sample] = country
    country_samples_dict[country].append(sample)
    sample_year_dict[sample] = year
    sample_qc_pass_dict[sample] = True if qc == 'True' else False

for sample in moz_samples:
    sample_country_dict[sample] = 'Mozambique'
    country_samples_dict['Mozambique'].append(sample)
    if sample in sample_date_dict:
        sample_year_dict[sample] = int('20' + sample_date_dict[sample].split('/')[2])
    else:
        sample_year_dict[sample] = 2017 # Assumption

# In[8]:


# Country HIV metadata
country_hiv_prev_dict = {}

f = open('%s/country_HIV_prev.tsv' % moz_data_dir, 'r')
header = f.readline()
for line in f:
    country, hiv_prev, hiv_count, hiv_deaths, year = line.strip().split('\t')
    if hiv_prev != '-':
        country_hiv_prev_dict[country] = float(hiv_prev[:-1])

# In[9]:


all_samples = moz_samples + pf6_samples

# # Load data at gene of interest

# In[10]:


def convert_format_dict(format_str, format_items_str):
    format_cats = format_str.split(':')
    format_items = format_items_str.split(':')
    return {cat: item for cat, item in zip(format_cats, format_items)}

def convert_info_dict(info_str):
    info_dict = {}
    for item in info_str.split(';'):
        parts = item.split('=')
        if len(parts) == 2:
            cat, val = parts
        else:
            cat = parts[0]; val = ''
        info_dict[cat] = val
    return info_dict

# ## Pickling below, skip

# In[11]:


# sample -> chrom -> pos -> record
sample_chrom_pos_record_dict = {sample: defaultdict(dict) for sample in all_samples}

# chrom -> pos -> SnpEff info
chrom_pos_anno_dict = defaultdict(dict)

# In[12]:


# Initial thresholds

mq_threshold = 60 # "MQ < 60.0" --filter-name "LowMQ" \
readposranksum_threshold = -8 # "ReadPosRankSum < -8.0" --filter-name "LowRPRS" \
quality_threshold = 500 # "QUAL < 500" --filter-name "LowQual" \
qd_threshold = 2 # "QD < 2" --filter-name "LowQD" \
sor_threshold = 4 # "SOR > 4.0" --filter-name "highSOR" \
dp_threshold = 7 # "DP < 7" --genotype-filter-name "LowFormatDP" \
fs_threshold = 60 # "FS > 60.0" --filter-name "highFS" \
mqranksum_threshold = -12.5 # "MQRankSum < -12.5" --filter-name "lowMQRankSum"

def check_thresholds(info_dict):
    if 'MQ' in info_dict and float(info_dict['MQ']) < mq_threshold:
        return False
    if 'ReadPosRankSum' in info_dict and float(info_dict['ReadPosRankSum']) < readposranksum_threshold:
        return False
    if 'QD' in info_dict and float(info_dict['QD']) < qd_threshold:
        return False
    if 'SOR' in info_dict and float(info_dict['SOR']) > sor_threshold:
        return False
    if 'DP' in info_dict and float(info_dict['DP']) < dp_threshold:
        return False
    if 'FS' in info_dict and float(info_dict['FS']) > fs_threshold:
        return False
    if 'MQRankSum' in info_dict and float(info_dict['MQRankSum']) < mqranksum_threshold:
        return False
    return True

# In[13]:


# Pf6

f = open('/storage/NFS/ROTATION_PROJECT/daisy/pf6/Pf_60_public_%s.ann.txt' % gene_chrom, 'r')

# Skip the meta-information lines
for line in f:
    if line[:2] !=  '##':
        break

samples = line.strip().split('\t')[9:] # ignore header    

for line in tqdm(f):
    items = line.strip().split('\t')
    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = items[:9]
    
    alleles = [REF] + ALT.split(',')
    POS = int(POS); QUAL = float(QUAL) if QUAL != '.' else QUAL
    
    if POS >= gene_start and POS <= gene_end:
        if FILTER == 'PASS':
            chrom_pos_anno_dict[CHROM][POS] = convert_info_dict(INFO)
            for tup, sample in zip(items[9:], samples):
                sample_chrom_pos_record_dict[sample][CHROM][POS] = (REF, ALT, QUAL, convert_format_dict(FORMAT, tup))
    elif POS > gene_end:
        break

# In[14]:


pickle.dump(sample_chrom_pos_record_dict, open('%s/AMA1_sample_chrom_pos_record_dict_v2.pkl' % moz_data_dir, 'wb'))
pickle.dump(chrom_pos_anno_dict, open('%s/AMA1_chrom_pos_anno_dict_v2.pkl' % moz_data_dir, 'wb'))

# In[15]:


sample_chrom_pos_record_dict = pickle.load(open('%s/AMA1_sample_chrom_pos_record_dict_v2.pkl' % moz_data_dir, 'rb'))
chrom_pos_anno_dict = pickle.load(open('%s/AMA1_chrom_pos_anno_dict_v2.pkl' % moz_data_dir, 'rb'))

# In[16]:


# Moz

moz_chrom_pos_anno_dict = defaultdict(dict)

f = open('/storage/NFS/MOZAMBIQUE/vcf_ann/Moz120_raw.snps.indels.ann.txt', 'r')

# Skip the meta-information lines
for line in f:
    if line[:2] !=  '##':
        break

samples = line.strip().split('\t')[9:] # ignore header    

for line in tqdm(f):
    items = line.strip().split('\t')
    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = items[:9]
    
    alleles = [REF] + ALT.split(',')
    POS = int(POS); QUAL = float(QUAL) if QUAL != '.' else QUAL
    
    if CHROM == gene_chrom:
        if POS >= gene_start and POS <= gene_end:            
            info_dict = convert_info_dict(INFO)
            if QUAL >= quality_threshold and check_thresholds(info_dict): # Initial filtering
                moz_chrom_pos_anno_dict[CHROM][POS] = convert_info_dict(INFO)
                for tup, sample in zip(items[9:], samples):
                    sample_chrom_pos_record_dict[sample][CHROM][POS] = (REF, ALT, QUAL, convert_format_dict(FORMAT, tup))
        elif POS > gene_end:
            break

# In[17]:


pickle.dump(sample_chrom_pos_record_dict, open('%s/AMA1_sample_chrom_pos_record_dict_v2.pkl' % moz_data_dir, 'wb'))
pickle.dump(moz_chrom_pos_anno_dict, open('%s/AMA1_moz_chrom_pos_anno_dict_v2.pkl' % moz_data_dir, 'wb'))

# ## Continue here

# In[18]:


sample_chrom_pos_record_dict = pickle.load(open('%s/AMA1_sample_chrom_pos_record_dict_v2.pkl' % moz_data_dir, 'rb'))
pf6_chrom_pos_anno_dict = pickle.load(open('%s/AMA1_chrom_pos_anno_dict_v2.pkl' % moz_data_dir, 'rb'))
moz_chrom_pos_anno_dict = pickle.load(open('%s/AMA1_moz_chrom_pos_anno_dict_v2.pkl' % moz_data_dir, 'rb'))

# In[19]:


# Thresholds for calling a SNV

alt_depth_threshold = 30
alt_freq_threshold = 0.3

# In[20]:


# Determine dN/dS per country

# country -> snv_type -> set of snvs of format (position, alt allele)
country_snvs_dict = {country: defaultdict(set) for country in country_region_dict}

for sample in sample_chrom_pos_record_dict:
    
    if sample in sample_qc_pass_dict and sample_qc_pass_dict[sample] == False: # Filter out samples
        continue
    
    country = sample_country_dict[sample]
    for pos in sample_chrom_pos_record_dict[sample][gene_chrom]:
        
        if pos < gene_trans_start or pos > gene_trans_end: # Skip UTRs
            continue
        
        REF, ALT, QUALITY, FITEMS = sample_chrom_pos_record_dict[sample][gene_chrom][pos]
        
        if FITEMS['DP'] == '.' or QUALITY < 500: # Bad coverage
            continue
        
        alt_alleles = ALT.split(',')
        alt_allele_depths = [int(val) for val in FITEMS['AD'].split(',')[1:]]
        total_depth = int(FITEMS['DP'])
        
        for alt_allele, alt_depth in zip(alt_alleles, alt_allele_depths):
            if REF in ['A', 'G', 'C', 'T'] and alt_allele in ['A', 'G', 'C', 'T']: # SNV
                if alt_depth >= alt_depth_threshold and (alt_depth/total_depth) >= alt_freq_threshold: # Confident SNV
                    pos_within_orig_nuc_sequence = chrom_pos_to_pos_within_orig_nuc_sequence(pos)
                    snv_type, codon_change, aa_change = classify_snv(orig_nuc_sequence, pos_within_orig_nuc_sequence, 
                                                                     comp_dict[alt_allele], verbose=True)
                    if country == 'Mozambique':
                        # print(comp_dict[REF], comp_dict[alt_allele], snv_type, codon_change, aa_change, moz_chrom_pos_anno_dict[gene_chrom][pos]['EFF'])
                        # print(FITEMS['AD'])
                        pass
                    
                    country_snvs_dict[country][snv_type].add((pos, alt_allele))

# In[21]:


country_dnds_dict = {}

for country in country_snvs_dict:
    n_count = 0; s_count = 0
    for variant_type in country_snvs_dict[country]:
        if variant_type in ['nonsynonymous', 'nonsense']:
            n_count += len(country_snvs_dict[country][variant_type])
        elif variant_type == 'synonymous':
            s_count += len(country_snvs_dict[country][variant_type])
    
    print('\t'.join([str(val) for val in [n_count, s_count, country]]))
    dnds = np.float64(n_count/n_denom)/(s_count/s_denom)
    country_dnds_dict[country] = dnds

# In[22]:


print_str = ''
countries_ordered_by_dnds = []
for country, ratio in sorted(country_dnds_dict.items(), key=lambda x: x[1]):
    print_str += ('%s (%.02f)' % (country, ratio)) + ', '
    countries_ordered_by_dnds.append(country)
print(print_str)

# In[26]:


# Compare gene dN/dS and HIV prevalence

fig, ax = plt.subplots()

dnds_list = []
hiv_prev_list = []
for country in country_dnds_dict:
    
    if country == 'Colombia': # Skip Colombia
        continue
    
    true_country = country
    if country == 'Viet Nam':
        true_country = 'Vietnam'
    elif country == 'Gambia':
        true_country = 'Gambia, The'
    elif country == 'Congo DR':
        true_country = 'Congo, Democratic Republic of'
    elif country == 'Ivory Coast':
        true_country = "Cote d'Ivoire"
    elif country in ['Lab', 'Bangladesh']:
        continue
    
    dnds = country_dnds_dict[country]
    hiv_prev = country_hiv_prev_dict[true_country]
    dnds_list.append(dnds)
    hiv_prev_list.append(hiv_prev)
    ax.text(hiv_prev+0.2, dnds-0.01, country)

m,b = np.polyfit(hiv_prev_list, dnds_list, 1)
xs = np.arange(min(hiv_prev_list), max(hiv_prev_list)+1)
ax.plot(xs, m*xs + b, '-', zorder=-1, color='gray')
r = np.corrcoef(hiv_prev_list, dnds_list)[0,1]
# ax.text(10, 2, "r = %.02f" % r, color='gray')

ax.plot(hiv_prev_list, dnds_list, '.')
ax.set_xlabel("HIV prevalence (%)")
ax.set_ylabel("AMA1 dN/dS")
plt.show()

# In[239]:


# Threshold

alt_depth_threshold = 50
alt_freq_threshold = 0.5
min_samples_threshold = 5

# In[13]:


# Determine dN/dS per country
# Impose requirement that at least 5 samples have SNV

# country -> snv_type -> (position, alt allele) -> number of samples
country_snvs_scount_dict = {country: {vtype: defaultdict(int) for vtype in ['synonymous', 'nonsynonymous', 'nonsense']} for country in country_region_dict}

for sample in sample_chrom_pos_record_dict:
    country = sample_country_dict[sample]
    for pos in sample_chrom_pos_record_dict[sample][gene_chrom]:
        
        if pos < gene_trans_start or pos > gene_trans_end: # Skip UTRs
            continue
        
        REF, ALT, QUALITY, FITEMS = sample_chrom_pos_record_dict[sample][gene_chrom][pos]
        
        if FITEMS['DP'] == '.' or QUALITY < 500: # Bad coverage
            continue
        
        alt_alleles = ALT.split(',')
        alt_allele_depths = [int(val) for val in FITEMS['AD'].split(',')[1:]]
        total_depth = int(FITEMS['DP'])
        
        for alt_allele, alt_depth in zip(alt_alleles, alt_allele_depths):
            if REF in ['A', 'G', 'C', 'T'] and alt_allele in ['A', 'G', 'C', 'T']: # SNV
                if alt_depth >= alt_depth_threshold and (alt_depth/total_depth) >= alt_freq_threshold: # Confident SNV
                    pos_within_orig_nuc_sequence = chrom_pos_to_pos_within_orig_nuc_sequence(pos)
                    snv_type, codon_change, aa_change = classify_snv(orig_nuc_sequence, pos_within_orig_nuc_sequence, 
                                                                     comp_dict[alt_allele], verbose=True)
                    country_snvs_scount_dict[country][snv_type][(pos, alt_allele)] += 1
                    if country == 'Mozambique' and snv_type == 'synonymous' and pos == 221977:
                        print(sample)

# In[15]:


FITEMS

# In[244]:


country_snvs_scount_dict['Mozambique']

# In[241]:


country_dnds_dict = {}

for country in country_snvs_scount_dict:
    n_count = 0.1; s_count = 0.1
    for variant_type in country_snvs_scount_dict[country]:
        count_snvs = 0
        for snv in country_snvs_scount_dict[country][variant_type]:
            if country_snvs_scount_dict[country][variant_type][snv] >= min_samples_threshold:
                count_snvs += 1
        
        if variant_type in ['nonsynonymous', 'nonsense']:
            n_count += count_snvs
        elif variant_type == 'synonymous':
            s_count += count_snvs
    
    print('\t'.join([str(val) for val in [n_count, s_count, country]]))
    dnds = np.float64(n_count/n_denom)/(s_count/s_denom)
    country_dnds_dict[country] = dnds

# In[242]:


print_str = ''
countries_ordered_by_dnds = []
for country, ratio in sorted(country_dnds_dict.items(), key=lambda x: x[1]):
    print_str += ('%s (%.02f)' % (country, ratio)) + ', '
    countries_ordered_by_dnds.append(country)
print(print_str)

# In[243]:


# Compare CSP dN/dS and HIV prevalence

fig, ax = plt.subplots()

dnds_list = []
hiv_prev_list = []
for country in country_dnds_dict:
    
    if country == 'Colombia': # Skip Colombia
        continue
    
    true_country = country
    if country == 'Viet Nam':
        true_country = 'Vietnam'
    elif country == 'Gambia':
        true_country = 'Gambia, The'
    elif country == 'Congo DR':
        true_country = 'Congo, Democratic Republic of'
    elif country == 'Ivory Coast':
        true_country = "Cote d'Ivoire"
    elif country in ['Lab', 'Bangladesh']:
        continue
    dnds = country_dnds_dict[country]
    hiv_prev = country_hiv_prev_dict[true_country]
    dnds_list.append(dnds)
    hiv_prev_list.append(hiv_prev)
    ax.text(hiv_prev+0.2, dnds-0.05, country)

m,b = np.polyfit(hiv_prev_list, dnds_list, 1)
xs = np.arange(min(hiv_prev_list), max(hiv_prev_list)+1)
ax.plot(xs, m*xs + b, '-', zorder=-1, color='gray')
r = np.corrcoef(hiv_prev_list, dnds_list)[0,1]
ax.text(10, 0.5, "r = %.02f" % r, color='gray')

ax.plot(hiv_prev_list, dnds_list, '.')
ax.set_xlabel("HIV prevalence (%)")
ax.set_ylabel("CSP dN/dS")
plt.show()
