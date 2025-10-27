# ====================================
# helper_utils.py
# Miscellaneous helper functions
# and small-size reference data
# ====================================

from collections import defaultdict
from urllib.parse import unquote
import os, sys, pickle, copy, gzip
import numpy as np, pandas as pd
from textwrap import wrap

# ============================
# String display
# ============================

# +
# From https://stackoverflow.com/questions/8924173/how-can-i-print-bold-text-in-python
class COLOR:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'

def get_color_list():
    print([name for name in dir(COLOR) if not name.startswith('_')])

def colored_character(char, color):
    return color + char + color + COLOR.END

nuc_color_dict = {'A': COLOR.RED, 'T': COLOR.PURPLE, 'C': COLOR.CYAN, 'G': COLOR.GREEN}

def format_DNA_string(seq):
    string = ''
    for i in range(len(seq)):
        nuc = seq[i]
        string += (nuc_color_dict[nuc] + nuc + nuc_color_dict[nuc])
        string += COLOR.END
    return string

def color_parts_of_string(orig_string, start_end_color_list):
    string = ''
    start_end_color_list = sorted(start_end_color_list) # Sort by start position
    start_end_color_list.append((len(orig_string), len(orig_string), '')) # Last dummy
    j = 0 # Index in start_end_color_list
    for i in range(len(orig_string)):
        start_i, end_i, color = start_end_color_list[j]
        next_start_i, next_end_i, next_color = start_end_color_list[j+1]
        
        if i < start_i: # Have not yet reached current interval
            string += orig_string[i]
        elif i >= start_i and i <= end_i: # In the current interval
            string += colored_character(orig_string[i], color)
        elif i > end_i: # After current interval
            if i >= next_start_i: # In next interval
                string += colored_character(orig_string[i], next_color)
                j += 1
            else: # Not in next interval
                string += orig_string[i]
    return string


# -

def print_blast_like_alignment(alignment, width=60):
    """Pretty-print alignment in BLAST-like format with word wrapping."""
    seq1, identity, seq2, _ = str(alignment).split('\n')
    wrapped_seq1 = wrap(seq1, width, break_on_hyphens=False)
    wrapped_seq2 = wrap(seq2, width, break_on_hyphens=False)
    wrapped_id = wrap(identity, width, break_on_hyphens=False)
    print(len(wrapped_seq1), len(wrapped_seq2), len(wrapped_id))
    
    # Print the alignment in blocks
    print(f"Score: {alignment.score:.1f}" + '\n')
    for block in range(len(wrapped_seq1)):
        start = block * width
        end = start + len(wrapped_seq1[block])
        print(f"Query  {start+1:<5} {wrapped_seq1[block]} {end}")
        print(f"     {' ' * 7} {wrapped_id[block]}")
        print(f"Sbjct  {start+1:<5} {wrapped_seq2[block]} {end}" + '\n')


# ============================
# Text parsing
# ============================

def get_format_data_dict(FORMAT, data):
    return {key: val for key, val in zip(FORMAT.split(':'), data.split(':'))}

def transcript_to_gene_ID(transcript_ID):
    return transcript_ID.split('.')[0]

def parse_info(INFO):
    items = INFO.split(';')
    info_dict = {}
    for item in items:
        try:
            key, value = item.split('=')
            info_dict[key] = value
        except:
            info_dict[item] = True
    return info_dict

def AD_to_AAF(AD):
    depths = np.array([int(d) for d in AD.split(',')])
    return sum(depths[1:])/sum(depths)

# +
def get_AD_for_single_alt_allele(ALT, orig_AD, desired_alt_allele):
    depths = [int(d) for d in orig_AD.split(',')]; alt_depths = depths[1:]
    for depth, alt_allele in zip(alt_depths, ALT.split(',')):
        if alt_allele == desired_alt_allele:
            return '%i,%i' % (depths[0], depth)

def get_AAF_for_single_alt_allele(ALT, orig_AD, desired_alt_allele):
    depths = [int(d) for d in orig_AD.split(',')]; alt_depths = depths[1:]
    desired_alt_allele_idx = ALT.split(',').index(desired_alt_allele)
    return alt_depths[desired_alt_allele_idx]/float(sum(depths))


# -

def get_major_alt_allele(AD, ALT):
    alt_depths = [int(d) for d in AD.split(',')][1:]
    alt_alleles = ALT.split(',')
    allele_depth_tups = sorted(zip(alt_alleles, alt_depths))
    major_alt_allele, major_alt_allele_depth = allele_depth_tups[-1]
    return major_alt_allele, alt_alleles.index(major_alt_allele)

def reverse_complement(dna):
    rc_dna = ''
    comp_nuc_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    for nuc in dna[::-1]:
        rc_dna += comp_nuc_dict[nuc]
    return rc_dna

# ============================
# Gene, effect querying
# ============================

HYPERVARIABLE_KEYWORDS = ['rifin','stevor','PfEMP1','SURFIN']

def desc_has_keyword(description, keywords):
    for keyword in keywords:
        if keyword in description:
            return True
    return False

def is_hypervariable(description):
    return desc_has_keyword(description, HYPERVARIABLE_KEYWORDS)


BAD_EFFECTS = ['intergenic_region', 'synonymous_variant', 'upstream_gene_variant', 'downstream_gene_variant']

gene_conversion_dict = {
    'mal_mito_1': 'PF3D7_MIT02100',
    'mal_mito_2': 'PF3D7_MIT02200',
    'mal_mito_3': 'PF3D7_MIT02300',
}

gene_reverse_conversion_dict = {g2: g1 for g1, g2 in gene_conversion_dict.items()}

# ============================
# Amino acid, mutation nomenclature
# ============================

ordered_aa1s = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

aa3_to_aa1_dict = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Glu": "E", "Gln": "Q", "Gly": "G", 
                   "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", 
                   "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "Ter"}

aa1_to_aa3_dict = {aa1: aa3 for aa3, aa1 in aa3_to_aa1_dict.items()}

AA3_to_aa1_dict = {aa3.upper(): aa1 for aa3, aa1 in aa3_to_aa1_dict.items()}

def aa_change_to_abbr(aa_change):
    if 'p.' not in aa_change:
        return aa_change
    aa_change_only = aa_change.split('/')[0].split('p.')[1]
    aa1_1 = aa3_to_aa1_dict[aa_change_only[:3]]
    if '_' in aa_change_only:
        return aa_change_only
    elif aa_change_only[-3:] in aa3_to_aa1_dict:
        aa1_2 = aa3_to_aa1_dict[aa_change_only[-3:]]
        return aa1_1 + aa_change_only[3:-3] + aa1_2
    else:
        return aa1_1 + aa_change_only[3:]

def parse_aa_change(aa_change):
    if 'p.' not in aa_change:
        return aa_change
    aa_change_only = aa_change.split('/')[0].split('p.')[1]
    
    for char_idx, char in enumerate(aa_change_only):
        if char in '1234567890':
            start_of_residue_number_idx = char_idx
            break
    
    for char_idx in range(start_of_residue_number_idx, len(aa_change_only)):
        if not aa_change_only[char_idx] in '1234567890':
            start_of_second_part_idx = char_idx
            break
    
    return (aa_change_only[:start_of_residue_number_idx], 
            int(aa_change_only[start_of_residue_number_idx:start_of_second_part_idx]),
            aa_change_only[start_of_second_part_idx:])

# ============================
# Sequence properties
# ============================

def calculate_repetitiveness(seq, k, measure='D_R2'):
    if measure == 'D_R2':
        N = len(seq)-k+1 # Number of k-mers
        kmer_count_dict = defaultdict(int)
        for i in range(N):
            kmer = seq[i:i+k]
            kmer_count_dict[kmer] += 1
        background = (N*((1/4)**3))**2
        numerator = 0
        for kmer in kmer_count_dict:
            X = kmer_count_dict[kmer]
            numerator += (X*(X-1) - background)
        return numerator / (N*(N-1))

def calculate_char_frequency(seq):
    char_count_dict = defaultdict(int)
    for char in seq:
        char_count_dict[char] += 1
    n = len(seq)
    char_freq_dict = {char: char_count_dict[char]/n for char in char_count_dict}
    return char_freq_dict

# +
def get_A_run_length(chrom_seq_dict, chrom, pos):
    left_length = 0; right_length = 0
    seq = chrom_seq_dict[chrom]
    orig_pos_idx = pos - 1 # Correct for 0-index
    if seq[orig_pos_idx] != 'A':
        return (left_length, right_length)
    pos_idx = orig_pos_idx
    while seq[pos_idx] == 'A' and pos_idx < len(seq):
        pos_idx += 1
        right_length += 1
    pos_idx = orig_pos_idx
    while seq[pos_idx] == 'A' and pos_idx >= 0:
        pos_idx -= 1
        left_length += 1
    return (left_length, right_length)

# Another idea: distribution of allele frequency in window
# Because the nucleotide where variant call is may not be part of the homopolymer directly
# Or it may there may be sporadic alternate nucleotides
# Or, use frequency of inconsistent indels (suggestive of alignment errors)
# Longest homopolymer run that starts, ends, or spans a 10 bp window

def get_homopolymer_run_length(chrom_seq_dict, chrom, pos):
    left_length = 0; right_length = 0
    seq = chrom_seq_dict[chrom]
    orig_pos_idx = pos - 1 # Correct for 0-index
    orig_nuc = seq[orig_pos_idx]
    pos_idx = orig_pos_idx
    pos_idx += 1
    while seq[pos_idx] == orig_nuc and pos_idx < len(seq):
        pos_idx += 1
        right_length += 1
    pos_idx = orig_pos_idx
    pos_idx -= 1
    while seq[pos_idx] == orig_nuc and pos_idx >= 0:
        pos_idx -= 1
        left_length += 1
    return (left_length, right_length, left_length+right_length+1)

def get_window(chrom_seq_dict, chrom, pos, window_size=50):
    seq = chrom_seq_dict[chrom]
    pos_idx = pos - 1
    half_window_size = window_size//2
    start = max(0, pos_idx-half_window_size)
    end = min(pos_idx+half_window_size, len(seq))
    return chrom_seq_dict[chrom][start:end]

def count_loci_in_window(focal_chrom, focal_pos, chrom_pos_list, window_size=50):
    half_window_size = window_size//2
    start = focal_pos-half_window_size
    end = focal_pos+half_window_size
    count = 0
    for chrom, pos in chrom_pos_list:
        if chrom == focal_chrom and pos >= start and pos <= end:
            count += 1
    return count


# -

# ============================
# Sequence comparison
# ============================

def hamming_distance(kmer1, kmer2):
    distance = 0
    for c1, c2 in zip(kmer1, kmer2):
        distance += (1 if c1 != c2 else 0)
    return distance

# ============================
# Statistics, scoring
# ============================

from scipy.stats import binomtest, fisher_exact, barnard_exact, boschloo_exact, chi2_contingency

def compare_child_to_parent(child_AD, parent_AD):
    pADs = [int(d) for d in parent_AD.split(',')]; pADs = [pADs[0], sum(pADs[1:])]
    cADs = [int(d) for d in child_AD.split(',')]; cADs = [cADs[0], sum(cADs[1:])]
    table = np.array([pADs, cADs]).T
    if len(pADs) == 2:
        res = fisher_exact(table, alternative='two-sided')
        return res.pvalue
    else:
        print('sus')
        res = chi2_contingency(table)
        return res.pvalue

def convert_value_to_score(value, min_value, max_value, min_score=0, max_score=1):
    if value < min_value:
        return min_score
    elif value > max_value:
        return max_score
    else:
        return min_score + ((value-min_value)/(max_value-min_value))*(max_score-min_score)

# ============================
# Sequence alignment
# ============================

def GlobalAlignment(s1, s2, MATCH_SCORE = 0, MISMATCH_SCORE = -1, GAP_SCORE = -1):
	score = [[0 for _ in range(len(s2)+1)] for _ in range(len(s1)+1)]
	backtrack = [['' for _ in range(len(s2)+1)] for _ in range(len(s1)+1)]
	
	for i in range(1, len(s1)+1):
		score[i][0] = score[i-1][0] + GAP_SCORE
		backtrack[i][0] = '|'
	
	for j in range(1, len(s2)+1):
		score[0][j] = score[0][j-1] + GAP_SCORE
		backtrack[0][j] = '-'
	
	for i in range(1, len(s1)+1):
		for j in range(1, len(s2)+1):
			down_score = score[i-1][j] + GAP_SCORE # Gap in s2
			right_score = score[i][j-1] + GAP_SCORE # Gap in s1
			diag_score = score[i-1][j-1] + (MATCH_SCORE if s1[i-1] == s2[j-1] else MISMATCH_SCORE)
			
			score[i][j] = max(down_score, right_score, diag_score)
			
			if score[i][j] == down_score:
				backtrack[i][j] = '|'
			elif score[i][j] == right_score:
				backtrack[i][j] = '-'
			elif score[i][j] == diag_score:
				backtrack[i][j] = '/'
	
	s1_aligned = ''; s2_aligned = ''
	i = len(s1); j = len(s2)
	
	while i != 0 or j != 0:
		if backtrack[i][j] == '|':
			s1_aligned += s1[i-1]
			s2_aligned += '-'
			i -= 1
		elif backtrack[i][j] == '-':
			s1_aligned += '-'
			s2_aligned += s2[j-1]
			j -= 1
		elif backtrack[i][j] == '/':
			s1_aligned += s1[i-1]
			s2_aligned += s2[j-1]
			i -= 1; j -= 1
	
	s1_aligned, s2_aligned = s1_aligned[::-1], s2_aligned[::-1]
	
	return score[len(s1)][len(s2)], s1_aligned, s2_aligned

# ============================
# BLOSUM and PAM matrices
# ============================

PAM250_dict = defaultdict(dict)
PAM250_matrix = np.zeros((20, 20))

f = open('/storage/NFS/ROTATION_PROJECT/daisy/REF_DATA/pam250.txt')
header = f.readline()
col_aas = header.strip().split()
for line in f:
    items = line.strip().split()
    row_aa = items[0]
    for col_aa, value in zip(col_aas, items[1:]):
        PAM250_dict[row_aa][col_aa] = int(value)

for idx1, aa_1 in enumerate(ordered_aa1s):
    for idx2, aa_2 in enumerate(ordered_aa1s):
        PAM250_matrix[idx1][idx2] = PAM250_dict[aa_1][aa_2]

BLOSUM62_dict = defaultdict(dict)
BLOSUM62_matrix = np.zeros((20, 20))

f = open('/storage/NFS/ROTATION_PROJECT/daisy/REF_DATA/blosum62.txt')
header = f.readline()
col_aas = header.strip().split()
for line in f:
    items = line.strip().split()
    row_aa = items[0]
    for col_aa, value in zip(col_aas, items[1:]):
        BLOSUM62_dict[row_aa][col_aa] = int(value)

for idx1, aa_1 in enumerate(ordered_aa1s):
    for idx2, aa_2 in enumerate(ordered_aa1s):
        BLOSUM62_matrix[idx1][idx2] = BLOSUM62_dict[aa_1][aa_2]

# The more positive, the more conserved
print("BLOSUM62 min and max:", np.min(BLOSUM62_matrix), np.max(BLOSUM62_matrix))
print("PAM250 min and max:", np.min(PAM250_matrix), np.max(PAM250_matrix))
