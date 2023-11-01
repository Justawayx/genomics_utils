import os, gzip, argparse
from collections import defaultdict
import numpy as np

parser = argparse.ArgumentParser(description='Creates an .html file visualizing FASTQ reads mapping to a given insert exactly (heuristic approach)')
parser.add_argument('-f', '--fastq', dest='fastq_paths', help='Path of fastq file(s) to examine', nargs='+')
parser.add_argument('-i', '--insert', dest='expected_insert', help='Path of text file containing DNA sequence of expected insert', required=True)
parser.add_argument('-o', '--output', dest='output_filename', help='Name of output file (default: reads_mapping_to_insert)', required=False)
parser.add_argument('-l', dest='l', type=int, default=50, help='Length of l-mers to split insert into for read filtering (default: 20)', required=False)
parser.add_argument('-k', dest='k', type=int, default=50, help='Length of k-mers used for "alignment" (default: 50)', required=False)

args = parser.parse_args()
fastq_paths = args.fastq_paths
expected_insert = args.expected_insert
output_filename = args.output_filename if args.output_filename is not None else "reads_mapping_to_insert"
l = args.l
k = args.k

html_head = '''
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PvSey1</title>
		<style>
		pre {margin: 0;}
		</style>
  </head>
  <body>
    <main>
'''
html_end = '''    <main>
  </body>
</html>'''

def hamming_distance(kmer1, kmer2):
    distance = 0
    for c1, c2 in zip(kmer1, kmer2):
        distance += (1 if c1 != c2 else 0)
    return distance

expected_insert_split_lmers = set()
for i in range(0, len(expected_insert), l):
    expected_insert_split_lmers.add(expected_insert[i:i+l])

expected_insert_kmer_idx_dict = {}
expected_insert_kmers = set()
for i in range(len(expected_insert)-k+1):
    expected_insert_kmers.add(expected_insert[i:i+k])
    expected_insert_kmer_idx_dict[expected_insert[i:i+k]] = i

#print(f"Using l={l}, k={k}")

candidate_reads = []

for fastq_path in fastq_paths:
	if fastq_path.endswith('.gz'):
		f = gzip.open(fastq_path, 'rt')
	else:
		f = open(fastq_path, 'rt')
	
	#print(f"Searching {fastq_path}...")
	for line in tqdm(f):
			read = f.readline().rstrip('\n'); f.readline(); f.readline()
			for lmer in expected_insert_split_lmers:
					if lmer in read:
							candidate_reads.append(read)
							break

#print(f"{len(candidate_reads)} reads found")

read_expected_insert_start_idx_dict = {}
read_start_idx_dict = {}
read_end_idx_dict = {}

for read in candidate_reads:
    start_match_flag = False
    for i in range(len(read)-k+1):
        read_kmer = read[i:i+k]
        if start_match_flag is True:
            new_expected_nucleotides = expected_insert[expected_insert_start_idx+k:expected_insert_start_idx+k+5]
            new_nucleotides = read[i+k-1:i+k-1+5]
            if expected_insert_start_idx+k > len(expected_insert):
                read_end_idx_dict[read] = i+k-2
                break    
            if hamming_distance(new_expected_nucleotides, new_nucleotides) <= 1:
                expected_insert_start_idx += 1
                continue
            else:
                read_end_idx_dict[read] = i+k-1+5-2
                break
        elif start_match_flag is False and read_kmer in expected_insert_kmers:

            expected_insert_start_idx = expected_insert_kmer_idx_dict[read_kmer]
            read_start_idx = i

            read_expected_insert_start_idx_dict[read] = expected_insert_start_idx
            read_start_idx_dict[read] = read_start_idx

            start_match_flag = True

# TODO
plasmid_anno_dict = {
    "Kpha-chr4-1589887": (45, 327),
    "ori": (389, 977),
    "AmpR": (1148, 2008),
    "HIS4": (2108, 4642),
    "GAP promoter": (5058, 5534),
    "*PvSEY1": (5540, 8208),
    "myc": (8209, 8238),
    "avitag": (8239, 8283),
}

expected_insert_anno = ''
for i in range(len(expected_insert)):
    in_anno = False
    for anno in plasmid_anno_dict:
        start, end = plasmid_anno_dict[anno]
        anno_abbr = anno[0]
        if i >= (start-1) and i < end:
            expected_insert_anno += anno_abbr
            in_anno = True
            break
    if in_anno is False:
        expected_insert_anno += ' '

#o = open(f"{test_dir}/PvSey1_reads_mapping_to_insert.html", 'w')

o.write(html_head)

BUFFER_LEN = 100
buffer = ' '*BUFFER_LEN

o.write('<div style="position:-webkit-sticky;position:sticky;top:0;">' + \
        '<code><pre><span style="background:white;">' + buffer + expected_insert_anno + '</span></code></pre>\n' + \
        '<code><pre><span style="background:yellow;">' +  buffer + expected_insert + '</span></code></pre>' + \
        '</div>' + '\n')

for read, expected_insert_start_idx in sorted(read_expected_insert_start_idx_dict.items(), key=lambda x: x[1]):
    if 'TGTCTCTT' in read: # Sus read
        continue
    read_start_idx = read_start_idx_dict[read]
    read_end_idx = read_end_idx_dict[read] if read in read_end_idx_dict else len(read)
    if expected_insert_start_idx == 0:
        o.write('<code><pre>' + ' '*(BUFFER_LEN-read_start_idx) + \
                '<span style="color:red;">' + read[:read_start_idx] + '</span>' + \
                read[read_start_idx:read_end_idx] + '<span style="color:red;">' + \
                read[read_end_idx:] + '</span>''</code></pre>' + '\n')
    else:
        o.write('<code><pre>' + buffer + ' '*(expected_insert_start_idx) + read[read_start_idx:read_end_idx] + \
                '<span style="color:red;">' + read[read_end_idx:] + '</span>' + '</code></pre>' + '\n')
    i += 1

o.write(html_end)
o.close()