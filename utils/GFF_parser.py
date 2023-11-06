# Imports
import gzip, json

# GFF/GTF format information:
# https://useast.ensembl.org/info/website/upload/gff.html

def parse_GFF_attribute(attribute_str):
	pairs = attribute_str.rstrip(';').split('; ')
	attribute_dict = {}
	for pair in pairs:
		key, value = pair.strip().split(' ')
		attribute_dict[key] = json.loads(value)
	return attribute_dict

fname = 'R64-1-1.82_GFF/Saccharomyces_cerevisiae.R64-1-1.82.gtf.gz'

if fname[-3:] == '.gz':
	f = gzip.open(fname, 'rt')
else:
	f = open(fname, 'rt')

gene_IDs = set()

for line in f:
	if line[0] == '#':
		print(line.rstrip('\n'))
		continue
	
	items = line.rstrip('\n').split('\t')
	seqname, source, feature, start, end, score, strand, frame, attribute = items
	start = int(start); end = int(end)
	
	if feature == 'gene':
		attribute_dict = parse_GFF_attribute(attribute)
		gene_ID = attribute_dict['gene_id']
		gene_IDs.add(gene_ID)

# Generate convenience files
with open('gene_list.txt', 'w') as o:
	for gene_ID in sorted(gene_IDs):
		o.write(gene_ID + '\n')

gene_ID_name_dict = {}

with open('summary.csv', 'r') as f:
	for line in f:
		items = line.rstrip('\n').strip('"').split('" "')
		gene_ID = items[0]
		gene_name = items[5]
		gene_ID_name_dict[gene_ID] = gene_name