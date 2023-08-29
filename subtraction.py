# Load gene annotations

gene_desc_dict = {}

f = open('/projects/winzeler/GENOME_RESOURCES/sc/S288C_reference_genome_R64-1-1_20110203/gene_association_R64-1-1_20110205.sgd', 'r')

for line in f:
	if line[0] == '!':
		continue
	items = line.strip().split('\t')
	gene = items[2]
	desc = items[9]
	gene_desc_dict[gene] = desc

'''
o = open('genedescript.txt', 'w')
f = open('hold14.txt', 'r')
for line in f:
	items = line.strip().split(' ')
	gene_item = items[10]
	if gene_item != '@':
		gene = gene_item[1:]
		parts = gene.split('_')
		if len(parts) == 3:
			gene = parts[0] + '(' + parts[1] + ')' + parts[2]
		if gene not in gene_desc_dict:
			desc = ''
			print("Couldn't find %s" % gene)
		else:
			desc = gene_desc_dict[gene]
	else:
		desc = ''
	o.write('@' + desc + '\n')
'''

f = open('2268_combined_snv_final_GATK-Filters.tsv', 'r')
with open('2268_combined_snv_final_GATK-Filters.new.tsv', 'w') as o:
	header = f.readline()
	o.write(header)
	for line in f:
		items = line.strip().split('\t')
		gene = items[2]
		parts = gene.split('_')
		if len(parts) == 3:
			gene = parts[0] + '(' + parts[1] + ')' + parts[2]
		
		if gene not in gene_desc_dict:
			desc = ''
			print("Couldn't find %s" % gene)
		else:
			items[3] = gene_desc_dict[gene]
		
		o.write('\t'.join(items) + '\n')
	
	f.close()