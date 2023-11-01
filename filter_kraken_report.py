import os
from collections import defaultdict

output_dir = '/oasis/tscc/scratch/dwc001/ghddi_dump_2023'
files = os.listdir(output_dir)

for file in sorted(files):
	if file.endswith("kraken.report"):
		sample = file.split("_kraken.report")[0]
		print(sample)
		
		f = open(f"{output_dir}/{file}", 'r')
		genus_abundance_dict = {}
		genus_species_dict = defaultdict(list)
		for line in f:
			abundance, frags1, frags2, rank_code, NCBI_ID, sci_name = line.strip('\n').split('\t')
			abundance = float(abundance)
			if abundance != 0 and rank_code == "G" or rank_code == 'U':
				genus_abundance_dict[sci_name.strip()] = abundance
			
			if abundance != 0 and rank_code == "S" or rank_code == 'U':
				genus = sci_name.strip().split(' ')[0]
				genus_species_dict[genus].append(sci_name.strip())
		
		low_ab_genus = []
		for genus, abundance in sorted(genus_abundance_dict.items(), key=lambda x: x[1], reverse=True):
				if abundance < 1:
					low_ab_genus.append(genus)
					continue
				species_str = "" if genus in ["Plasmodium", "unclassified"] else f" ({', '.join(genus_species_dict[genus])})"
				print(str(abundance) + '\t' + genus + species_str)
		
		print("")