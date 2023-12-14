import sys, os, gzip, argparse
from collections import defaultdict

parser = argparse.ArgumentParser(
	description="Various interactive command line tools for dealing with files.",
	epilog="Available tools: rename (Renames list of files)")
parser.add_argument('action')

args = parser.parse_args()
action = args.action

VALID_ACTIONS = ["rename", "combine_vcf"]

if action not in VALID_ACTIONS:
	sys.exit("'%s' is not a valid tool" % action)

if action.lower() == "combine_vcf":
	
	file_dir = input("Enter directory in which files are located: ").strip()
	
	if not os.path.exists(file_dir): # Check if directory exists
		sys.exit("ERROR: Invalid directory %s" % file_dir)
	
	vcfs_path = input("Enter path to a text file that contains ordered list of VCFs to combine, first to last: ")
	
	if not os.path.exists(vcfs_path):
		sys.exit("ERROR: VCF list %s does not exist" % file_path)
	
	ordered_vcfs = [line.strip() for line in open(vcfs_path, 'r')]
	
	ordered_site_intervals = []
	
	for vcf in ordered_vcfs:
		vcf_path = '%s/%s' % (file_dir, vcf)
		if not os.path.exists(vcf_path):
			sys.exit("ERROR: %s does not exist" % vcf_path)
		
		f = open(vcf_path, 'r')
		first_site = None; last_site = None
		for line in f:
			if not line.startswith('#'):
				CHROM, POS = line.split('\t')[:2]
				first_site = (CHROM, POS)
				break
		
		for line in f:
			CHROM, POS = line.split('\t')[:2]
			last_site = (CHROM, POS)
		
		ordered_site_intervals.append((first_site, last_site))
	
	print("Combining")
	for first_site, last_site in ordered_site_intervals:
		print('\t'.join(first_site) +'\t' + '\t'.join(last_site))
	
	answer = input("Proceed? y/n")
	while answer.lower() not in ['y', 'n']:
		answer = input("Please enter y or n: ")
	
	if answer.lower() == 'n':
		sys.exit("Aborting")
	
	combined_vcf_path = '%s/%s_COMBINED' % (file_path, vcfs[0])
	o = open(combined_vcf_path, 'w')
	
	for i, vcf in enumerate(vcfs):
		vcf_path = '%s/%s' % (file_dir, vcf)
		first_site, last_site = ordered_site_intervals[i]
		if i == (len(ordered_site_intervals)-1):
			next_start_chrom, next_start_pos = ordered_site_intervals[i+1][0]
			f = open(vcf_path, 'r')
			for line in f:
				if i == 0:
					o.write(line + '\n')
				if line.startswith('#CHROM'):
					break
			
			for line in f:
				CHROM, POS = line.split('\t')[:2]
				if CHROM == next_start_chrom and POS == next_start_pos:
					break
				o.write(line + '\n')
			
		else:
			f = open(vcf_path, 'r')
			for line in f:
				if line.startswith('#CHROM'):
					break
			
			for line in f:
				o.write(line + '\n')
	
	print("Done! Check %s" % combined_vcf_path)
    

if action.lower() == "rename": # Rename
	
	file_dir = input("Enter directory in which files are located: ").strip()
	
	if not os.path.exists(file_dir): # Check if directory exists
		sys.exit("ERROR: Invalid directory %s" % file_dir)
	
	map_path = input("Enter path to a text file that contains two tab separated columns, left column with the ORIGINAL file name, right column with the NEW file name: ")
	
	if not os.path.exists(map_path):
		sys.exit("ERROR: Mapping file %s does not exist" % file_path)
	
	orig_fnames = []
	new_fnames = []
	f = open(map_path, 'r')
	for line in f:
		if line.rstrip('\n') == '':
			continue
		try:
			orig_fname, new_fname = line.rstrip('\n').split('\t')
			orig_fnames.append(orig_fname); new_fnames.append(new_fname)
		except:
			sys.exit("ERROR: Mapping file is formatted improperly")
	f.close()
	
	if len(orig_fnames) != len(set(orig_fnames)):
		sys.exit("ERROR: Original file name list has duplicates")
	
	orig_fpaths = []
	for orig_fname in orig_fnames:
		file_path = "%s/%s" % (file_dir, orig_fname)		
		if not os.path.exists(file_path):
			sys.exit("ERROR: File %s does not exist" % file_path)
		orig_fpaths.append(file_path)
	
	if len(new_fnames) != len(set(new_fnames)):
		sys.exit("ERROR: New file name list has duplicates")
	elif len(new_fnames) != len(orig_fnames):
		sys.exit("ERROR: Original and new filename list lengths don't match")
	
	new_fpaths = []
	for new_fname in new_fnames:
		file_path = "%s/%s" % (file_dir, new_fname)		
		if os.path.exists(file_path):
			sys.exit("ERROR: File %s exists, would be overwritten" % file_path)
		new_fpaths.append(file_path)
	
	command_list = []
	for orig_fpath, new_fpath in zip(orig_fpaths, new_fpaths):
		command = "mv %s %s" % (orig_fpath, new_fpath)
		command_list.append(command)
		os.system(command)
	
	with open("rename_log.txt", 'w') as f:
		for command in command_list:
			f.write(command + '\n')
	
	print("Done! Commands run are in rename_log.txt")
