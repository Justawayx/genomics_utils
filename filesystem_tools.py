import sys, os, gzip, argparse
from collections import defaultdict

parser = argparse.ArgumentParser(
	description="Various interactive command line tools for dealing with files.",
	epilog="Available tools: rename (Renames list of files)")
parser.add_argument('action')

args = parser.parse_args()
action = args.action

VALID_ACTIONS = ["rename"]

if action not in VALID_ACTIONS:
	sys.exit("'%s' is not a valid tool" % action)

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
