import sys, os, re

answer = input("Has fastq_dir in config.cfg been updated? y/n")
if answer.lower != 'y':
	sys.exit("Make sure to update fastq_dir in config.cfg")

def get_env_value_from_config(env_variable):
	CONFIG_PATH = os.path.abspath('config.cfg')
	for line in open(CONFIG_PATH, 'r'):
		if line.startswith('#'):
			continue
		if '=' in line:
			key, value = line.rstrip('\n').split('=')
			if key == env_variable:
				return value
	print("Environmental variable %s not found in config.cfg" % env_variable)
	return False

fastq_dir = get_env_value_from_config('fastq_dir')
species_abbr = get_env_value_from_config('species_abbr')

answer = input("Do the fastq files in fastq_dir need to be renamed (and are from IGM)? y/n")
if answer.lower != 'y':
	
	# Current IGM format:
	# {sample}_Sx_Lxxx_Rx_001.fastq.gz
	
	rename_commands = []
	samples = []
	
	for orig_fastq_name in os.listdir(fastq_dir):
		if 'fastq' in orig_fastq_name:
			try:
				sample, suffix = re.split("_S[0-9]+_L[0-9]+_", orig_fastq_name)
			except:
				print(f"Skipping {orig_fastq_name}")
				continue
			samples.append(sample)
			new_fastq_name = f"{species_abbr}_{sample}_{suffix[:2]}.fastq.gz"
			rename_commands.append(f"mv {fastq_dir}/{orig_fastq_name} {fastq_dir}/{new_fastq_name}")
	
	if len(samples) > 0:
		f = open("samples.txt", 'w')
		f.write('\n'.join(sorted(set(samples))))
		f.close()
	
	for command in rename_commands:
		print(command)
		os.system(command)
	
	print("---------------------------------------------")
	print(f"{len(rename_commands)} fastq files have been renamed in {fastq_dir}")
	print(f"{len(set(samples))} samples are included in samples.txt")
	print("Double check samples.txt before proceeding")
	print("---------------------------------------------")

else:	
	print("Make samples.txt in this directory (does not include species prefix or read suffix)")