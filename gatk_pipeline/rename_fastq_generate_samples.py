import os, re

CONFIG_PATH = os.path.abspath('config.cfg')

for line in open(CONFIG_PATH, 'r'):
	if line.startswith("fastq_dir="):
		fastq_dir = line.split("fastq_dir=")[1].strip()
		break

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
		new_fastq_name = f"{sample}_{suffix[:2]}.fastq.gz"
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