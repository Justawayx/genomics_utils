samples = []
with open('../sample_list.tsv', 'r') as f:
    for line in f:
        samples.append(line.strip('\n'))

for i in range(1, 14+1):
    chromosome = 'chr%i' % i
    with open('gvcf_%s_list.tsv' % chromosome, 'w') as f:
        for sample in samples:
            f.write('\t'.join([sample, '/opt/input/gVCF_dir/%s/%s.%s.g.vcf' % (chromosome, sample, chromosome)]) + '\n')
