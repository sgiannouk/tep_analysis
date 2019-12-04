import os, csv
import subprocess

stats = {}
SraRunTable = "SraRunTable.txt"
sample_dir = "/home/stavros/playground/tep_analysis/analysis/preprocessed_data"


with open(SraRunTable) as srin:
	for i, line in enumerate(srin):
		if not line.startswith("AvgSpotLen"): 
			biosample = line.strip().split("\t")[1]  # for verification
			srr_ID = line.strip().split("\t")[7]  # id we have
			gsm_ID = line.strip().split("\t")[9]  # for verification
			sra_ID = line.strip().split("\t")[8]  # for verification ## USE FOR RENAMING

			if (biosample, gsm_ID, sra_ID) in stats:
				stats[(biosample, gsm_ID, sra_ID)].append(srr_ID)
			else:
				stats[(biosample, gsm_ID, sra_ID)] = [(srr_ID)]

i = 0
j = 0 
for aa, bb in stats.items():
	if len(bb) == 1:
		i+=1
		sample = os.path.join(sample_dir, "{0}.trimmed.fq.gz".format(bb[0]))
		if os.path.exists(sample):
			# subprocess.run('mv {0} {1}/{2}.trimmed.fq.gz'.format(sample, sample_dir, aa[1]), shell=True)
			print('mv {0} {1}/{2}.trimmed.fq.gz'.format(sample, sample_dir, aa[1]))
	if len(bb) == 2:
		j+=1	
		first = os.path.join(sample_dir, "{0}.trimmed.fq.gz".format(bb[0]))
		second = os.path.join(sample_dir ,"{0}.trimmed.fq.gz".format(bb[1]))
		if os.path.exists(first) and os.path.exists(second):
			# subprocess.run('cat {0} {1} > {2}/{3}.trimmed.fq.gz'.format(first, second, sample_dir, aa[1]), shell=True)
			print('cat {0} {1} > {2}/{3}.trimmed.fq.gz'.format(first, second, sample_dir, aa[1]))


