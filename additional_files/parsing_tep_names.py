import os, csv

stats = {}
SraRunTable = "Best_data_info.txt"

with open(SraRunTable) as srin:
	for i, line in enumerate(srin):
		if not line.startswith("!"): 
			gsm_ID = line.strip().split("\t")[0]
			sample_title = line.strip().split("\t")[1]
			sample_description = line.strip().split("\t")[3]
			sra_ID = line.strip().split("\t")[6]
			if sample_description == 'nonCancer':
				stats[gsm_ID] = [(sra_ID, 'Non-Cancer', 'Non-Cancer', sample_title)]
			elif sample_description == 'NSCLC':
				stats[gsm_ID] = [(sra_ID, 'NSCLC', sample_description, sample_title)]
			else:		
				stats[gsm_ID] = [(sra_ID, 'Non-Cancer', sample_description, sample_title)]

stats = dict(sorted(stats.items()))
# for a, b in stats.items():
# 	print(a, '\t'.join(*b))

with open("data_info_new.txt", "w") as outfile:
	outfile.write('!GEO_Accession\tSRA_Accession\tClassification_Group\tSample_Characteristics\tSample_Title\n')
	for a , b in stats.items():
		outfile.write('{0}\t{1}\n'.format(a, '\t'.join(*b)))

