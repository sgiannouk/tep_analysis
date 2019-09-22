###Stavros Giannoukakos### 

#Version of the program
__version__ = "0.1.0"

import argparse
from Bio import SeqIO
import subprocess, pysam
import multiprocessing as mp
from operator import itemgetter
from multiprocessing import Pool
import random, shutil, time, glob, csv, sys, os, re


# tep_data =  "/shared/projects/tom_teps_umc"
tep_data = "/home/stavros/playground/tep_analysis/test_data"

healthy_individuals = "/home/stavros/elba/preprocessed_files/data_info.txt"
refTranscGRCh38 = "/home/stavros/references/reference_transcriptome/GRCh38_gencode.v31.transcripts.fa.gz"
refAnnot = "/home/stavros/references/reference_annotation/GRCh38_gencode.v31.primAssembly_psudo_trna.annotation.gtf"
customTranscriptome = os.path.join(os.path.dirname(sys.argv[0]), "/cstm_ref/custom_ref_transcriptome.fasta")


usage = "tep_analysis [options]"
epilog = " -- January 2019 | Stavros Giannoukakos -- "
description = "DESCRIPTION"

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, usage=usage, description=description, epilog=epilog)
# Number of threads/CPUs to be used
parser.add_argument('-t', '--threads', dest='threads', default=str(40), metavar='', 
                	help="Number of threads to be used in the analysis")
# Number of threads/CPUs to be used
parser.add_argument('-p', '--partialAssignment', required=False, metavar='', 
                	help="Multi-mapping reads will partially assigned in the transcripts.\nFor example, if one read is aligned in 2 different transcripts\nit will be counted as 0.5. (Default: full assignemnt)")
# Display the version of the pipeline 
parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(__version__))
# Get the options and return them
args = parser.parse_args()


current_dir = os.getcwd()

# Main folder hosting the analysis
analysis_dir = os.path.join(current_dir, "tep_data_analysis")

# Subfolder 
reports_dir = os.path.join(analysis_dir, "reports")
prepr_dir = os.path.join(analysis_dir, "preprocessed_data")  # Save processed fastq files
preprocessing_reports_dir = os.path.join(reports_dir, "preprocessing_reports")
trimming_reports_dir = os.path.join(preprocessing_reports_dir, "trimming_reports")
alignments_dir = os.path.join(analysis_dir, "alignments") 
alignment_reports_dir = os.path.join(reports_dir, "alignment_reports")
expression_analysis_dir = os.path.join(analysis_dir, "expression_analysis")
var_calling_dir = os.path.join(analysis_dir, "variant_calling")


expression_matrix = {}  # Dictionary hosting all expressions


def preprocessing_samples():
	""" Function that calls 'TrimGalore' with default parameters. 'MultiQC' will also collect 
	and  summarise  all 'FastQC' results in once final report called 'summarised_report'. One 
	can get an overview of the data before the alignment starts. """

	### Run Trimgalore with default settings
	if not os.path.exists(prepr_dir): os.makedirs(prepr_dir)  # Creating necessary directory
	if not os.path.exists(trimming_reports_dir): os.makedirs(trimming_reports_dir)  # Creating necessary directory
	print("Preprocessing the input data...")
	raw_data = [f for f in glob.glob(os.path.join(tep_data, "*.gz"))]
	for i, files in enumerate(raw_data, 1):
		print("{0}/{1}. Preprocessing {2}".format(i, len(raw_data), os.path.basename(files)))
		trimGalore = ' '.join([
		"trim_galore",  # Call Trimgalore to preprocess the raw data
		"--output_dir", prepr_dir,  # Output directory of processed files
		" --fastqc_args",  # Input of additional FastQC arguments
		"\"--threads", args.threads,  # Number of threads to use
		"--outdir", trimming_reports_dir,"\"",
		files,  # Output directory of FastQC reports
		"2>>", os.path.join(trimming_reports_dir, "trimgalore_report.txt")])
		subprocess.run(trimGalore, shell=True) 
	
	### Moving reports to "trimming_reports" directory
	os.system('mv %s/*_trimming_report.txt %s' %(prepr_dir, trimming_reports_dir))  # Moving trimming reports to qc_reports folder for MultiQC
	
	### Run MultiQC to summarise the QC reports from all samples into a summary report
	runMultiQC = " ".join([
	"multiqc",  # Call MultiQC
	"--quiet",  # Print only log warnings
	"--outdir", trimming_reports_dir,  # Create report in the FastQC reports directory
	"--filename", "preliminary_summarised_report",  # Name of the output report 
	trimming_reports_dir,  # Directory where all FastQC and Cutadapt reports reside
	"2>>", os.path.join(trimming_reports_dir, "multiQC-report.txt")])  # Output multiQC report
	subprocess.run(runMultiQC, shell=True)
	
	os.system('mv %s/*summarised_report.html %s' %(trimming_reports_dir, preprocessing_reports_dir))  # Moving trimming reports to trimming_reports folder
	return

def create_cstm_transcriptome(healthy_individuals):
	""" To obtain the most abundant transcripts, 'Bowtie2' is aligning 20 random healthy individuals 
	from the preprocessed data against the reference transcriptome (GENCODE v31). The output SAM files 
	are getting piped into 'Samtools view' to be converted to sorted BAM format. 'Samtools merge' is
	then merging all the generated BAM files into one file. From this file, the list of the most abundant 
	transcripts of each gene is being  generated. The reference transcriptome and the list of most abundant 
	transcripts per gene are then being used as a guide to construct a new custom reference transcriptome. 
	'Bowtie2' is then used to align all samples against the custom transcriptome. """

	### Creating the directory that will host the analysis
	temp = os.path.join(analysis_dir, "temp_directory")
	if not os.path.exists(temp): os.makedirs(temp)
	### Picking up 20 healthy individuals completely randomly.
	random_healthyIndividuals = []  # List that will host the random healthy individuals
	### Reading the project's sample status to obtain the healthy individuals
	with open(healthy_individuals) as fin:
		for line in fin:
			if not line.startswith("!") and line.split()[1] == "Healthy_Control":
				if str(line.split()[0]+"_trimmed.fq.gz") in os.listdir(prepr_dir):
					random_healthyIndividuals.append(os.path.join(prepr_dir, line.split()[0]+"_trimmed.fq.gz"))

	###  Picking up 20 individuals randomly and saving their path to random_healthyIndividuals list
	if len(random_healthyIndividuals) >= 20:
		random_healthyIndividuals = random.sample(random_healthyIndividuals, 20)
	else:
		random_healthyIndividuals = random_healthyIndividuals

	### Calling Bowtie2 to align the preprocessed reads against the reference  transcriptome. Bowtie2 is running with '--very-sensitive-local' parameters
	if not os.path.exists(alignment_reports_dir): os.makedirs(alignment_reports_dir)
	print("Bowtie2 - aligning {0} random healthy individuals against the reference transcriptome...".format(len(random_healthyIndividuals)))
	for files in random_healthyIndividuals:
		bowtie2_cstm = ' '.join([
		"bowtie2",  # Calling 'Bowtie2' to map against the reference transcriptome (GENCODE v31)
		"--threads", args.threads,  # Number of threads to be used in Bowtie2 
		"--very-sensitive-local",  # Parameters of alignment
		"-x", os.path.join(os.path.basename(refTranscGRCh38),"bowtie2_GRCh38_gencode.v31.transcripts/GRCh38_gencode.v31.transcripts"),  # Input the transcriptome's indexes
		"-U", files,  # Input query reads for mapping
  		"|", "samtools sort",  # Calling 'samtools sort' to sort the output alignment file
		"--threads", args.threads,  # Number of threads to be used by 'samtools sort'
		"--output-fmt BAM",  # Output in BAM format
  		"-o", os.path.join(temp, os.path.basename(files).replace("_trimmed.fq.gz", "_trAligned.bam")), "-",  # Sorted output  BAM file
  		"2>>", os.path.join(alignment_reports_dir, "bt2_cstm_ref-report.txt")])  # Directory where bowtie2 reports reside
		subprocess.run(bowtie2_cstm, shell=True)

	### Generating a list with the most abundant
	# Merging all the BAM files generated from the alignment against the reference transcriptome
	aligned_files = " ".join(glob.glob("{0}/*_trAligned.bam".format(temp)))
	print("Merging the output BAM files...")
	samtools_merge = " ".join([
	"samtools merge",  # Calling 'samtools merge' to merge the BAM files - coming from both reference genome and transcriptome
	"--threads", args.threads,  # Number of additional threads to use
	os.path.join(temp ,"total_transcripts.bam"), # Output file that is located in the temp directory
	aligned_files, # List of the input BAM files
	"2>>", os.path.join(alignment_reports_dir, "samtools_merge-report.txt")])  # Samtools merge reports
	subprocess.run(samtools_merge, shell=True)

	### Obtaining all transcripts along with their counts in order to eventually get the most abundant from each gene
	transcripts_dict = {}
	print("Obtaining the most abundant transcript per gene...")
	bamfile = pysam.AlignmentFile(os.path.join(temp ,"total_transcripts.bam"))  # Reading the merged BAM file
	for read in bamfile:
		ref = bamfile.get_reference_name(read.reference_id)  # Obtaining the reference names from each entry
		if ref != None:  # Skipping not mapped reds 
			gene_id = ref.split("|")[1].strip()  # Obtaining the gene id from each header
			transcript_id = ref.split("|")[0].strip()  # Obtaining the transcript id from each header
			# Creating a tuple that saves [gene] = [(transcript1, count), ... , (transcriptN, count)]
			if gene_id in transcripts_dict:
				if [(key, value) for key, value in transcripts_dict[gene_id] if key  == transcript_id]:
					transcripts_dict[gene_id] = ([(key, value+1) for key, value in transcripts_dict[gene_id] if key  == transcript_id])
				else:
					transcripts_dict[gene_id].append((transcript_id, 1))
			else:
				transcripts_dict[gene_id] = [(transcript_id, 1)]
	bamfile.close()  # Closing the BAM file after reading

	###  Creating a dictionary from the reference transcriptome 
	fasta_transcriptome = SeqIO.to_dict(SeqIO.parse(refTranscGRCh38, "fasta"))
	### saving the most abundant transcripts to a new file called "custom_ref_transcriptome.fasta"
	print("Generating the custom reference transcriptome...")
	with open(customTranscriptome, "a") as outRef:
		for gene, transc_n_count in transcripts_dict.items():
			if len(transc_n_count) == 1:
				ID = str(transc_n_count[0][0]+"|"+gene+"|")
				record = [(key, value) for key, value in fasta_transcriptome.items() if ID in key]
				SeqIO.write(record[0][1], outRef, "fasta")
			else:
				most_abundant = max(transc_n_count, key = itemgetter(1))[0]
				ID = str(most_abundant+"|"+gene+"|")
				record = [(key, value) for key, value in fasta_transcriptome.items() if ID in key]
				SeqIO.write(record[0][1], outRef, "fasta")

	transcripts_dict.clear()  # Emptying transcripts_dict dictionary
	fasta_transcriptome.clear()  # Emptying fasta_transcriptome dictionary
	
	### Cleaning refTranscriptome folder from pre-aligned files
	shutil.rmtree(temp)
	return 

def aligning_cstm_transcritome():
	""" Creating Bowtie2's indexes for the custom reference transcriptome. """
	if not os.path.exists("{0}.1.bt2".format(customTranscriptome[:-6])):
		print("Building Bowtie2's indexes for the custom reference transcriptome...")
		bowtie2_build = ' '.join([
		"bowtie2-build",  # Calling 'bowtie2-build' function to construct the necessary for the alignments_dir indexes
		"--threads", args.threads,  # Number of threads to be used by Bowtie2-built
		customTranscriptome,  # Input the custom fasta file
		customTranscriptome[:-6],  # Prefix of the output indexes
		"2>>", os.path.join(alignment_reports_dir, "bt2-built_cstm_ref-report.txt")])  # Directory where bowtie2 reports reside
		subprocess.run(bowtie2_build, shell=True)
		time.sleep(60)  # Wait for 1' just to make sure that the indexes have been generated

	""" Calling Bowtie2 to align the preprocessed reads against the custom reference 
	transcriptome. Bowtie2 is running with '--very-sensitive-local' parameters. """
	preprocesed_data = glob.glob("{0}/*trimmed.fq.gz".format(prepr_dir))
	if not os.path.exists(alignments_dir): os.makedirs(alignments_dir)
	print("Bowtie2 - aligning in total {0} samples against the custom reference transcriptome...".format(len(preprocesed_data)))
	for files in preprocesed_data:
		bowtie2_align = ' '.join([
		"bowtie2",  # Calling 'Bowtie2' to map against the custom reference transcriptome (generated based on GENCODE v29)
		"--threads", args.threads,  # Number of threads to be used in Bowtie2 
		"--all",  # report all alignments_dir; very slow, MAPQ not meaningful
		"--un-gz", os.path.join(alignments_dir, os.path.basename(files).replace("_trimmed", ".unaligned")),  # Write unpaired reads that didn't align in fq.gz format.
		"--sensitive-local",  # Default parameters of local alignment
		"-x", customTranscriptome[:-6],  # Input the transcriptome's indexes
		"-U", files,  # Input query reads for mapping
  		"|", "samtools sort",  # Calling 'samtools sort' to sort the output alignment file
		"--threads", args.threads,  # Number of threads to be used by 'samtools sort'
		"--output-fmt BAM",  # Output in BAM format
  		"-o", os.path.join(alignments_dir, os.path.basename(files).replace("_trimmed.fq.gz", ".aligned.bam")), "-"])  # Sorted output  BAM file
		subprocess.run(bowtie2_align, shell=True)
	### Generating basic mapping stats
	mapping_quality_control()
	return

def mapping_quality_control():
	### EXPORTING ALIGNMENT STATS
	print(">>> RSeQC and Picard - Generating alignment stats: in progress ..")
	aligned_data = glob.glob("{0}/*.aligned.bam".format(alignments_dir))
	for file in aligned_data:
		file_name = os.path.basename(file).split(".")[0]
		
		samtools_index = " ".join([
		"samtools index",  # Indexing the concat_samples.bam file
		"-@", args.threads,  # Number of threads to be used
		file])  # Input BAM file
		subprocess.run(samtools_index, shell=True)

		# BAM stats
		bam_stat = ' '.join([
		"bam_stat.py",
		"-i", file,  # Input BAM file
		"> {0}/{1}.bamstat.txt".format(alignment_reports_dir, file_name),  # Output file
		"2>>", os.path.join(postanalysis_dir, "bamstats-report.txt")])
		subprocess.run(bam_stat, shell=True)

		# Picard CollectAlignmentSummaryMetrics
		CollectAlignmentSummaryMetrics = ' '.join([
		"picard CollectAlignmentSummaryMetrics",  # Call picard CollectAlignmentSummaryMetrics
		"INPUT= {0}".format(file),  # Input BAM file
		"OUTPUT= {0}/{1}.{2}.alignment_metrics.txt".format(postanalysis_dir, file_name, aligned_to),  # Output
		"REFERENCE_SEQUENCE= {0}".format(customTranscriptome),  # Reference sequence file
		"2>>", os.path.join(alignment_reports_dir, "CollectAlignmentSummaryMetrics-report.txt")])
		subprocess.run(CollectAlignmentSummaryMetrics, shell=True)

	print(">>> multiQC - Summarising all QC reports: in progress ..")
	multiQC = " ".join([
	"multiqc",  # Call MultiQC
	"--quiet",  # Print only log warnings
	"--outdir", alignment_reports_dir,  # Create report in the FastQC reports directory
	"--filename", "post-alignment_summarised_report",  # Name of the output report 
	alignment_reports_dir,  # Directory where all FastQC and Cutadapt reports reside
	"2>>", os.path.join(alignment_reports_dir, "mapping_multiQC-report.txt")])  # Output multiQC report
	subprocess.run(multiQC, shell=True)
	return

def expression_analysis(file):
	""" Summarisation and construction of the expression matrix on transcript level. For this task each 
	alignment file (BAM) is being read and all the necessary information is being extracted. In the end 
	all info will be summasised in an expression matrix. """

	### Creating the directory that will host the analysis
	if not os.path.exists(expression_analysis_dir): os.makedirs(expression_analysis_dir)

	annotation = []
	### Obtaining all the annotation used for the construction of the reference transcriptome
	with open(customTranscriptome) as fin:
		for line in fin:
			if line.startswith(">"):
				annotation.append(line.strip().strip(">").split("|")[0])

	matrix = [0] * len(annotation)
	bamfile = pysam.AlignmentFile(file, threads=int(args.threads)) # Reading the BAM file
	for read in bamfile:
		if not read.is_reverse:  # Getting only the reads mapping to the forward strand
			ref = bamfile.get_reference_name(read.reference_id)  # Obtaining the reference names from each entry
			if ref != None:  # Skipping not mapped reads 
				transcript_id = ref.split("|")[0].strip()  # Obtaining the transcript id from each header
				if not args.partialAssignment:
					matrix[annotation.index(transcript_id)] += 1				
	bamfile.close()  # Closing the BAM file after reading
	
	### Save all the info in the disctionary
	expression_matrix[os.path.basename(file)[:-12]] = matrix

	### For each read BAM file, then extracted counts will be written in a csv file
	with open("{0}/{1}_expression_matrix.csv".format(expression_analysis_dir, os.path.basename(file)[:-12]), "w") as mat_out:
		writer = csv.writer(mat_out, delimiter='\t')
		writer.writerows(zip(annotation, matrix))
	return

def generate_expression_matirx():
	""" Merging all expression matrices from each sample in order to generate 
	the final expression matrix containing all samples. """
	annot = []
	### Obtaining all info from each file
	expr_files = glob.glob("{0}/*_expression_matrix.csv".format(expression_analysis_dir)) # Obtaining the expression matrices
	for i, file in enumerate(expr_files, 1):
		if i == 1:
			with open(file) as fin:
				for line in fin:
					annot.append(line.strip().split("\t")[0])
	for i, file in enumerate(expr_files, 1):
		reads = []
		with open(file) as fin:
			for line in fin:
				reads.append(line.strip().split("\t")[1])
		expression_matrix[os.path.basename(file.replace("_expression_matrix.csv", ""))] = reads
	### Writing the expression matrix in a csv file
	with open(os.path.join(expression_analysis_dir, "expression_matrix.csv"), "w") as fout:
		fout.write("gene_name\t{0}\n".format("\t".join(expression_matrix.keys())))
		writer = csv.writer(fout, delimiter='\t')
		writer.writerows(zip(annot,*expression_matrix.values()))
	return

def variant_detection():

	if not os.path.exists(var_calling_dir): os.makedirs(var_calling_dir)
	aligned_files = glob.glob("{0}/*.aligned.bam".format(refTranscriptome))[:2]
	
	# Indexing the custom reference transcriptome
	if not os.path.exists("{0}.fai".format(customTranscriptome)):
		subprocess.run("samtools faidx {0}".format(customTranscriptome) ,shell=True)
	
	# variant_calling = ' '.join([
	# "parallel",
	# "-j", args.threads,  # Run n jobs in parallel
	# "\"" ,"bcftools mpileup",
	# "-Ou",  # Output type uncompressed BCF
	# "--fasta-ref", customTranscriptome,  # faidx indexed reference sequence file
	# {},
	# "|", "bcftools call",  # Piping the output to bcftools call
	# "--threads", int(args.threads),  # Number of threads to be used
	# "-m",  # Alternative model for multiallelic and rare-variant calling
	# "-Oz", 
	# "-o", os.path.join(var_calling_dir, {}.replace("cutrAligned.bam" ,"variants.vcf")), "\"",  # Output vcf 
	# ":::", " ".join(aligned_files),
	# "|", "tee", "--append", os.path.join(alignment_reports_dir, "variantCalling_report.txt")])
	# # subprocess.run(variant_calling, shell=True)
	# print(variant_calling)

	# Calling variants in the aligend files
	# for file in aligned_files:
	# 	variant_calling = ' '.join([
	# 	"bcftools mpileup",
	# 	"--threads", int(args.threads),  # Number of threads to be used
	# 	"--fasta-ref", customTranscriptome,  # faidx indexed reference sequence file
	# 	file,  # Input aligned file
	# 	"|", "bcftools call",  # Piping the output to bcftools call
	# 	"--threads", int(args.threads),  # Number of threads to be used
	# 	"-Ov",  # Output type uncompressed VCF
	# 	"-m",  # Alternative model for multiallelic and rare-variant calling
	# 	"-o", os.path.join(var_calling_dir, file.replace("cutrAligned.bam" ,"variants.vcf")),  # Output vcf
	# 	"|", "tee", "--append", os.path.join(alignment_reports_dir, "variantCalling_report.txt")])
	# 	subprocess.run(variant_calling, shell=True)

	for file in aligned_files:
		print(file)
		variant_calling = ' '.join([
		"bcftools mpileup",
		"--threads", args.threads,  # Number of threads to be used
		"--fasta-ref", customTranscriptome,  # faidx indexed reference sequence file
		"-Ov",
		"-o", os.path.join(var_calling_dir, file.replace("cutrAligned.bam" ,"variants.vcf")),  # Output vcf
		file,
		"|", "tee", "--append", os.path.join(alignment_reports_dir, "variantCalling_report.txt")])
		subprocess.run(variant_calling, shell=True)


		# file,  # Input aligned file
		# "|", "bcftools call",  # Piping the output to bcftools call
		# "--threads", int(args.threads),  # Number of threads to be used
		# "-Ov",  # Output type uncompressed VCF
		# "-m",  # Alternative model for multiallelic and rare-variant calling
		# "-o", os.path.join(var_calling_dir, file.replace("cutrAligned.bam" ,"variants.vcf")),  # Output vcf
		# "|", "tee", "--append", os.path.join(alignment_reports_dir, "variantCalling_report.txt")])
		# subprocess.run(variant_calling, shell=True)

	return

def RNA_editing():

	return

def alternative_splicing():

	return

def novel_transcripts():

	return

def clean_up():
	## REMOVING UNNECESSARY FILES & REPORTS (reports directory)
	for path, subdir, folder in os.walk(reports_dir):
		for name in folder:
			file = os.path.join(path, name)
			if os.stat(file).st_size == 0 or\
			(name.endswith("multiQC-report.txt") and os.stat(file).st_size == 583) or\
	  		 name.endswith("fastqc.zip"):
				os.remove(file)

	return

def main():
	
	# Performing preprocessing of the data if it isn't already being done
	if len(glob.glob("{0}/*_trimmed.fq.gz".format(prepr_dir))) == 0: 
		preprocessing_samples()  # Preprocessing the raw reads

	# Generating the custom reference transcriptome, if not already done
	if not os.path.exists(customTranscriptome):
		create_cstm_transcriptome(healthy_individuals)

	# Aligning all samples against the custom ref. transcriptome
	aligning_cstm_transcritome():

	# Generating the expression matrix of all samples aligned against the custom ref. transcriptome
	arg = [file for file in glob.glob("{0}/*.aligned.bam".format(alignments_dir))]
	pool = mp.Pool(processes=int(args.threads))
	pool.map(expression_analysis, arg)
	pool.close()

	generate_expression_matirx()

	# variant_detection()
	# RNA_editing()
	
if __name__ == "__main__": main()
