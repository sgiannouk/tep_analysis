###Stavros Giannoukakos### 
import fnmatch
from Bio import SeqIO
import subprocess, pysam
from pprint import pprint
from operator import itemgetter
import random, shutil, time, glob, sys, os, re

thrds = str(60)
# inputFolder = "/home/stavros/playground/test_folder"
inputFolder =  "/home/stavros/shared/elba/"


refTranscriptome_v29 = "/home/stavros/playground/progs/reference_files/reference_transcriptome/ref_transcriptome_v29"

refTranscriptome_v29_fasta = "/home/stavros/playground/progs/reference_files/reference_transcriptome/gencode.v29.pc_transcripts.fa"
# refTranscriptome_v29_fasta = "/home/stavros/playground/progs/reference_files/reference_transcriptome/gencode.v29.1000.fa"

healthy_individuals = "/home/stavros/elba/preprocessed_files/data_info.txt"

# refGenome_GRCh38 = "/home/stavros/playground/progs/reference_files/reference_genome/GRCh38_primAssembly"
reference_annotation = "/home/stavros/playground/progs/reference_files/gene_annotation/gencode.v29.primary_assembly.annotation.gtf"
customTranscriptome = os.path.join(os.path.dirname(refTranscriptome_v29_fasta), "custom_reference_transcriptome.fasta")
resmRef = os.path.join(os.path.dirname(refTranscriptome_v29_fasta), "RESM_customReference")
kallistoRef = os.path.join(os.path.dirname(refTranscriptome_v29_fasta), "kallisto_customReference")

""" All the necessary folders that will host the analysis are being created. 
These include 'preprocessed_files' that will host the filtered and quality 
controlled data. """
# Main folder hosting the analysis
analysisDir = os.path.join(inputFolder, "analysis")

# Main subfolder 
preprocessedFiles = os.path.join(analysisDir, "preprocessed_files")  # Save processed fastq files
alignments = os.path.join(analysisDir, "alignments")  # Saving all alignments 
reportsFolder = os.path.join(analysisDir, "reports")  # Reports folders

# Secondary subfolders
refTranscriptome = os.path.join(alignments, "reference_transcriptome")
expressionAnalysisDir = os.path.join(alignments, "expression_analysis")
## Reports
preprocessingReports = os.path.join(reportsFolder, "preprocessing_reports")
trimmingReports = os.path.join(preprocessingReports, "trimming_reports")
alignmentReports = os.path.join(reportsFolder, "alignment_reports")

temp = os.path.join(preprocessingReports, "temp")

# Generation of the folders
for files in [preprocessedFiles, preprocessingReports, alignmentReports, trimmingReports, temp]:
	if not os.path.exists(files): os.makedirs(files)

""" Function that calls 'TrimGalore' with default parameters. 'MultiQC' will also collect and summarise all 'FastQC' results in
once final report called 'summarised_report'. One can get an overview of the data before the alignment starts. """
def preprocessing_samples():
	""" Run Trimgalore with default settings. """
	print("Preprocessing the input data...")
	raw_data = [f for ext in ["*.fastq.gz", "*.fq.gz"] for f in glob.glob(os.path.join(inputFolder, ext))]
	for i, files in enumerate(raw_data):
		print("{0}/{1}. Preprocessing {2}".format((i+1), len(raw_data), os.path.basename(files)))
		trimGalore = ' '.join([
		"trim_galore",  # Call Trimgalore to preprocess the raw data
		"--output_dir", preprocessedFiles,  # Output directory of processed files
		" --fastqc_args",  # Input of additional FastQC arguments
		"\"--threads", thrds,  # Number of threads to use
		"--outdir", reportsFolder,"\"",
		files])  # Output directory of FastQC reports
		subprocess.run(trimGalore, shell=True) 
	
	os.system('mv %s/*_trimmed_fastqc* %s' %(reportsFolder, preprocessingReports))  # Moving (after trimming) QC reports to qc_reports folder
	os.system('mv %s/*_trimming_report.txt %s' %(preprocessedFiles, preprocessingReports))  # Moving trimming reports to qc_reports folder for MultiQC
	
	""" Run MultiQC to summarise the QC reports from all samples into a summary report """
	runMultiQC = " ".join([
	"multiqc",  # Call MultiQC
	"--quiet",  # Print only log warnings
	"--outdir", reportsFolder,  # Create report in the FastQC reports directory
	"--filename", "summarised_report",  # Name of the output report 
	preprocessingReports])  # Directory where all FastQC and Cutadapt reports reside
	subprocess.run(runMultiQC, shell=True)
	
	os.system('mv %s/*_trimming_report.txt %s' %(preprocessingReports, trimmingReports))  # Moving trimming reports to trimming_reports folder
	os.system('mv %s/*fastqc.zip %s' %(preprocessingReports, temp))  # Moving all 'fastqc.zip' temporary files to temp folder
	os.system('mv  %s/summarised_report_data %s' %(reportsFolder, temp))  # Moving MultiQC temporary files to temp folder
	os.system("chmod 755 -R {0}".format(preprocessedFiles))
	return

""" To obtain the most abundant transcripts, 'Bowtie2' is aligning 20 random healthy individuals from the preprocessed data against the 
reference transcriptome (GENCODE v29). The output SAM files are getting piped into 'Samtools view' to be converted to sorted BAM format. 
'Samtools merge' is then merging all the generated BAM files into one file. From this file, the list of the most abundant transcripts of 
each gene is being  generated. The reference transcriptome and the list of most abundant transcripts per gene are then being used as a 
guide to construct a new custom reference transcriptome. 'Bowtie2' is then used to align all samples against the custom transcriptome. """
def generate_abundant_transcripts(healthy_individuals):
	# Creating the directory that will host the analysis
	if not os.path.exists(refTranscriptome): os.makedirs(refTranscriptome)
	""" Picking up 20 healthy individuals completely randomly. """
	random_healthyIndividuals = []  # List that will host the random healthy individuals
	# Reading the project's sample status to obtain the healthy individuals
	with open(healthy_individuals) as fin:
		for line in fin:
			if not line.startswith("!") and line.split()[1] == "Healthy_Control":
				if str(line.split()[0]+"_trimmed.fq.gz") in os.listdir(preprocessedFiles):
					random_healthyIndividuals.append(os.path.join(preprocessedFiles, line.split()[0]+"_trimmed.fq.gz"))

	#  Picking up 20 individuals randomly and saving their path to random_healthyIndividuals list
	if len(random_healthyIndividuals) >= 20:
		random_healthyIndividuals = random.sample(random_healthyIndividuals, 20)

	""" Calling Bowtie2 to align the preprocessed reads against the reference 
	transcriptome. Bowtie2 is running with '--very-sensitive-local' parameters. """
	print("Bowtie2 - aligning {0} random healthy individuals against the reference transcriptome...".format(len(random_healthyIndividuals)))
	for i, files in enumerate(random_healthyIndividuals):
		bowtie2 = ' '.join([
		"bowtie2",  # Calling 'Bowtie2' to map against the reference transcriptome (GENCODE v29)
		"--threads", thrds,  # Number of threads to be used in Bowtie2 
		"--met-file", os.path.join(alignmentReports, os.path.basename(files).replace("_trimmed.fq.gz", "_createReferenceTranscriptome_report.txt")),  #  Export metrics to reports folder
		"--very-sensitive-local",  # Parameters of alignment
		"-x", refTranscriptome_v29,  # Input the transcriptome's indexes
		"-U", files,  # Input query reads for mapping
  		"|", "samtools view",  # Calling 'samtools view' to compress the Bowtie's output file to BAM
  		"--threads", thrds,  # Number of threads to be used by Samtools in the conversion of the SAM files to BAM
  		"-S -u1",  # Input format is auto-detected, the output file should be in BAM format, and use fast compression
  		"-",  # Piping the input file
  		"|", "samtools sort",  # Calling 'samtools sort' to sort the output alignment file
  		"--threads", thrds,  # Number of threads to be used by 'samtools sort'
  		"-o", os.path.join(refTranscriptome, os.path.basename(files).replace("_trimmed.fq.gz", "_trAligned.bam")), "-"])  # Sorted output  BAM file
		# subprocess.run(bowtie2, shell=True)

	""" Generating a list with the most abundant. """
	# Merging all the BAM files generated from the alignment against the reference transcriptome
	print("Merging the output BAM files...")
	aligned_files = " ".join(glob.glob("{0}/*_trAligned.bam".format(refTranscriptome)))
	samtools_merge = " ".join([
	"samtools merge",  # Calling 'samtools merge' to merge the BAM files - coming from both reference genome and transcriptome
	"--threads", thrds,  # Number of additional threads to use
	os.path.join(refTranscriptome ,"total_transcripts.bam"), # Output file that is located in the alignments/merged_files directory
	aligned_files]) # List of the input BAM files
	# subprocess.run(samtools_merge, shell=True)

	""" Obtaining all transcripts along with their counts in order to eventually get the most abundant from each gene. """
	print("Obtaining the most abundant transcript per gene...")
	transcripts_dict = {}
	bamfile = pysam.AlignmentFile(os.path.join(refTranscriptome ,"total_transcripts.bam"))  # Reading the merged BAM file
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
	
	# for gene, transc_n_count in transcripts_dict.items():
	# 	print(gene, transc_n_count)

	#  Creating a dictionary from the reference transcriptome 
	fasta_transcriptome = SeqIO.to_dict(SeqIO.parse(refTranscriptome_v29_fasta, "fasta"))
	# saving the most abundant transcripts to a new file called "custom_reference_transcriptome.fasta"
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
	
	#Cleaning refTranscriptome folder from pre-aligned files
	## shutil.rmtree(refTranscriptome)

	""" Creating Bowtie2's indexes for the custom reference transcriptome. """
	print("Building Bowtie2's indexes for the custom reference transcriptome...")
	bowtie2_build = ' '.join([
	"bowtie2-build",  # Calling 'bowtie2-build' function to construct the necessary for the alignments indexes
	customTranscriptome,  # Input the custom fasta file
	"{0}/custom_reference_transcriptome".format(os.path.dirname(customTranscriptome))])  # Prefix of the output indexes
	subprocess.run(bowtie2_build, shell=True)

	time.sleep(60)  # Wait for 1' just to make sure that the indexes have been generated

	""" Calling Bowtie2 to align the preprocessed reads against the custom reference 
	transcriptome. Bowtie2 is running with '--very-sensitive-local' parameters. """
	preprocesed_data = glob.glob("{0}/*trimmed.fq.gz".format(preprocessedFiles))
	print("Bowtie2 - aligning all {0} samples against the custom reference transcriptome...".format(len(preprocesed_data)))
	for i, files in enumerate(preprocesed_data):
		bowtie2_align = ' '.join([
		"bowtie2",  # Calling 'Bowtie2' to map against the custom reference transcriptome (generated based on GENCODE v29)
		"--threads", thrds,  # Number of threads to be used in Bowtie2 
		"--un-gz", os.path.join(refTranscriptome, os.path.basename(files).replace("trimmed", "trUnaligned")),  # Write unpaired reads that didn't align in fq.gz format.
		"--met-file", os.path.join(alignmentReports, os.path.basename(files).replace("_trimmed.fq.gz", "_custTranscriptome_report.txt")),  #  Export metrics to reports folder
		"--very-sensitive-local",  # Parameters of alignment
		"-x", customTranscriptome[:-6],  # Input the transcriptome's indexes
		"-U", files,  # Input query reads for mapping
  		"|", "samtools view",  # Calling 'samtools view' to compress the Bowtie's output file to BAM
  		"--threads", thrds,  # Number of threads to be used by Samtools in the conversion of the SAM files to BAM
  		"-S -u1",  # Input format is auto-detected, the output file should be in BAM format, and use fast compression
  		"-",  # Piping the input file
  		"|", "samtools sort",  # Calling 'samtools sort' to sort the output alignment file
  		"--threads", thrds,  # Number of threads to be used by 'samtools sort'
  		"-o", os.path.join(refTranscriptome, os.path.basename(files).replace("_trimmed.fq.gz", "_trAligned1.bam")), "-"])  # Sorted output  BAM file
		subprocess.run(bowtie2_align, shell=True)
	return

""" RSEM is designed to work with reads aligned to transcript sequences, as opposed to whole genome sequences.
	script handles both the alignment of reads against reference transcript sequences and the calculation of relative abundances. """
def expression_analysis():
	# Creating the directory that will host the analysis
	if not os.path.exists(expressionAnalysisDir): os.makedirs(expressionAnalysisDir)
	aligned_files = glob.glob("{0}/*_cutrAligned.bam".format(refTranscriptome))[:2] # Obtaining the aligned files
	
	# Create the necessary indexes for Kallisto quantification  
	if not os.path.exists(kallistoRef):
		subprocess.run("kallisto index --index={0} {1} | tee --append {2}" \
					   .format(kallistoRef ,customTranscriptome,os.path.join(alignmentReports, "generateRefKallisto_report.txt")), shell=True)

	for file in aligned_files:
		if not os.path.exists(file.replace("_cutrAligned.bam", ".fastq.gz")):
			# Using samtools bam2fq tool to convert the aligned BAM files back to fastq which later on will be feed to kallisto
			print(">> Converting {0} in fastQ format...".format(os.path.basename(file)))
			subprocess.run("samtools bam2fq --threads {0} -c 2 -0 {1} {2}".format(thrds, file.replace("_cutrAligned.bam", ".fastq.gz"), file), shell=True)

		# Running Kallisto to obtain the counts per transcript and per sample
		fileDir = os.path.join(expressionAnalysisDir, os.path.basename(file)[:-16])
		if not os.path.exists(fileDir): 
			os.makedirs(fileDir)
			print(">> Analysing {0}".format(os.path.basename(file)))
			kallisto = ' '.join([
			"kallisto quant",  # Calling kallisto quant pipeline
			"--single",  # Quantify single-end reads
			"--plaintext",  # Output in plain text instead of HDF5
			"--fragment-length=200.0",  # Estimated average fragment length
			"--sd=20.0",  # Estimated standard deviation of fragment length
			"--index={0}".format(kallistoRef),  # kallisto index to be used for quantification
			"--threads={0}".format(thrds),  # Number of threads to be used
			"--output-dir={0}".format(fileDir),  # Directory hosting the output files
			file.replace("_cutrAligned.bam", ".fastq.gz"),  # The input fastq file that will be analysed
			"|", "tee", "--append", os.path.join(alignmentReports, "kallisto_report.txt")])  
			# subprocess.run(kallisto, shell=True)

	ref_annot = {}
	with open(reference_annotation) as rin:
		for i, line in enumerate(rin, 1):
			if not line.startswith("#") and "transcript_id" in line:
				gene_id = re.sub('"', '', line.strip().split(";")[0].split("gene_id")[1].strip())
				transcript_id = re.sub('"', '', line.strip().split(";")[1].split("transcript_id")[1].strip())
				transcript_name = re.sub('"', '', line.strip().split(";")[5].split("transcript_name")[1].strip())
				if not "{0}|{1}".format(transcript_id, gene_id) in ref_annot:
					ref_annot["{0}|{1}".format(transcript_id, gene_id)] = transcript_name
				else:
					if not transcript_name in ref_annot["{0}|{1}".format(transcript_id, gene_id)]:
						ref_annot["{0}|{1}".format(transcript_id, gene_id)] = transcript_name
	
	# Exporting useful info regarding the relation between ids and names of the ref. transcriptome 
	with open(os.path.join(expressionAnalysisDir ,"refAnnot_Transcriptome_info.csv"), "w") as refannot_out:
		refannot_out.write("#transcript_id|gene_id,transcript_name\n")
		for ids, names in ref_annot.items():
			refannot_out.write("{0},{1}\n".format(ids, names))
	
	sample_names = []
	expression_matrix = {}
	transcriptome_stats = {}
	# Generating the expression matrix
	expression_files = glob.glob("{0}/*/abundance.tsv".format(expressionAnalysisDir))
	sample_names.append("#transcript_id|gene_id")
	for i, file in enumerate(expression_files, 1):
		sample_names.append(file.split("/")[-2])
		with open(file) as fin:
				for line in fin:
					if not line.startswith("target_id"):
						if i == 1:
							transcriptome_stats["|".join(line.strip().split()[0].split("|")[:2])] = [line.strip().split()[1], line.strip().split()[2]]
							expression_matrix["|".join(line.strip().split()[0].split("|")[:2])] = [line.strip().split()[3]]
						else:
							if ("|".join(line.strip().split()[0].split("|")[:2])) in expression_matrix:
								expression_matrix["|".join(line.strip().split()[0].split("|")[:2])].append(line.strip().split()[3])
							else:
								print("NA")

	# Exporting basic statistics regarding the reference transcriptome 
	with open(os.path.join(expressionAnalysisDir ,"refTranscriptome_stats.csv"), "w") as stats_out:
		stats_out.write("#transcript_id|gene_id,length_of_transcript,effective_length_of_transcript\n")
		for ids, metrics in transcriptome_stats.items():
			stats_out.write("{0},{1},{2}\n".format(ids, metrics[0], metrics[1]))
	# Exporting the expression matrix 
	with open(os.path.join(expressionAnalysisDir ,"expression_matrix.csv"), "w") as exprmat_out:
		exprmat_out.write("{0}\n".format(",".join(sample_names)))
		for ids, counts in expression_matrix.items():
			exprmat_out.write("{0},{1}\n".format(ids,",".join(counts)))
	# Exporting the expression matrix but "transcript_id|gene_id" has been replaced with common "transcript_name"
	with open(os.path.join(expressionAnalysisDir ,"expression_matrix_transcript_names.csv"), "w") as exprmatnames_out:
		sample_names[0] = "#transcript_name"
		exprmatnames_out.write("{0}\n".format(",".join(sample_names)))
		for ids, counts in expression_matrix.items():
			exprmatnames_out.write("{0},{1}\n".format(ref_annot[ids],",".join(counts)))

	return

def variant_detection():

	return

def RNA_editing():

	return

def alternative_splicing():

	return

def novel_transcripts():

	return

def main():
	
	# Performing preprocessing of the data if it isn't already being done
	# if len(glob.glob("{0}/*_trimmed.fq.gz".format(preprocessedFiles))) == 0:
	# 	preprocessing_samples()  # Preprocessing the raw reads

	# generate_abundant_transcripts(healthy_individuals)

	expression_analysis()

if __name__ == "__main__": main()