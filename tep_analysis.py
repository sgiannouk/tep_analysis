# Stavros Giannoukakos

#Version of the program
__version__ = "0.1.3"

import json, csv
import numpy as np
import pandas as pd
import argparse, subprocess
import multiprocessing as mp
from datetime import datetime
from operator import itemgetter
from multiprocessing import Pool
from natsort import index_natsorted
import random, shutil, glob, sys, os
startTime = datetime.now()



### INPUT DATA
raw_data = "/home/stavros/playground/platelets_analysis/data"
### CLINICAL DATA
clinical_data = f"{os.path.dirname(os.path.realpath(__file__))}/additional_files/clinical_data.tsv"
### ADAPTERS AND PHIX SEQUENCES
adapter_sequences = f"{os.path.dirname(os.path.realpath(__file__))}/additional_files/adapter_sequences.fasta"
phix_sequences = f"{os.path.dirname(os.path.realpath(__file__))}/additional_files/phix.fasta"
### REFERENCE FILES
refGenomeGRCh38 = "/home/stavros/references/reference_genome/GRCh38_GencodeV31_primAssembly/GRCh38.primary_assembly.genome.fa"
refTranscriptomeGRCh38 = "/home/stavros/references/reference_transcriptome/GRCh38_GencodeV35/GRCh38_gencode.v35.transcripts.fa"
star_indexes = "/home/stavros/references/reference_genome/GRCh38_GencodeV31_primAssembly/GRCh38.primary_assembly.genome.STARindexes"
refAnnot = "/home/stavros/references/reference_annotation/GRCh38_gencode.v35.primary_assembly.annotation.gtf"
### FUSION DETECTION
arriba = "/home/stavros/playground/progs/arriba_v2.1.0/./arriba"
arriba_visualisation = "/home/stavros/playground/progs/arriba_v2.1.0/./draw_fusions.R"
cytobandsGRCh38 = "/home/stavros/playground/progs/arriba_v2.1.0/database/cytobands_hg38_GRCh38_v2.1.0.tsv"
blacklistGRCh38 = "/home/stavros/playground/progs/arriba_v2.1.0/database/blacklist_hg38_GRCh38_v2.1.0.tsv.gz"
knownfusionsGRCh38 = "/home/stavros/playground/progs/arriba_v2.1.0/database/known_fusions_hg38_GRCh38_v2.1.0.tsv.gz"
proteindomainsGRCh38 = "/home/stavros/playground/progs/arriba_v2.1.0/database/protein_domains_hg38_GRCh38_v2.1.0.gff3"
### ISOFORM DETECTION
salmon = "/home/stavros/playground/progs/salmon/src/salmon"
ref_salmon = "/home/stavros/references/reference_annotation/GRCh38_gencode.v35.gffread.fa"
### ALTERNATIVE SPLICING
asgal = "python3 /home/stavros/playground/progs/galig/asgal"
### GATK FILES
gatk_known_sites = "/home/stavros/references/gatk_subfiles/known_sites"
### RNA EDITING
reditools = "python3 /home/stavros/playground/progs/reditools2.0/src/cineca/reditools.py"




usage = "platelets_processing.py [options]"
epilog = " -- October 2020 | Stavros Giannoukakos -- "
description = "DESCRIPTION"

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, usage=usage, description=description, epilog=epilog)
# Adapter sequence
parser.add_argument('-a', '--adapter', dest='adapter', default=adapter_sequences, metavar='', 
                	help="Adapter sequence to be trimmed from the input data")
# PhiX sequence
parser.add_argument('-p', '--phix', dest='phix', default=phix_sequences, metavar='', 
                	help="PhiX sequence to be removed from the input data")
# Number of threads/CPUs to be used
parser.add_argument('-t', '--threads', dest='threads', default=str(30), metavar='', 
                	help="Number of threads to be used in the analysis")
# Allow "low" confidence in gene fusions
parser.add_argument('-l', '--low', dest='low', action="store_true", 
                	help="During gene fusion matrix construction,\nallow \"low\" confidence fusions (default False)")
# Display the version of the pipeline
parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(__version__))
# Get the options and return them
args = parser.parse_args()



current_dir = os.path.dirname(os.path.realpath(__file__))
# Main folder hosting the analysis
analysis_dir = os.path.join(current_dir, "analysis")
# Reporting subdirectories
reports_dir = os.path.join(analysis_dir, "reports")
pipeline_reports = os.path.join(analysis_dir, "reports/pipeline_reports")
preprocessing_reports_dir = os.path.join(analysis_dir, "reports/preprocessing_reports")
alignment_reports_dir = os.path.join(analysis_dir, "reports/alignment_reports")
genefusion_reports_dir = os.path.join(analysis_dir, "reports/gene_fusion_reports")
varcall_reports_dir = os.path.join(analysis_dir, "reports/variant_calling_reports")
# Main subdirectories hosting the individual analysis
preprocessing_dir = os.path.join(analysis_dir, "preprocessed_data")
alignments_dir = os.path.join(analysis_dir, "alignments")
isoform_quant_dir = os.path.join(analysis_dir, "isoform_expression")
genefusion_dir = os.path.join(analysis_dir, "gene_fusion")
expression_analysis_dir = os.path.join(analysis_dir, "gene_expression")
alternativesplicing_dir = os.path.join(analysis_dir, "alternative_splicing")
altexpression_dir = os.path.join(analysis_dir, "alt_isoform_expression")
varcall_dir = os.path.join(analysis_dir, "mutation_profiling")
funcannot_dir = os.path.join(analysis_dir, "mutation_profiling/functional_annotation")
rnaediting_dir = os.path.join(analysis_dir, "rna_editing")

temp = os.path.join(varcall_dir, "temp")
if not os.path.exists(pipeline_reports): os.makedirs(pipeline_reports)


def preprocessing_samples(): 
	""" Function that calls BBDuk to remove non-biological sequences from the input data.
	It will also remove low quality reads and short sequences. FastQC will perform a 
	thorough quality control on the preprocessed reads. 'MultiQC' will also collect 
	and  summarise  all the results in one final report (summarised_report'). 
	One can get an overview of the data before the alignment starts. """
	print(f'\t{datetime.now().strftime("%d.%m.%Y %H:%M")} PREPROCESSING THE INPUT SAMPLES')
	
	
	if not os.path.exists(preprocessing_dir): os.makedirs(preprocessing_dir)  # Creating necessary directory
	if not os.path.exists(preprocessing_reports_dir): os.makedirs(preprocessing_reports_dir)  # Creating necessary directory

	os.chdir(preprocessing_dir)
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  BBDuk (v38.22) - Preprocessing the input data...')
	data = glob.glob(os.path.join(raw_data, "*.fastq.gz"))  ### Obtaining the raw data
	### Run BBDuk with custom settings ###
	for i, file in enumerate(data, 1):
		print(f'\n{i}/{len(data)}. Preprocessing {os.path.basename(file)}')
		startTime = datetime.now()  # Recording the start time
		sample_name = os.path.basename(file).split(".")[0]
		runBBDuk = ' '.join([
		"bbduk.sh",  # Call BBDuk (BBTools) to preprocess the raw data
		f"threads={args.threads}",  # Set number of threads to use
		f"in={file}",  # Input of the forward file
		"out={0}".format(os.path.join(preprocessing_dir, f"{sample_name}.trimmed.fq.gz")),  # Export filtered read to file
		"trimpolya=20",  # Trim poly-A or poly-T tails of at least 20nt on either end of reads
		"trimpolyg=20",  # Trim poly-G of at least 20nt on either end of reads
		"entropymask=f",  # Discard low-entropy reads
		"maxns=1",  # Reads with at least 1 N (after trimming) will be discarded
		"ktrim=r",  # Trim reads on the right to remove bases matching reference kmers
		"qtrim=r",  # Trim read right ends to remove bases with quality below trimq
		"trimq=20",  # Regions with average quality BELOW this will be trimmed
		"minlength=40",  # Reads shorter than this after trimming will be discarded
		"minavgquality=20",  # Reads with average quality (after trimming) below this will be discarded (40)
		"ziplevel=7",  #  Compression level 7
		"k=25",  # Setting the kmer size we want to search for
		"mink=8",  # Look for shorter kmers at read tips down to this length
		"hammingdistance=1",  # Maximum Hamming distance for ref kmers
		f"ref={args.adapter},{args.phix}",  # contaminant files
		# f"enthist={preprocessing_reports_dir}/{sample_name}.enthist",  # Read entropy histogram
		f"stats={preprocessing_reports_dir}/{sample_name}.stats",  # Write statistics about which contaminants were detected
		"2>>", os.path.join(pipeline_reports, "preprocessing_bbduk-report.txt")])  # Output report
		subprocess.run(runBBDuk, shell=True)
		print(f'\tPreprocessing finished in {datetime.now() - startTime}')


	# print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  FastQC (v0.11.8) - Quality Control reports are being generated...')
	# runFastQC = ' '.join([
	# "fastqc",  # Call fastQC to quality control all processed data
	# "--threads", args.threads,  # Number of threads to use
	# "--quiet",  # Print only log warnings
	# "--outdir", preprocessing_reports_dir,  # Create all output files in this specified output directory
	# ' '.join(data),  # String containing all samples that are about to be checked
	# "2>>", os.path.join(pipeline_reports, "preprocessing_fastqc-report.txt")])  # Output fastQC report
	# subprocess.run(runFastQC, shell=True)
	
	# ### Run MultiQC to summarise the QC reports from all samples into a summary report
	# print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  MultiQC (v1.8) - Summarising the preprocessed data...')
	# runMultiQC = " ".join([
	# "multiqc",  # Call MultiQC
	# "--quiet",  # Print only log warnings
	# # "--template", "simple",  # Report template to use
	# "--title", "\'Overall preprocessing summary report\'",
	# "--outdir", preprocessing_reports_dir,  # Create report in the FastQC reports directory
	# "--filename", "preliminary_summarised_report",  # Name of the output report 
	# preprocessing_reports_dir,  # Directory where all FastQC and Cutadapt reports reside
	# "2>>", os.path.join(pipeline_reports, "preprocessing_multiqc-report.txt")])  # Output multiQC report
	# subprocess.run(runMultiQC, shell=True)

	# ### Moving multiqc summarised report to report directory and removing fastqc.zip files
	# subprocess.run(f'mv {preprocessing_reports_dir}/preliminary_summarised_report.html {reports_dir}', shell=True)
	# subprocess.run(f'rm {preprocessing_reports_dir}/*.zip', shell=True)
	return

class geneNisoform_level_analysis:
	
	def __init__(self):

		print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} GENE AND ISOFORM QUANTIFICATION')
		

		preprocesed_data = glob.glob(f"{preprocessing_dir}/*trimmed.fq.gz")
		for current_num, sample in enumerate(preprocesed_data, 1):
			sample_name = os.path.basename(sample).split(".")[0]
			print(f'\n{current_num}/{len(preprocesed_data)} | Processing sample {sample_name}...')
			
			# ### Calling the different methods
			# self.aligning_against_refGenome(sample, sample_name)
			# self.detect_gene_fusions(sample_name)
			# self.isoform_quantification(sample_name)
			self.alternative_splicing_events(sample_name)
		return

	def aligning_against_refGenome(self, sample, sample_name):
		""" Calling STAR to align the preprocessed reads against the reference genome """		
		
		if not os.path.exists(alignments_dir): os.makedirs(alignments_dir)
		if not os.path.exists(alignment_reports_dir): os.makedirs(alignment_reports_dir)
		if not os.path.exists(expression_analysis_dir): os.makedirs(expression_analysis_dir)
		os.chdir(alignments_dir)  # Changing working directory cause STAR output several files in the wd
		

		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  STAR (v2.7.6a) - Aligning against the reference genome...')
		startTime = datetime.now()  # Recording the start time
		STAR_align = " ".join([
		"STAR",  # Calling STAR aligner
		"--runMode alignReads",  # Map reads against the reference genome
		"--genomeDir", star_indexes,  # Load the  ref. genome
		"--runThreadN", args.threads,  # Number of threads to be used by STAR
		"--readFilesIn", sample,  # Input read
		"--sjdbOverhang 100",  # Length of the donor/acceptor sequence on each side of the junctions
		"--quantMode GeneCounts TranscriptomeSAM",  # Count reads per gene and output in transcriptome coordinates
		"--quantTranscriptomeBan IndelSoftclipSingleend",  # Prohibit indels, soft clipping and single-end alignments - compatible with RSEM
		"--twopassMode Basic",  # Basic 2-pass mapping, with all 1st pass junctions inserted into the genome indexes on the fly
		"--genomeLoad", "NoSharedMemory",  # Load genome into shared and remove it after run
		"--outSAMtype", "BAM SortedByCoordinate",  # Sort the output BAM file by Coordinate
		"--limitBAMsortRAM", "50000000000",  # Maximum available RAM (bytes) for sorting BAM (50GB)
		"--outSAMunmapped Within",  # Output the unmapped reads in the SAM file
		"--readFilesCommand gunzip -c",  # Unzipping the input .gz files
		"--chimSegmentMin 12",  # Check fusion genes
		"--chimJunctionOverhangMin 10",  # Minimum overhang for a chimeric junction
		"--chimOutJunctionFormat 0",  # Formatting type for the Chimeric.out.junction file
		"--alignSJstitchMismatchNmax 5 -1 5 5",  # Maximum number of mismatches for stitching of the splice junctions
		"--chimScoreDropMax 30",  # Max drop (difference) of chimeric score (the sum of scores of all chimeric segments) from the read length
		"--chimScoreJunctionNonGTAG 0",  # Penalty for a non-GT/AG chimeric junction
		"--chimScoreSeparation 1",  # Minimum difference (separation) between the best chimeric score and the next one
		"--chimSegmentReadGapMax 3",  # Maximum gap in the read sequence between chimeric segments
		"--chimMultimapNmax 50",  # Maximum number of chimeric multi-alignments
		"--chimOutType WithinBAM Junctions",  # Output old SAM into separate Chimeric.out.sam file
		"--outSAMattrRGline", f"ID:{sample_name} SM:{sample_name} PL:Illumina",  # Add SAM read group line
		"--outFileNamePrefix", os.path.join(alignments_dir, os.path.basename(sample)[:-13]),  # Output BAM files
		"--outReadsUnmapped Fastx", os.path.join(alignments_dir, os.path.basename(sample)[:-13]),  # Output of unaligned reads in fastq format
		"2>>", os.path.join(pipeline_reports, "aligning_star-report.txt")])  # Input files
		subprocess.run(STAR_align, shell=True)
		# subprocess.run(f'mv {alignments_dir}/{sample_name}.Chimeric.out.sam {alignments_dir}/{sample_name}.chimeric.sam', shell=True)  # Renaming chimeric read files		
		subprocess.run(f'mv {alignments_dir}/{sample_name}.Aligned.sortedByCoord.out.bam {alignments_dir}/{sample_name}.aligned.genome.bam', shell=True)  # Renaming mapped (gemone) read files
		subprocess.run(f'mv {alignments_dir}/{sample_name}.Aligned.toTranscriptome.out.bam {alignments_dir}/{sample_name}.aligned.transcriptome.bam', shell=True)  # Renaming mapped (gemone) read files		
		subprocess.run(f'mv {alignments_dir}/{sample_name}.ReadsPerGene.out.tab {expression_analysis_dir}/{sample_name}.gene-expression.tsv', shell=True)  # Moving and renaming the expression files		
		subprocess.run(f'mv {alignments_dir}/{sample_name}.Unmapped.out.mate1 {alignments_dir}/{sample_name}.unmapped.fastq', shell=True)  # Renaming unmapped read files
		subprocess.run(f'pigz -9 --fast --processes {args.threads} {alignments_dir}/{sample_name}.unmapped.fastq', shell=True)  # Compress the unmapped reads
		print(f'\tAlignment and Gene Fusion detection finished in {datetime.now() - startTime}')
		
		### Cleaning 'alignments_dir' directory from unnecessary files
		subprocess.run(f'mv {alignments_dir}/*.out {alignment_reports_dir}', shell=True)  # Remove STAR directories
		subprocess.run(f'rm -r {alignments_dir}/*STAR*', shell=True)  # Remove STAR directories
		return

	def detect_gene_fusions(self, sample_name):
		""" Using the Arriba gene fusion pipeline to detect gene fusions 
		from the RNA-Seq data of all the samples. """

		if not os.path.exists(genefusion_dir): os.makedirs(genefusion_dir)
		if not os.path.exists(genefusion_reports_dir): os.makedirs(genefusion_reports_dir)


		print(f'\n{datetime.now().strftime("%d.%m.%Y %H:%M")}  Arriba (v2.0.0) - Detection of gene fusions...')
		startTime = datetime.now()  # Recording the start time
		subprocess.run(f'samtools index -@ {args.threads} {alignments_dir}/{sample_name}.aligned.genome.bam', shell=True)
		runArriba = " ".join([
		arriba,  # Calling Arriba pipeline
		"-S 5",  # The 'min_support' filter discards all fusions with fewer than this many supporting reads
		# "-c", f"{alignments_dir}/{sample_name}.chimeric.sam",  # File in SAM format with chimeric alignments as generated by STAR (Chimeric.out.sam)
		"-x", f"{alignments_dir}/{sample_name}.aligned.genome.bam",  # File in BAM format with main alignments as generated by STAR (Aligned.out.bam)
		"-g", refAnnot,  # GTF file with gene annotation
		"-a", refGenomeGRCh38,  # FastA file with the reference genome sequence (assembly)
		"-b", blacklistGRCh38,  # File containing blacklisted events (recurrent artifacts and transcripts observed in healthy tissue)
		"-k", knownfusionsGRCh38,  # File containing known/recurrent fusions. Boosting sensitivity for entities characterized by fusions between the same pair of genes
		"-p", proteindomainsGRCh38,  # File in GFF3 format containing coordinates of the protein domains of genes
		"-o", f"{genefusion_dir}/{sample_name}.fusion.tsv",  # Output file with fusions that have passed all filters
		"2>>", os.path.join(pipeline_reports, "genefusion_arriba-report.txt")])  # Output report
		subprocess.run(runArriba, shell=True)

		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Arriba Vis (v2.0.0) - Visualisation of detected gene fusions...')
		runArribaVis = " ".join([
		arriba_visualisation,  # Calling Arriba visualisation pipeline
		f"--fusions={genefusion_dir}/{sample_name}.fusion.tsv",  # Input file with all fusions that have passed all filters
		f"--alignments={alignments_dir}/{sample_name}.aligned.genome.bam",  # File in BAM format with main alignments as generated by STAR (Aligned.out.bam)
		f"--output={genefusion_reports_dir}/{sample_name}.fusions.pdf",  # Output as .pdf
		f"--annotation={refAnnot}",  # GTF file with gene annotation
		f"--cytobands={cytobandsGRCh38}",  #
		f"--proteinDomains={proteindomainsGRCh38}",  # File in GFF3 format containing coordinates of the protein domains of genes
		"2>>", os.path.join(pipeline_reports, "genefusion_arribaVis-report.txt")])  # Output report
		# subprocess.run(runArribaVis, shell=True)
		print(f'\tGene fusions detection finished in {datetime.now() - startTime}')
		return

	def isoform_quantification(self, sample_name):
		""" Calling RSEM after STAR to perform isoform estimation """
		
		if not os.path.exists(isoform_quant_dir): os.makedirs(isoform_quant_dir)


		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Salmon (v1.3.0) - Estimating isoform expression...')
		startTime = datetime.now()  # Recording the start time
		run_salmon = " ".join([
		salmon, "quant",  # Calling salmon quant
		"--quiet",  # Be quiet while doing quantification
		"--useVBOpt",  # Use the Variational Bayesian EM
		"--seqBias",  # Perform sequence-specific bias correction
		"--gcBias",  # Perform fragment GC bias correction
		"--gencode",  # The transcript fasta is in GENCODE format
		"--libType U",  # Definition of the strandedness of the RNA-Seq reads
		### Bootstraps are required for estimation of technical variance
		"--numBootstraps 50",  # Number of bootstrap samples to generate
		"--threads", args.threads,  # Number of threads to be used
		"--targets", ref_salmon,  # FASTA format file containing target transcripts
		"--output", f"{isoform_quant_dir}/{sample_name}",  #  Output quantification directory
		"--alignments", f"{alignments_dir}/{sample_name}.aligned.transcriptome.bam",  # Input BAM file
		"2>>", os.path.join(pipeline_reports, "isoformquant_salmon-report.txt")])  # Input files
		subprocess.run(run_salmon, shell=True)
		print(f'\tIsoform expression estimation (salmon) finished in {datetime.now() - startTime}')
		return
		
	def alternative_splicing_events(self, sample_name):
		""" Calling ASGAL (Alternative Splicing Graph ALigner) to detect alternative 
		splicing events in out datasets """
		
		if not os.path.exists(alternativesplicing_dir): os.makedirs(alternativesplicing_dir)


		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  ASGAL - Estimating alternative splicing events...')
		startTime = datetime.now()  # Recording the start time
		# output_dir = os.path.join(altexpression_dir, sample_name)
		# os.makedirs(output_dir)
		runAsgal = " ".join([
		asgal,  # Calling ASGAL pipeline
		"--multi",  # Run genome-wide mode
		"--support 5",  # Minimum intron coverage
		"--threads", args.threads,  # Number of threads to use for salmon mapping and parallel gene computation
		"--genome", refGenomeGRCh38,  # Path to reference genome
		"--annotation", refAnnot,  # Path to reference annotation
		"--transcripts", refTranscriptomeGRCh38,  # Reference transcriptome
		"--sample", f"{preprocessing_dir}/{sample_name}.trimmed.fq.gz",  # Input preprocessed sample
		"--output", alternativesplicing_dir,  # Path of the output directory 
		# "2>>", os.path.join(pipeline_reports, "altsplicing_asgal-report.txt")
		])  # Output report file
		subprocess.run(runAsgal, shell=True)
		print(f'\tAlternative splicing events detection finished in {datetime.now() - startTime}')
		return

def mapping_quality_control():

	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} GENERATING ALIGNMENT STATS')


	### EXPORTING ALIGNMENT STATS
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  RSeQC and Picard - Generating alignment stats...')
	aligned_data = glob.glob(f"{alignments_dir}/*.aligned.genome.bam")
	for file in aligned_data:
		file_name = os.path.basename(file).split(".")[0]
		# RSeQC - BAM stats 
		bam_stat = ' '.join([
		"bam_stat.py",
		"-i", file,  # Input BAM file
		f"> {alignment_reports_dir}/{file_name}.bamstat.txt",  # Output file
		"2>>", os.path.join(pipeline_reports, "postalignment_bamstats-report.txt")])
		subprocess.run(bam_stat, shell=True)

		# Picard CollectMultipleMetrics
		collectMultipleMetrics = ' '.join([
		"picard-tools CollectMultipleMetrics",  # Call picard-tools CollectMultipleMetrics
		f"INPUT= {file}",  # Input BAM file
		f"OUTPUT= {alignment_reports_dir}/{file_name}.",  # Output
		f"REFERENCE_SEQUENCE= {refGenomeGRCh38}",  # Reference sequence file
		"2>>", os.path.join(pipeline_reports, "postalignment_CollectMultipleMetrics-report.txt")])
		subprocess.run(collectMultipleMetrics, shell=True)

	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  multiQC - Summarising all QC reports...')
	multiQC = " ".join([
	"multiqc",  # Call MultiQC
	"--quiet",  # Print only log warnings
	"--outdir", alignment_reports_dir,  # Create report in the FastQC reports directory
	# "--template", "simple",  # Report template to use
	"--title", "\'Overall post-alignment summary report\'",
	"--filename", "post-alignment_summarised_report",  # Name of the output report 
	alignment_reports_dir,  # Directory where all FastQC and Cutadapt reports reside
	"2>>", os.path.join(pipeline_reports, "postalignment_mutiqc-report.txt")])  # Output multiQC report
	subprocess.run(multiQC, shell=True)

	### Moving multiqc summarised report to report directory
	subprocess.run(f'mv {alignment_reports_dir}/*summarised_report.html {reports_dir}', shell=True)
	return

def gene_expression_matrix():
	""" Summarisation and construction of the expression matrix on gene level. For this task each 
	alignment file (BAM) is being read and all the necessary information is being extracted. In the end 
	all info will be summasised in an expression matrix. """
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} GENERATING THE EXPRESSION MATRIX (Gene Level)')


	os.chdir(expression_analysis_dir)  # Changing working directory 
	files = [os.path.basename(file).split(".")[0] for file in sorted(os.listdir()) if file.endswith(".gene-expression.tsv")]
	files.insert(0, "GeneID")
	files = "\t".join(files)
	subprocess.run("paste *.gene-expression.tsv | awk 'BEGIN {OFS=\"\t\"; FS=\"\t\"}; NR>4{{j=$1; for (i=2; i<=NF; i+=4){j=j FS $i} print j}}' > gene_expression_matrix.tab", shell=True)
	subprocess.run(f"sed -i '1i {files}' gene_expression_matrix.tab", shell=True)
	subprocess.run("paste *.gene-expression.tsv | awk 'BEGIN {OFS=\"\t\"; FS=\"\t\"}; NR<5{{j=$1; for (i=2; i<=NF; i+=4){j=j FS $i} print j}}' > sample_stats.tab", shell=True)
	subprocess.run(f"sed -i '1i {files}' sample_stats.tab", shell=True)
	return

def transcript_expression_matrix():
	""" Summarisation and construction of the expression matrix on gene level. For this task each 
	alignment file (BAM) is being read and all the necessary information is being extracted. In the end 
	all info will be summasised in an expression matrix. """
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} GENERATING THE EXPRESSION MATRIX (Isoform Level)')


	isoform_quant = sorted(glob.glob(f"{isoform_quant_dir}/*" + os.path.sep))
	sample_names = [paths.split("/")[-2] for paths in isoform_quant]
	run_salmon_merge = " ".join([
	salmon, "quantmerge",  # Calling salmon quantmerge
	"--quants", " ".join(isoform_quant),  # List of quantification directories
	"--names", " ".join(sample_names),  # List of names to give to the samples
	"--column tpm",  # The name of the column that will be merged; options are {len, elen, tpm, numreads}
	"--output", f"{isoform_quant_dir}/isoform_expression_matrix.tab",  #  Output quantification directory
	"2>>", os.path.join(pipeline_reports, "isoformquant_salmon_merge-report.txt")])  # Input files
	subprocess.run(run_salmon_merge, shell=True)
	return

def gene_fusion_matrix():
	""" Summarisation and construction of the gene fusion expression matrix 
	based on the fusion calls that were made with the Arriba software """ 
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} GENERATING THE EXPRESSION MATRIX (Gene Fusions)')

	
	fusion_data = glob.glob(f"{genefusion_dir}/*.fusion.tsv")

	# Obtaining all the gene fusions
	fusion_list = set()  # Keeping only unique fusions, not duplicates
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Obtaining the list of fusions detected from all samples...')
	for i, file in enumerate(fusion_data, 1):
		with open(file) as finf:
			for line in finf:
				if not line.startswith("#"):
					# Sorting the detected fusions alphabetically to avoid duplications
					fusions = sorted([line.strip().split("\t")[0], line.strip().split("\t")[1]])
					fusions = f"{fusions[0]}_{fusions[1]}"
					conf = line.strip().split("\t")[14]
					if args.low:  # Allow low confident fusions
						fusion_list.add(fusions)
					else:  # Do now allow low confident fusions
						if conf != "low":
							fusion_list.add(fusions)


	# Using the list of the obtained fusions to make a 
	# predefined dictionary (per sample) and fill if 
	# with values while reading each file.
	overall_dict = {}
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Analysing all {len(fusion_data)} samples...')
	for i, file in enumerate(fusion_data, 1):
		sample = os.path.basename(file)[:-11]
		with open(file) as fin:
			fusion_dict = {key:0 for key in sorted(fusion_list)}
			for line in fin:
				if not line.startswith("#"):
					fusions = sorted([line.strip().split("\t")[0], line.strip().split("\t")[1]])
					fusions = f"{fusions[0]}_{fusions[1]}"
					conf = line.strip().split("\t")[14]
					if args.low:
						if i==1:print("ATTENTION: Low confident fusions are allowed in the fusion matrix construction!")
						if fusions in fusion_dict:
							fusion_dict[fusions] += 1
					else:
						if conf != "low":
							fusion_dict[fusions] += 1
			overall_dict[sample] = fusion_dict.values()

	overall_dict = {k: overall_dict[k] for k in sorted(overall_dict)}  # Sort dictionary

	# Exporting the detected fusions 
	with open(f"{genefusion_dir}/gene_fusion_expression_matrix.tab", "w") as mat_out:
		mat_out.write("fusion\t{0}\n".format("\t".join(list(overall_dict.keys()))))
		writer = csv.writer(mat_out, delimiter='\t')
		writer.writerows(zip(sorted(fusion_list), *overall_dict.values()))
	return

def alternative_expression_matrix():
	""" To obtain the most abundant transcripts, 20 random healthy individuals will be picked up from 
	the clinical data file. Those 20 chosen samples will be used as a guide to elect the list of the 
	most abundant transcripts of each gene. Once we have this list, it will be used to transform the 
	transcript expression matrix into an alternative isoform expression matrix, where only the selected
	most abundant transcripts will be divided by the total tpms of all transcripts originating from the
	same gene. This way we will obtain the frequency of the alternative expression. """
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} GENERATING THE EXPRESSION MATRIX (Alternative Isoform Expression)')


	if not os.path.exists(altexpression_dir): os.makedirs(altexpression_dir)
	json_abundtranscripts = f"{altexpression_dir}/abundant_transcripts_db.json"

	
	# Obtaining the gene origin of each transcript
	database = {}
	with open(refAnnot) as refin:
		for line in refin:
			if not line.startswith("#"):
				if line.strip().split("\t")[2] == "transcript":
					database[line.strip().split("transcript_id \"")[1].split(";")[0].strip("\"")] = line.strip().split("gene_id \"")[1].split(";")[0].strip("\"")


	if not os.path.exists(json_abundtranscripts):
		most_abundant_transcripts = calc_abundant_transcripts(database)
		json.dump(most_abundant_transcripts, open(json_abundtranscripts,'w'))
	else:
		most_abundant_transcripts = json.load(open(json_abundtranscripts))
	
	# for a, b in most_abundant_transcripts.items():
	# 	# if a == "ENSG00000154518.9":
	# 	# if a == "ENSG00000000003.15":
	# 	if b in ["ENSG00000000003.15", "ENSG00000000005.6", "ENSG00000000419.13", "ENSG00000000457.14", 
	# 			 "ENSG00000000460.17", "ENSG00000000938.13", "ENSG00000000971.16", "ENSG00000001036.14", 
	# 			 "ENSG00000001084.13", "ENSG00000001167.14", "ENSG00000001460.18", "ENSG00000001461.17", 
	# 			 "ENSG00000001497.17", "ENSG00000001561.7"]:
	# 		print(b, a)
	

	# Import the isoform expression matrix with pandas as data frame
	expr_mat = pd.read_csv(f'{isoform_quant_dir}/isoform_expression_matrix.tab', delimiter='\t')
	# Adding the gene ID column in the data frame
	expr_mat['GeneID'] = expr_mat.apply(lambda row: database[row['Name']], axis=1)
	cols = list(expr_mat.columns.values)  # Obtaining the column names
	# Moving GeneID as first column in the data frame
	cols = [cols[-1]] + cols[0:-1]
	expr_mat = expr_mat[cols]  # Applying changes
	expr_mat = expr_mat.rename(columns={'Name': 'TranscriptID'})  # Renaming Name to transcript column
	expr_mat = expr_mat.sort_values(['GeneID','TranscriptID'], ascending=True).reset_index(drop=True)  # Sorting the data frame by GeneID and IsoformID
	# Set the index to ['GeneID', 'TranscriptID']
	expr_mat.set_index(['GeneID', 'TranscriptID'], inplace=True)

	# pd.set_option("display.max_rows", 100)
	# print(expr_mat.head(84))

	# v = expr_mat.head(84)
	# print("\n\n\n")

	# Obtaining the fraction of the most abundant transcript per gene 
	# divided by the sum of the transcripts belonging to that gene
	results = pd.DataFrame([])  # Saving the results to a new data frame
	for group_name, df_group in expr_mat.groupby('GeneID'):
		for index, row in df_group.iterrows():
			if index[1] in most_abundant_transcripts:  # Transcripts found in the most abundant list
				results = results.append(df_group.groupby('GeneID').apply(lambda x: row/x.sum()))

	# print(results)
	
	# Saving the alternative splicing matrix to a tab-delimited file
	results.round(6).to_csv(f'{altexpression_dir}/alternative_isoform_expression_matrix.tab', na_rep='NaN', sep='\t', index=True)  
	return

class mutation_profiling:

	def __init__(self):
		
		print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} GENE MUTATION PROFILING')
		

		if not os.path.exists(varcall_dir):os.makedirs(varcall_dir)
		if not os.path.exists(rnaediting_dir):os.makedirs(rnaediting_dir)
		if not os.path.exists(varcall_reports_dir):os.makedirs(varcall_reports_dir)

		aligned_data = glob.glob(f"{alignments_dir}/*.aligned.genome.bam")
		startTime_mut = datetime.now()
		for current_num, sample in enumerate(aligned_data, 1):
			sample_name = os.path.basename(sample).split(".")[0]
			print(f'\n{current_num}/{len(aligned_data)} | Processing sample {sample_name}...')
			# if current_num == 1:

			# Calling the different methods
			self.remove_duplicates(sample, sample_name)
			self.base_recalibration(sample, sample_name)
			self.variant_calling(sample, sample_name)
			self.vc_stats(sample, sample_name)
			# self.filterout_snps()
		
		# Running RediTools in a parallel fashion
		print(f'{datetime.now().strftime("%d.%m %H:%M")} 5/5 | REDItools (v2.0) - Analysing all {len(aligned_data)} samples for RNA editing events...')
		startTime_rnaedit = datetime.now()
		rnaediting_args = [(sample, os.path.basename(sample).split(".")[0])for sample in aligned_data]
		pool = mp.Pool(processes=int(args.threads))
		pool.map(self.RNA_editing, rnaediting_args)
		pool.close()
		print(f'\tRNA editing profiling finished after {datetime.now() - startTime_rnaedit}')

		print(f'\tMutation profiling finished after {datetime.now() - startTime_mut}')
		return

	def remove_duplicates(self, sample, sample_name):
		""" Removing duplicates from the aligned to the genome files """
		
		if not os.path.exists(temp):os.makedirs(temp)


		print(f'{datetime.now().strftime("%d.%m %H:%M")} 1/5 | Samtools markdup - Removing duplicate reads...')
		startTime = datetime.now()  # Recording the start time
		# Marking duplicate reads
		rmDuplicates = " ".join([
		"samtools markdup",  # Call samtools markdup
		"--threads", args.threads,  # Number of threads to be used
		"--output-fmt BAM",  # Output in BAM format
		"-r",  # Remove duplicate reads
		"-s",  # Report stats
		sample,
		f"{temp}/{sample_name}.rmdup.bam",
		"2>>", os.path.join(pipeline_reports, "removeDuplicates_report.txt")])
		subprocess.run(rmDuplicates, shell=True)

		# Indexing
		bam_index = " ".join([
		"samtools index",  # Call samtools markdup
		"-@", args.threads,  # Number of threads to be used
		f"{temp}/{sample_name}.rmdup.bam",
		"2>>", os.path.join(pipeline_reports, "postprocessing_bamindex-report.txt")])
		subprocess.run(bam_index, shell=True)
		print(f'\tRemoving duplicates finished in {datetime.now() - startTime}')
		return

	def base_recalibration(self, sample, sample_name):
		print(f'{datetime.now().strftime("%d.%m %H:%M")} 2/5 | GATK BaseRecalibrator - Base Recalibration...')
		known_sites = [f for f in glob.glob("{0}/*hg38.vcf".format(gatk_known_sites))]
		startTime = datetime.now()  # Recording the start time

		runBaseRecalibrator = " ".join([
		"gatk BaseRecalibrator",  # Call gatk BaseRecalibrator
		"--input", f"{temp}/{sample_name}.rmdup.bam",
		# "--intervals", targeted_regions_bed,  # Genomic intervals over which to operate
		"--known-sites", known_sites[0],  # databases of known polymorphic sites
		"--known-sites", known_sites[1],  # databases of known polymorphic sites
		"--known-sites", known_sites[2],  # databases of known polymorphic sites
		"--reference", refGenomeGRCh38,  # Reference sequence file
		"--output", f"{temp}/{sample_name}_recal_data.table",
		"2>>", os.path.join(pipeline_reports, "postprocessing_prebaserecalibrator-report.txt")])
		subprocess.run(runBaseRecalibrator, shell=True)

		# Generation of the recalibrated reads
		runApplyBQSR = " ".join([
		"gatk ApplyBQSR",  # Call gatk ApplyBQSR
		"--input", f"{temp}/{sample_name}.rmdup.bam",  # Input file containing sequence data
		# "--intervals", targeted_regions_bed,  # Genomic intervals over which to operate
		"--bqsr-recal-file", f"{temp}/{sample_name}_recal_data.table",  # File with base quality score recalibration
		"--output", f"{temp}/{sample_name}.rmdup.recal.bam",  # Output file
		"2>>", os.path.join(pipeline_reports, "postprocessing_applybqsr-report.txt")])
		subprocess.run(runApplyBQSR, shell=True)

		# Indexing
		bam_index = " ".join([
		"samtools index",  # Call samtools markdup
		"-@", args.threads,  # Number of threads to be used
		f"{temp}/{sample_name}.rmdup.recal.bam",
		"2>>", os.path.join(pipeline_reports, "postprocessing_bamindex-report.txt")])
		subprocess.run(bam_index, shell=True)
		print(f'\tBase recalibration finished in {datetime.now() - startTime}')

		### Removing unnecessary files
		subprocess.run(f"rm {temp}/*.rmdup.bam", shell=True)
		subprocess.run(f"rm {temp}/*.rmdup.bam.bai", shell=True)
		subprocess.run(f"rm {temp}/*data.table", shell=True)
		return

	def variant_calling(self, sample, sample_name):
		""" Calling Bcftools to perform InDel realignment on the fly and 
		  perform the  variant calling. Bcftools is a variant calling
		program for SNV, MNV, indels (<50 bp), and complex variants  """

		print(f'{datetime.now().strftime("%d.%m %H:%M")} 3/5 | Bcftools - Calling somatic variants ({sample_name})...')
		startTime = datetime.now()

		## 2. Run BcfTools to call variants
		bcftools = " ".join([
		"bcftools mpileup",
		"--min-MQ 20",  # Skip alignments with mapQ smaller than 20
		"--min-BQ 25",  # Skip bases with baseQ/BAQ smaller than 25
		"--redo-BAQ",  # Recalculate BAQ on the fly
		"--per-sample-mF",  # Apply -m and -F per-sample for increased sensitivity
		"--min-ireads 2",  # Minimum number gapped reads for indel candidates
		# "--regions-file", targeted_regions_bed,  # Analysis only on specific target genes
		"--annotate FORMAT/AD,FORMAT/DP,INFO/AD",  # Optional tags to output
		"-f", refGenomeGRCh38,  # faidx indexed reference sequence file
		"-Ou",  # Output 'u' uncompressed BCF
		os.path.join(temp, f"{sample_name}.rmdup.recal.bam"),  # Input data
		"| bcftools call",  # Calling Bcftools call
		"--threads", args.threads,  # Number of threads to use
		"--ploidy GRCh38",  # Ploidy (2)
		"--output-type z",  # Output as 'z'; compressed VCF 
		"--variants-only",  # Output variant sites only
		"--multiallelic-caller",  # alternative model for multiallelic and rare-variant calling
		"--output", f"{varcall_dir}/{sample_name}.bcftools.vcf.gz",
		"2>>", os.path.join(pipeline_reports, "variantcalling_bcftools_report.txt")])
		subprocess.run(bcftools, shell=True)
		subprocess.run(f"tabix -p vcf {varcall_dir}/{sample_name}.bcftools.vcf.gz", shell=True)
		print(f'\tBcfTools finished after {datetime.now() - startTime}')
		return
		
	def vc_stats(self, sample, sample_name):
		""" Generating basic statistics and plots """
		print(f'{datetime.now().strftime("%d.%m %H:%M")} 4/5 | Bcftools stats - Generating basic stats ({sample_name})...')

		bcfstats = " ".join([
		"bcftools stats",  # Calling bcfstats
		"--threads", args.threads,  # Number of threads to use for decompressing
		"--fasta-ref", refGenomeGRCh38,  # Indexed reference sequence file to determine INDEL context
		# "--regions-file", targeted_regions_bed,  # Analysis only on specific target genes
		f"{varcall_dir}/{sample_name}.bcftools.vcf.gz",  # VCF file
		">", f"{varcall_reports_dir}/{sample_name}.vcf.stats.txt",
		"2>>", os.path.join(pipeline_reports, "variantcalling_bcfstats-report.txt")])
		subprocess.run(bcfstats, shell=True)
		
		# vcfplot = " ".join([
		# "plot-vcfstats",  # Calling plot-vcfstats
		# "--prefix", f"{varcall_reports_dir}/{sample_name}.",  # Output directory/samplename
		# f"{varcall_reports_dir}/{sample_name}.vcf.stats.txt",
		# "2>>", os.path.join(pipeline_reports, "variantcalling_vcfplot-report.txt")])
		# subprocess.run(vcfplot, shell=True)
		return
	
	def RNA_editing(self, func_args):
		""" Generating basic statistics and plots """
		sample = func_args[0]  # Obtaining the sample
		sample_name = func_args[1]  # Obtaining the sample name
		print(f'{datetime.now().strftime("%d.%m %H:%M")} 5/5 | REDItools (v2.0) - Analysing sample {sample_name} for RNA editing events...')
		

		runReditools = " ".join([
		reditools,  # Calling REDItools
		"--file", f"{temp}/{sample_name}.rmdup.recal.bam",  # The bam file to be analyzed
		"--strict",  # Activate strict mode: only sites with edits will be included in the output
		"--strand 0",  # Strand specific protocol: 0 an unstranded protocol was followed
		"--reference", refGenomeGRCh38,  # The reference FASTA file
		"--min-edits 10",  # Minimum edit per position
		"--min-read-quality 20",  # Reads whose mapping quality is below this value will be discarded
		"--min-base-quality  25",  # Bases whose quality is below 25 will not be included in the analysis
		"--output-file", f"{rnaediting_dir}/{sample_name}.out.table.txt",  # The output statistics file
		"2>>", os.path.join(pipeline_reports, "rnaediting_reditools-report.txt")])
		subprocess.run(runReditools, shell=True)
		return 

def snp_matrix():
	""" Summarising all SNPs found during the mutation profiling in 
	a single overall matrix (expression matrix) """ 
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} GENERATING THE EXPRESSION MATRIX (SNP analysis)')

	
	snp_data = ' '.join(sorted(glob.glob(f"{varcall_dir}/*.bcftools.vcf.gz")))
	# Merging all vcfs to one file 
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  BcfTools merge (v1.7) - Merging all VCF files into one and indexing...')
	subprocess.run(f'bcftools merge --output-type z --threads {args.threads} --output {varcall_dir}/merged_profiles.vcf.gz {snp_data}',shell = True)
	
	# Indexing the merged file
	subprocess.run(f'tabix -p vcf {varcall_dir}/merged_profiles.vcf.gz',shell = True)

	# Convert the merged vcf to tab-delimited matrix
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  gatk VariantsToTable (v4.1.9.0) - Converting the merged VCF file to a tab-delimited file...')
	# Obtaining the AD: Allelic depths for the ref and alt alleles
	# Explanation: Allele specific depth. i.e. if ref is 'A' and alt is 'G' and AD is '6,9' you got 6 A reads and 9 G reads
	# Genotype as output instead of Allelic depth
	subprocess.run(f'gatk VariantsToTable -V {varcall_dir}/merged_profiles.vcf.gz -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT -O {varcall_dir}/variant_table.tab',shell = True)
	
	# Removing unnecessary files
	subprocess.run(f'rm {varcall_dir}/merged_profiles.vcf.gz {varcall_dir}/merged_profiles.vcf.gz.tbi', shell=True)
	return

def RNA_editing_matrix():
	""" Summarisation and construction of the RNA editing expression matrix 
	based on the calls that were made with the RediTools software """ 
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} GENERATING THE EXPRESSION MATRIX (RNA editing events)')

	
	rnaediting_data = sorted(glob.glob(f"{rnaediting_dir}/*.out.table.txt"))

	rna_edit_expression = []
	for file in rnaediting_data:
		samplename = os.path.basename(file)[:-14]
		# Reading each RNA editing table file and keeping events with frequency higher than 0.80 and mean quality higher than 20
		df = (pd.read_csv(file, sep='\t', usecols=['Region','Position','AllSubs','MeanQ','Frequency']) [lambda x:  (x['Frequency'] >=0.80) & (x['MeanQ'] >= 25)])
		df = df.drop(['MeanQ', 'Frequency'], axis = 1)  # Removing MeanQ and Frequency column
		df = df.rename(columns={'Region': 'Chromosome', 'AllSubs': 'Substitution'})  # Renaming 2 columns
		df = df.set_index(['Chromosome','Position','Substitution'])  # Setting the index columns
		df[samplename] = 1  # Setting 1 to all the RNA editing occurrences found in the file
		rna_edit_expression.append(df)  # Appending the recorded events to a list
	
	# Concatenating all per file data to a single matrix
	rna_edit_expression = pd.concat(rna_edit_expression, axis=1).fillna(0)

	rna_edit_expression = rna_edit_expression.sort_values(by="Chromosome", key=lambda x: np.argsort(index_natsorted(x)))
	# rna_edit_expression = rna_edit_expression.sort_values(['Chromosome', 'Position'], ascending=[True, True])  # Sorting the data frame

	# Saving the RNA-editing matrix to a file
	rna_edit_expression.to_csv(f'{rnaediting_dir}/RNAediting_expression_matrix.tab', sep='\t', index=True)
	rna_edit_expression = rna_edit_expression.reset_index()  # Resetting indexes to be able to export those columns in a tsv
	rna_edit_expression.to_csv(f'{varcall_dir}/RNAediting_events_to_exclude.tsv', sep='\t', columns = ['Chromosome','Position','Substitution'], index=False)
	return

def final_snp_matrix():
	""" Removing all RNA editing events from the SNP expression table
	and generating the final SNP expression matrix  """ 
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} GENERATING THE EXPRESSION MATRIX\nGenerating the final SNP expression matrix by filtering out all RNA editing events')
	

	snps_to_remove = []
	# Reading the RNAediting_events_to_exclude and save them to a list
	with open(f'{varcall_dir}/RNAediting_events_to_exclude.tsv', 'r') as snpin:
		for line in snpin:
			if not line.startswith("Chromosome"):
				chrom = line.strip().split("\t")[0]
				pos = line.strip().split("\t")[1]
				subs = line.strip().split("\t")[2]
				elm = f'{chrom}_{pos}_{subs}'
				snps_to_remove.append(elm)


	# Reading the SNP expression matrix and rewriting it by excluding 
	# the rows that were found in the RNA editing matrix
	with open(f'{varcall_dir}/variant_table.tab', 'r') as matin, open(f'{varcall_dir}/snp_expression_matrix.tab', 'w') as matout:
		for line in matin:
			if line.startswith("CHROM"):
				header = line.replace(".GT","").replace("CHROM","Chromosome").replace("POS","Position")
				header = header.replace("REF","Ref").replace("ALT","Alt").replace("TYPE","Type")
				matout.write(header)
			else:
				chrmosome = line.strip().split("\t")[0]
				position = line.strip().split("\t")[1]
				ref = line.strip().split("\t")[2]
				alt = line.strip().split("\t")[3]
				element = f'{chrmosome}_{position}_{ref}{alt}'
				if not element in snps_to_remove:
					matout.write(line.replace("./.","0"))
				# else:
				# 	print(element)
	return


def calc_abundant_transcripts(database):
	""" Helper function
		Calculating the most abundant transcript per gene based on 
		20 random Healthy Control samples. We will use these transcripts
		to calculate the alternative isoform expression frequency. """
	
	random_healthyIndividuals = []  # List that will host the healthy individuals
	### Reading the project's clinical data to obtain the healthy individuals
	with open(clinical_data) as fin:
		for line in fin:
			line = line.strip()
			if not line.startswith("Experiment"):
				if line.split("\t")[-1] == "HC":
					sample_id = "{0}_{1}".format(line.split("\t")[0], line.split("\t")[-1])
					random_healthyIndividuals.append(f"{isoform_quant_dir}/{sample_id}/quant.sf")

	###  Picking up 20 individuals randomly and saving their path to random_healthyIndividuals list
	if len(random_healthyIndividuals) >= 20:
		random_healthyIndividuals = random.sample(random_healthyIndividuals, 20)
	else:
		random_healthyIndividuals = random_healthyIndividuals

	# random_healthyIndividuals = ["/home/stavros/playground/platelets_analysis/analysis/isoform_expression/SRX2370887_HC/quant.sf",
	# 							 "/home/stavros/playground/platelets_analysis/analysis/isoform_expression/SRX2370888_HC/quant.sf",
	# 							 "/home/stavros/playground/platelets_analysis/analysis/isoform_expression/SRX2370889_HC/quant.sf"]

	# Calculating the collective tpm per transcript per sample
	abundant_transcripts = {}
	for qunat_file in random_healthyIndividuals:
		with open(qunat_file, "r") as qin:
			for line in qin:
				if not line.startswith("Name"):
					transcript_id = str(line.strip().split("\t")[0])
					gene_id = str(database[transcript_id])
					tpm = float(line.strip().split("\t")[3])

					if gene_id in abundant_transcripts:
						if [update_value(abundant_transcripts, gene_id, idx, tpm) for idx, (isoform, value) in enumerate(abundant_transcripts[gene_id]) if transcript_id == isoform]:
							pass
						else:	
							abundant_transcripts[gene_id].append([transcript_id,tpm])
					else:
						abundant_transcripts[gene_id] = [[transcript_id,tpm]]
	
	# for a, b in abundant_transcripts.items():
	# 	# if a == "ENSG00000154518.9":
	# 	if a == "ENSG00000000003.15":
	# 		print(a,b)

	# Obtaining the most abundant transcript from each gene 
	# and generating a list with the winner transcript.
	most_abundant_transcripts = {}
	for gene, transc_n_count in abundant_transcripts.items():
			if len(transc_n_count) == 1:
				most_abundant_transcripts[transc_n_count[0][0]] = gene
			else:
				most_abundant = max(transc_n_count, key = itemgetter(1))[0]
				most_abundant_transcripts[most_abundant] = gene
	return most_abundant_transcripts

def update_value(abundant_dict, gene_id, index, tpm):
	abundant_dict[gene_id][index][1] += tpm
	return abundant_dict

def summary():
	## REMOVING UNNECESSARY FILES & REPORTS
	for path, subdir, folder in os.walk(pipeline_reports):
		for name in folder:
			file = os.path.join(path, name)
			if os.stat(file).st_size == 0:
				os.remove(file)
	shutil.rmtree(temp)
	return 



def main():

	# ### Performing preprocessing of the data if it isn't already being done
	# preprocessing_samples()  # Preprocessing the raw reads

	# ### Performing multilevel analysis including gene and isoform level quantification
	# geneNisoform_level_analysis()

	# mapping_quality_control()

	# gene_expression_matrix()

	# transcript_expression_matrix()

	# gene_fusion_matrix()

	# alternative_expression_matrix()

	# ### Variant profiling analysis
	# mutation_profiling()

	# snp_matrix()

	RNA_editing_matrix()

	# final_snp_matrix()
	# summary()

	print(f'The pipeline finished after {datetime.now() - startTime}')
	
if __name__ == "__main__": main()
