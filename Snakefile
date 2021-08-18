"""
Juno pipeline
Authors: Ernst Hamer, Alejandra Hernandez-Segura, Dennis Schmitz, Robert Verhagen, Diogo Borst, Tom van Wijk, Maaike van der Beld
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium Surveillance (IDS), Bacteriologie (BPD)
Date: 09-10-2020

Documentation: https://github.com/DennisSchmitz/BAC_gastro


Snakemake rules (in order of execution):
    1 fastQC        # Asses quality of raw reads.
    2 trimmomatic   # Trim low quality reads and adapter sequences.
    3 fastQC        # Asses quality of trimmed reads.
    4 spades        # Perform assembly with SPAdes.
    5 quast         # Run quality control tool QUAST on contigs/scaffolds.
    6 checkM        # Gives scores for completeness, contamination and strain heterogeneity (optional).
    7 picard        # Determines library fragment lengths.
    8 bbmap         # Generate scaffold alignment metrics.
    9 multiQC       # Summarize analysis results and quality assessments in a single report 

"""

#################################################################################
##### Import config file, sample_sheet and set output folder names          #####
#################################################################################

from pandas import *
import pathlib
import pprint
import os
import yaml
import json

#################################################################################
##### Load samplesheet, load genus dict and define output directory         #####
#################################################################################

# SAMPLES is a dict with sample in the form sample > read number > file. E.g.: SAMPLES["sample_1"]["R1"] = "x_R1.gz"
sample_sheet = config["sample_sheet"]
SAMPLES = {}
with open(sample_sheet) as sample_sheet_file:
    SAMPLES = yaml.safe_load(sample_sheet_file) 

# OUT defines output directory for most rules.
OUT = config["out"]

#@################################################################################
#@#### 				Processes                                    #####
#@################################################################################

    #############################################################################
    ##### Data quality control and cleaning                                 #####
    #############################################################################
include: "bin/rules/fastqc_raw_data.smk"
include: "bin/rules/trimmomatic.smk"
include: "bin/rules/fastqc_clean_data.smk"

    #############################################################################
    ##### De novo assembly                                                  #####
    #############################################################################
include: "bin/rules/run_spades.smk"

    #############################################################################
    ##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
    #############################################################################
include: "bin/rules/run_quast.smk"

include: "bin/rules/run_checkm.smk"
include: "bin/rules/parse_checkm.smk"

include: "bin/rules/picard_insert_size.smk"
include: "bin/rules/generate_contig_metrics.smk"
include: "bin/rules/parse_bbtools.smk"
include: "bin/rules/parse_bbtools_summary.smk"

include: "bin/rules/multiqc.smk"


#@################################################################################
#@#### The `onstart` checker codeblock                                       #####
#@################################################################################

onstart:
    print("Checking if all specified files are accessible...")
    important_files = [ config["sample_sheet"],
                    'files/trimmomatic_0.36_adapters_lists/NexteraPE-PE.fa' ]
    for filename in important_files:
        if not os.path.exists(filename):
            raise FileNotFoundError(filename)

#@################################################################################
#@#### These are the conditional cleanup rules                               #####
#@################################################################################

onsuccess:
    shell("""
        echo -e "\tGenerating Snakemake report..."
        snakemake --config sample_sheet={sample_sheet} \
                    --configfile config/pipeline_parameters.yaml config/user_parameters.yaml \
                    --cores 1 --unlock
        snakemake --config sample_sheet={sample_sheet} \
                    --configfile config/pipeline_parameters.yaml config/user_parameters.yaml \
                    --cores 1 --report '{OUT}/audit_trail/snakemake_report.html'
        echo -e "Finished"
    """)


#################################################################################
##### Specify final output:                                                 #####
#################################################################################

localrules:
    all


rule all:
    input:
        expand(OUT + "/qc_raw_fastq/{sample}_{read}_fastqc.zip", sample = SAMPLES, read = ['R1', 'R2']),   
        expand(OUT + "/clean_fastq/{sample}_{read}.fastq.gz", sample = SAMPLES, read = ['pR1', 'pR2', 'uR1', 'uR2']),
        expand(OUT + "/qc_clean_fastq/{sample}_{read}_fastqc.zip", sample = SAMPLES, read = ['pR1', 'pR2']),
        expand(OUT + "/de_novo_assembly/{sample}/scaffolds.fasta", sample = SAMPLES),   
        OUT + "/qc_de_novo_assembly/quast/report.tsv",
        expand(OUT + "/qc_de_novo_assembly/bbtools_scaffolds/per_sample/{sample}_MinLenFiltSummary.tsv", sample = SAMPLES),
        OUT + "/qc_de_novo_assembly/bbtools_scaffolds/bbtools_scaffolds.tsv",
        OUT + "/qc_de_novo_assembly/bbtools_scaffolds/bbtools_summary_report.tsv",
        OUT + "/multiqc/multiqc.html"
