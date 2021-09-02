"""
Juno_assembly pipeline
Authors: Alejandra Hernandez-Segura, Robert Verhagen, Ernst Hamer, Dennis Schmitz, Diogo Borst, Tom van Wijk, Maaike van der Beld
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium Surveillance (IDS), Bacteriologie (BPD)
Date: 26-08-2021

Documentation: https://rivm-bioinformatics.github.io/ids_bacteriology_man/juno-assembly.html

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

# SAMPLES is a dict with sample in the form:
# sample > read number > file. E.g.: SAMPLES["sample_1"]["R1"] = "x_R1.gz"
sample_sheet = config["sample_sheet"]
SAMPLES = {}
with open(sample_sheet) as sample_sheet_file:
    SAMPLES = yaml.safe_load(sample_sheet_file) 

OUT = config["out"]
IN = config["input_dir"]

#@################################################################################
#@#### 				Processes                                    #####
#@################################################################################

    #############################################################################
    ##### Data quality control and cleaning                                 #####
    #############################################################################
include: "bin/rules/fastqc_raw_data.smk"
include: "bin/rules/clean_fastq.smk"
include: "bin/rules/fastqc_clean_data.smk"

    #############################################################################
    ##### De novo assembly                                                  #####
    #############################################################################
include: "bin/rules/de_novo_assembly.smk"

    #############################################################################
    ##### Species identification                                            #####
    #############################################################################
include: "bin/rules/identify_species.smk"

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
                    'files/trimmomatic_0.36_adapters_lists/NexteraPE-PE.fa',
                    'files/accepted_genera_checkm.txt',
                    'files/multiqc_config.yaml' ]
    for filename in important_files:
        if not os.path.exists(filename):
            raise FileNotFoundError(filename)

#@################################################################################
#@#### These are the conditional cleanup rules                               #####
#@################################################################################

# onsuccess:
#     shell("""
#     """)


#################################################################################
##### Specify final output:                                                 #####
#################################################################################

localrules:
    all


rule all:
    input:
        expand(OUT + "/qc_raw_fastq/{sample}_{read}_fastqc.zip", sample = SAMPLES, read = ['R1', 'R2']),   
        expand(OUT + "/clean_fastq/{sample}_{read}.fastq.gz", sample = SAMPLES, read = ['pR1', 'pR2']),
        expand(OUT + "/qc_clean_fastq/{sample}_{read}_fastqc.zip", sample = SAMPLES, read = ['pR1', 'pR2']),
        expand(OUT + "/de_novo_assembly/{sample}/scaffolds.fasta", sample = SAMPLES),   
        OUT + "/qc_de_novo_assembly/quast/report.tsv",
        expand(OUT + "/qc_de_novo_assembly/bbtools_scaffolds/per_sample/{sample}_MinLenFiltSummary.tsv", sample = SAMPLES),
        OUT + "/qc_de_novo_assembly/bbtools_scaffolds/bbtools_scaffolds.tsv",
        OUT + "/qc_de_novo_assembly/bbtools_scaffolds/bbtools_summary_report.tsv",
        expand(OUT + '/identify_species/{sample}_species_content.txt', sample = SAMPLES),
        expand(OUT + '/identify_species/{sample}_bracken_species.kreport2', sample = SAMPLES),
        OUT + "/multiqc/multiqc.html"
