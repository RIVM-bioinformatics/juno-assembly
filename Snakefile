"""
Bac_Gastro pipeline
Authors: Ernst Hamer, Dennis Schmitz, Robert Verhagen, Diogo Borst, Tom van Wijk
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: IDS - BPD - Bacteriology
Date: 24-02-2020

Changelog, examples, installation guide and explanation on:
    https://github.com/ELAHamer/BAC_gastro


Snakemake rules (in order of execution):
    1 fastQC        # Asses quality of raw reads.
    2 trimmomatic   # Trim low quality reads and adapter sequences.
    3 fastQC        # Asses quality of trimmed reads.
    4 spades        # Perform assembly with SPAdes.
    5 quast         # Run quality control tool QUAST on contigs/scaffolds.
    6 checkM        # Gives scores for completeness, contamination and strain heterogeneity.
    7 picard        # Determines library fragment lengths.
    8 bbmap         # Generate scaffold alignment metrics.
    9 multiQC       # Summarize analysis results and quality assessments in a single report 

Custom configuration options (passed via config.yaml or via 
    `--config` command line option):
        * sourcedata (Directory where input files can be found)
        * runsheet (YAML file with sample info, see format below)
        * out (Directory where output is written to)
        * excel file with genus for each sample
"""

#################################################################################
##### Import config file, sample_sheet and set output folder names          #####
#################################################################################

configfile: "profile/pipeline_parameters.yaml"
configfile: "profile/variables.yaml"

from pandas import *
import pathlib
import pprint
import os
import yaml
import json


yaml.warnings({'YAMLLoadWarning': False}) # Suppress yaml "unsafe" warnings

#################################################################################
##### Load samplesheet, load genus dict and define output directory         #####
#################################################################################

# SAMPLES is a dict with sample in the form sample > read number > file. E.g.: SAMPLES["sample_1"]["R1"] = "x_R1.gz"
SAMPLES = {}
with open(config["sample_sheet"]) as sample_sheet_file:
    SAMPLES = yaml.load(sample_sheet_file) 

# OUT defines output directory for most rules.
OUT = pathlib.Path(config["out"])

# make a list of all genera supported by the current version of CheckM
genuslist = []
with open(config["genuslist"]) as file_in:
    for line in file_in:
        if "genus" in line:
            genuslist.append((line.split()[1].lower().strip()))

#GENUS added to samplesheet dict (for CheckM)
xls = ExcelFile(pathlib.Path(config["genus_file"]))
df1 = xls.parse(xls.sheet_names[0])[['Monsternummer','genus']]
genus_dict = dict(zip(df1['Monsternummer'].values.tolist(), df1['genus'].values.tolist()))
genus_dict = json.loads(json.dumps(genus_dict), parse_int=str) # Convert all dict values and keys to strings


#################################################################################
##### Catch sample and genus errors, when not specified by the user         #####
#################################################################################

error_samples_genus = []
error_samples_sample = []

#search samples in genus dict and add genus to the SAMPLES dict
for sample, value in SAMPLES.items():
    if str(sample) in genus_dict:
        if str(genus_dict[sample]).lower() in genuslist:
            SAMPLES[sample] = [value,genus_dict[sample]]

        # Genus not recognized by checkM
        else: 
            error_samples_genus.append(sample)

    # Sample not found in Excel file
    else: 
        error_samples_sample.append(sample)



if error_samples_sample:
    print(f""" \n\nERROR: The sample(s):\n\n{chr(10).join(error_samples_sample)} \n
    Not found in the Excel file: {pathlib.Path(config["genus_file"])}. 
    Please insert the samples with it’s corresponding genus in the Excel file before starting the pipeline.
    When the samples are in the Excel file, checkM can asses the quality of the microbial genomes. \n
    It is also possible to remove the sample that causes the error from the samplesheet, and run the analysis without this sample.\n\n""")
    sys.exit(1)
        
if error_samples_genus:
    print(f""" \n\nERROR:  The genus supplied with the sample(s):\n\n{chr(10).join(error_samples_genus)}\n\nWhere not recognized by CheckM\n
    Please supply the sample row in the Excel file {pathlib.Path(config["genus_file"])}
    with a correct genus. If you are unsure what genera are accepted by the current
    version of the pipeline, please consult the checkm_taxon_list.txt for available genera.\n\n""")
    sys.exit(1)
        

#@################################################################################
#@#### 				Processes                                    #####
#@################################################################################

    #############################################################################
    ##### Data quality control and cleaning                                 #####
    #############################################################################
include: "bin/rules/qc_raw_data.smk"
include: "bin/rules/clean_data.smk"
include: "bin/rules/qc_clean_data.smk"
include: "bin/rules/cat_unpaired.smk"
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
include: "bin/rules/fragment_length_analysis.smk"
include: "bin/rules/generate_contig_metrics.smk"
include: "bin/rules/parse_bbtools.smk"
include: "bin/rules/parse_bbtools_summary.smk"
include: "bin/rules/multiqc_report.smk"



#@################################################################################
#@#### The `onstart` checker codeblock                                       #####
#@################################################################################

onstart:
    try:
        print("Checking if all specified files are accessible...")
        for filename in [ config["sample_sheet"],
                         config["genus_file"],
                         config["genuslist"],
                         'files/trimmomatic_0.36_adapters_lists/NexteraPE-PE.fa' ]:
            if not os.path.exists(filename):
                raise FileNotFoundError(filename)
    except FileNotFoundError as e:
        print("This file is not available or accessible: %s" % e)
        sys.exit(1)
    else:
        print("\tAll specified files are present!")
    shell("""
        mkdir -p {OUT}
        mkdir -p {OUT}/results
        echo -e "\nLogging pipeline settings..."
        echo -e "\tGenerating methodological hash (fingerprint)..."
        echo -e "This is the link to the code used for this analysis:\thttps://github.com/DennisSchmitz/BAC_gastro/tree/$(git log -n 1 --pretty=format:"%H")" > '{OUT}/results/log_git.txt'
        echo -e "This code with unique fingerprint $(git log -n1 --pretty=format:"%H") was committed by $(git log -n1 --pretty=format:"%an <%ae>") at $(git log -n1 --pretty=format:"%ad")" >> '{OUT}/results/log_git.txt'
        echo -e "\tGenerating full software list of current Conda environment (\"bac_gastro_master\")..."
        conda list > '{OUT}/results/log_conda.txt'
        echo -e "\tGenerating config file log..."
        rm -f '{OUT}/results/log_config.txt'
        for file in profile/*.yaml
        do
            echo -e "\n==> Contents of file \"${{file}}\": <==" >> '{OUT}/results/log_config.txt'
            cat ${{file}} >> '{OUT}/results/log_config.txt'
            echo -e "\n\n" >> '{OUT}/results/log_config.txt'
        done
    """)

#@################################################################################
#@#### These are the conditional cleanup rules                               #####
#@################################################################################

#onerror:
 #   shell("""""")


onsuccess:
    shell("""
        echo -e "\tGenerating HTML index of log files..."
        echo -e "\tGenerating Snakemake report..."
        snakemake --unlock
        snakemake --report '{OUT}/results/snakemake_report.html'
        echo -e "Finished"
    """)



#################################################################################
##### Specify final output:                                                 #####
#################################################################################

localrules:
    all,
    cat_unpaired


rule all:
    input:
        expand(str(OUT / "FastQC_pretrim/{sample}_{read}_fastqc.zip"), sample = SAMPLES, read = ['R1', 'R2']),   
        expand(str(OUT / "trimmomatic/{sample}_{read}.fastq"), sample = SAMPLES, read = ['pR1', 'pR2', 'uR1', 'uR2']),
        expand(str(OUT / "FastQC_posttrim/{sample}_{read}_fastqc.zip"), sample = SAMPLES, read = ['pR1', 'pR2', 'uR1', 'uR2']),
        expand(str(OUT / "SPAdes/{sample}/scaffolds.fasta"), sample = SAMPLES),   
        str(OUT / "QUAST/report.tsv"),
        expand(str(OUT / "CheckM/per_sample/{sample}/CheckM_{sample}.tsv"), sample=SAMPLES, genus = "genus"),
        str(OUT / "CheckM/CheckM_combined/CheckM_report.tsv"),   
        expand(str(OUT / "scaffolds_filtered/{sample}_sorted.bam"), sample = SAMPLES),
        expand(str(OUT / "bbtools_scaffolds/per_sample/{sample}_MinLenFiltSummary.tsv"), sample = SAMPLES),
        str(OUT / "bbtools_scaffolds/bbtools_combined/bbtools_scaffolds.tsv"),
        str(OUT / "bbtools_scaffolds/bbtools_combined/bbtools_summary_report.tsv"),
        str(OUT / "MultiQC/multiqc.html") 







    
    
    
    


        










