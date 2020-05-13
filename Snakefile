"""
Bac_Gastro pipeline
Authors: Ernst Hamer, Robert Verhagen, Diogo Borst, Tom van Wijk
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
    7 piccard       # Determines library fragment lengths.
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
import yaml


#################################################################################
##### Load samplesheet, load genus dict and define output directory         #####
#################################################################################

# SAMPLES is a dict with sample in the form sample > read number > file. E.g.: SAMPLES["sample_1"]["R1"] = "x_R1.gz"
SAMPLES = {}
with open(config["sample_sheet"]) as sample_sheet_file:
    SAMPLES = yaml.load(sample_sheet_file) 

# OUT defines output directory for most rules.
OUT = pathlib.Path(config["out"])
genuslist = (config["genuslist"]).split(",")

#GENUS added to samplesheet dict (for CheckM)
xls = ExcelFile(pathlib.Path(config["genus_file"]))
df1 = xls.parse(xls.sheet_names[0])[['Monsternummer','genus']]
genus_dict = dict(zip(df1['Monsternummer'].values.tolist(), df1['genus'].values.tolist())) 

#search samples in genus dict and add genus to the SAMPLES dict
for sample, value in SAMPLES.items():
    if int(sample) in genus_dict:
        if genus_dict[int(sample)].lower() in genuslist:
            SAMPLES[sample] = [value,genus_dict[int(sample)]]
        # if genus is not listeria, shigella, salmonella or escherichia:
        else:
            SAMPLES[sample] = [value,"nan"]
            print(f"Sample {sample} is not supplied with a genus supported by this version of the pipeline")
    # if 
    else:
        SAMPLES[sample] = [value,"nan"]
        print(f"Sample {sample} is not supplied with a genus supported by this version of the pipeline")
        


#################################################################################
##### Specify final output:                                                 #####
#################################################################################

localrules:
    all,
    cat_unpaired,


rule all:
    input:
        expand(str(OUT / "FastQC_pretrim/{sample}_{read}_fastqc.zip"), sample = SAMPLES, read = ['R1', 'R2']),   
        expand(str(OUT / "trimmomatic/{sample}_{read}.fastq"), sample = SAMPLES, read = ['pR1', 'pR2', 'uR1', 'uR2']),
        expand(str(OUT / "FastQC_posttrim/{sample}_{read}_fastqc.zip"), sample = SAMPLES, read = ['pR1', 'pR2', 'uR1', 'uR2']),
        expand(str(OUT / "SPAdes/{sample}/scaffolds.fasta"), sample = SAMPLES),   
        expand(str(OUT / "QUAST/per_sample/{sample}/report.html"), sample=SAMPLES),
        str(OUT / "QUAST/combined/report.tsv"),
        expand(str(OUT / "CheckM/per_sample/{sample}/CheckM_{sample}.tsv"), sample=SAMPLES, genus = "genus"),
        str(OUT / "CheckM/CheckM_combined/CheckM_report.tsv"),   
        expand(str(OUT / "scaffolds_filtered/{sample}_sorted.bam"), sample = SAMPLES),
        expand(str(OUT / "bbtools_scaffolds/per_sample/{sample}_MinLenFiltSummary.tsv"), sample = SAMPLES),
        str(OUT / "bbtools_scaffolds/bbtools_combined/bbtools_scaffolds.tsv"),
        str(OUT / "bbtools_scaffolds/bbtools_combined/bbtools_summary_report.tsv"),
        str(OUT / "MultiQC/multiqc.html"),


#################################################################################
##### sub-processes                                                         #####
#################################################################################

    #############################################################################
    ##### Data quality control and cleaning                                 #####
    #############################################################################

rule QC_raw_data:
    input:
        lambda wildcards: SAMPLES[wildcards.sample][0][wildcards.read],
    output:
        html=str(OUT / "FastQC_pretrim/{sample}_{read}_fastqc.html"),
        zip=str(OUT / "FastQC_pretrim/{sample}_{read}_fastqc.zip")
    conda:
        "environments/QC_and_clean.yaml"
    benchmark:
        str(OUT / "log/benchmark/QC_raw_data_{sample}_{read}.txt")
    threads: 1
    log:
        str(OUT / "log/fastqc/QC_raw_data_{sample}_{read}.log")
    params:
        output_dir=str(OUT / "FastQC_pretrim/")
    shell:
        """
        bash bin/fastqc_wrapper.sh {input} {params.output_dir} {output.html} {output.zip} {log} 
        """

rule Clean_the_data:
    input:
        lambda wildcards: (SAMPLES[wildcards.sample][0][i] for i in ("R1", "R2"))
    output:
        r1=str(OUT / "trimmomatic/{sample}_pR1.fastq"),
        r2=str(OUT / "trimmomatic/{sample}_pR2.fastq"),
        r1_unpaired=str(OUT / "trimmomatic/{sample}_uR1.fastq"),
        r2_unpaired=str(OUT / "trimmomatic/{sample}_uR2.fastq"),
    conda:
        "environments/QC_and_clean.yaml"
    benchmark:
        str(OUT / "log/benchmark/Clean_the_data_{sample}.txt")
    threads: config["threads"]["Clean_the_data"]
    log:
        str(OUT / "log/trimmomatic/Clean_the_data_{sample}.log")
    params:
        adapter_removal_config=config["Trimmomatic"]["adapter_removal_config"],
        quality_trimming_config=config["Trimmomatic"]["quality_trimming_config"],
        minimum_length_config=config["Trimmomatic"]["minimum_length_config"],
    shell:
        """
trimmomatic PE -threads {threads} \
{input[0]:q} {input[1]:q} \
{output.r1} {output.r1_unpaired} \
{output.r2} {output.r2_unpaired} \
{params.adapter_removal_config} \
{params.quality_trimming_config} \
{params.minimum_length_config} > {log} 2>&1
touch -r {output.r1} {output.r1_unpaired}
touch -r {output.r2} {output.r2_unpaired}
        """
# touch the output in the case 100% of the reads are trimmed and no file is created

rule QC_clean_data:
    input:
        str(OUT / "trimmomatic/{sample}_{read}.fastq")
    output:
        html=str(OUT / "FastQC_posttrim/{sample}_{read}_fastqc.html"),
        zip=str(OUT / "FastQC_posttrim/{sample}_{read}_fastqc.zip")
    conda:
        "environments/QC_and_clean.yaml"
    benchmark:
        str(OUT / "log/benchmark/QC_clean_data_{sample}_{read}.txt")
    threads: 1
    log:
        str(OUT / "log/fastqc/QC_clean_data_{sample}_{read}.log")
    params:
        output_dir=str(OUT / "FastQC_posttrim/")
    shell:
        """
if [ -s "{input}" ] # If file exists and is NOT empty (i.e. filesize > 0) do...
then
    fastqc --quiet --outdir {params.output_dir} {input} > {log}
else
    touch {output.html}
    touch {output.zip}
fi
    """

# Create a concatenated unpaired fastq file to support SPAdes de novo assembly
rule cat_unpaired:
    input:
        r1_unpaired=str(OUT / "trimmomatic/{sample}_uR1.fastq"),
        r2_unpaired=str(OUT / "trimmomatic/{sample}_uR2.fastq"),
    output:
        str(OUT / "trimmomatic/{sample}_unpaired_joined.fastq")
    shell:
        """
        cat {input.r1_unpaired} {input.r2_unpaired} > {output}
        """


    #############################################################################
    ##### De novo assembly                                                  #####
    #############################################################################

rule run_SPAdes:
    input:
        r1=str(OUT / "trimmomatic/{sample}_pR1.fastq"),        
        r2=str(OUT / "trimmomatic/{sample}_pR2.fastq"),
        fastq_unpaired=str(OUT / "trimmomatic/{sample}_unpaired_joined.fastq")
    output:
        all_scaffolds=str(OUT / "SPAdes/{sample}/scaffolds.fasta"),
        filt_scaffolds=str(OUT / "scaffolds_filtered/{sample}_scaffolds_ge500nt.fasta")
    conda:
        "environments/de_novo_assembly.yaml"
    benchmark:
        str(OUT / "log/benchmark/De_novo_assembly_{sample}.txt")
    threads: config["threads"]["De_novo_assembly"]
    params:
        output_dir = str(OUT / "SPAdes/{sample}"),
        max_GB_RAM="100",
        kmersizes=config["SPAdes"]["kmersizes"],
        minlength=config["scaffold_minLen_filter"]["minlen"],
    log:
        str(OUT / "log/spades/{sample}_SPAdes_assembly.log")
    shell:
        """
        spades.py --isolate\
            -1 {input.r1:q} -2 {input.r2:q} \
            -s {input.fastq_unpaired} \
            -o {params.output_dir:q} \
            -k {params.kmersizes} \
            -m {params.max_GB_RAM} > {log:q}
            seqtk seq {output.all_scaffolds} 2>> {log} |\
            gawk -F "_" '/^>/ {{if ($4 >= {params.minlength}) {{print $0; getline; print $0}};}}' 2>> {log} 1> {output.filt_scaffolds} 

        """
    
    
    #############################################################################
    ##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
    #############################################################################

rule run_QUAST:
    input:
        str(OUT / "SPAdes/{sample}/scaffolds.fasta")
    output:
        str(OUT / "QUAST/per_sample/{sample}/report.html")
    conda:
        "environments/QUAST.yaml"
    threads: 4
    params:
        output_dir = str(OUT / "QUAST/per_sample/{sample}"),
    log:
        str(OUT / "log/quast/{sample}_QUAST_quality.log")
    shell:
        """
        quast --threads {threads} {input} --output-dir {params.output_dir} > {log:q}
        """


rule run_QUAST_combined:
    input:
        expand(str(OUT / "scaffolds_filtered/{sample}_scaffolds_ge500nt.fasta"), sample=SAMPLES)
    output:
        str(OUT / "QUAST/combined/report.tsv")
    conda:
        "environments/QUAST.yaml"
    threads: 4
    params:
        output_dir = str(OUT / "QUAST/combined"),
    log:
        str(OUT / "log/quast/quast_combined_quality.log")
    shell:
        """
        quast --threads {threads} {input:q} --output-dir {params.output_dir:q} > {log:q}
        """

rule run_CheckM:
    input:
        expand(str(OUT / "SPAdes/{sample}/scaffolds.fasta"), sample=SAMPLES)
    output:
        str(OUT / "CheckM/per_sample/{sample}/CheckM_{sample}.tsv"),
    conda:
        "environments/CheckM.yaml"
    threads: 4
    params:
        input_dir=str(OUT / "SPAdes/{sample}/"),
        output_dir=str(OUT / "CheckM/per_sample/{sample}"),
        genus = lambda wildcards: SAMPLES[wildcards.sample][1],
    log:
        str(OUT / "log/checkm/run_CheckM_{sample}.log")
    benchmark:
        str(OUT / "log/benchmark/CheckM_{sample}.txt")
    shell:
        """
        checkm taxonomy_wf genus "{params.genus}" {params.input_dir} {params.output_dir} -t {threads} -x scaffolds.fasta > {output}
        mv {params.output_dir}/checkm.log {log}
        """
        
rule parse_CheckM:
    input:
        expand(str(OUT / "CheckM/per_sample/{sample}/CheckM_{sample}.tsv"), sample=SAMPLES)
    output:
        str(OUT / "CheckM/CheckM_combined/CheckM_report.tsv")
    threads: 1
    log:
        str(OUT / "log/checkm/CheckM_combined.log")
    script:
        "bin/parse_checkM.py"


rule Fragment_length_analysis:
    input:
        fasta=str(OUT / "scaffolds_filtered/{sample}_scaffolds_ge500nt.fasta"),
        pR1=str(OUT / "trimmomatic/{sample}_pR1.fastq"),
        pR2=str(OUT / "trimmomatic/{sample}_pR2.fastq"),
    output:
        bam=str(OUT / "scaffolds_filtered/{sample}_sorted.bam"),
        bam_bai=str(OUT / "scaffolds_filtered/{sample}_sorted.bam.bai"),
        txt=str(OUT / "scaffolds_filtered/{sample}_insert_size_metrics.txt"),
        pdf=str(OUT / "scaffolds_filtered/{sample}_insert_size_histogram.pdf")
    conda:
        "environments/scaffold_analyses.yaml"
    log:
        str(OUT / "log/Fragment_length_analysis/Fragment_length_analysis_{sample}.log")
    benchmark:
        str(OUT / "log/benchmark/Fragment_length_analysis_{sample}.txt")
    threads: config["threads"]["Fragment_length_analysis"]
    shell:
        """
bwa index {input.fasta} > {log} 2>&1
bwa mem -t {threads} {input.fasta} \
{input.pR1} \
{input.pR2} 2>> {log} |\
samtools view -@ {threads} -uS - 2>> {log} |\
samtools sort -@ {threads} - -o {output.bam} >> {log} 2>&1
samtools index -@ {threads} {output.bam} >> {log} 2>&1
picard -Dpicard.useLegacyParser=false CollectInsertSizeMetrics \
-I {output.bam} \
-O {output.txt} \
-H {output.pdf} >> {log} 2>&1
        """

rule Generate_contigs_metrics:
    input:
        bam=str(OUT / "scaffolds_filtered/{sample}_sorted.bam"),
        fasta=str(OUT / "scaffolds_filtered/{sample}_scaffolds_ge500nt.fasta"),
    output:
        summary=str(OUT / "bbtools_scaffolds/per_sample/{sample}_MinLenFiltSummary.tsv"),
        perScaffold=str(OUT / "bbtools_scaffolds/per_sample/{sample}_perMinLenFiltScaffold.tsv"),
    conda:
        "environments/scaffold_analyses.yaml"
    log:
        str(OUT / "log/contigs_metrics/Generate_contigs_metrics_{sample}.log")
    benchmark:
        str(OUT / "log/benchmark/Generate_contigs_metrics_{sample}.txt")
    threads: 1
    shell:
        """
pileup.sh in={input.bam} \
ref={input.fasta} \
out={output.perScaffold} \
secondary=f \
samstreamer=t > {log} 2>&1
cp {log} {output.summary}
        """

rule parse_bbtools:
    input:
        expand(str(OUT / "bbtools_scaffolds/per_sample/{sample}_perMinLenFiltScaffold.tsv"), sample=SAMPLES),
    output:
        str(OUT / "bbtools_scaffolds/bbtools_combined/bbtools_scaffolds.tsv"),
    threads: 1
    log:
        str(OUT / "log/contigs_metrics/Generate_contigs_metrics_combined.log")
    script:
        "bin/parse_bbtools.py"

rule parse_bbtools_summary:
    input:
        expand(str(OUT / "bbtools_scaffolds/per_sample/{sample}_MinLenFiltSummary.tsv"), sample=SAMPLES)
    output:
        str(OUT / "bbtools_scaffolds/bbtools_combined/bbtools_summary_report.tsv")
    threads: 1
    log:
        str(OUT / "log/contigs_metrics/Generate_contigs_metrics_combined.log")
    script:
        "bin/parse_bbtools_summary.py"


rule MultiQC_report:
    input:
        expand(str(OUT / "FastQC_pretrim/{sample}_{read}_fastqc.zip"), sample = SAMPLES, read = "R1 R2".split()),
        expand(str(OUT / "FastQC_posttrim/{sample}_{read}_fastqc.zip"), sample = SAMPLES, read = "pR1 pR2 uR1 uR2".split()),
        str( OUT / "QUAST/combined/report.tsv"),
        str( OUT / "CheckM/CheckM_combined/CheckM_report.tsv"),
        expand(str(OUT / "log/trimmomatic/Clean_the_data_{sample}.log"), sample = SAMPLES),
        expand(str(OUT / "scaffolds_filtered/{sample}_insert_size_metrics.txt"), sample = SAMPLES),
    output:
        str(OUT / "MultiQC/multiqc.html"),
    conda:
        "environments/MultiQC_report.yaml"
    benchmark:
        str(OUT / "log/benchmark/MultiQC_report.txt")
    threads: 1
    params:
        config_file="files/multiqc_config.yaml",
        output_dir=str(OUT / "MultiQC")
    log:
        str(OUT / "log/multiqc/MultiQC_report.log")
    shell:
        """
multiqc --force --config {params.config_file} \
-o {params.output_dir} -n multiqc.html {input} > {log} 2>&1
    """