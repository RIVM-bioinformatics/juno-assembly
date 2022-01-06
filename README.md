[![C.I.](https://github.com/RIVM-bioinformatics/Juno_pipeline/actions/workflows/juno_assembly_test.yaml/badge.svg)](https://github.com/RIVM-bioinformatics/Juno_pipeline/actions/workflows/juno_assembly_test.yaml)

<div align="center">
    <h1>Juno-assembly</h1>
    <br />
    <h2>Pipeline to process bacterial raw sequencing data up to de-novo assembly and the accompanying statistics.</h2>
    <br />
    <img src="files/juno_assembly_lightbg.png" alt="pipeline logo">
</div>

## Pipeline information

* **Author(s):**            Alejandra Hernández Segura, Robert Verhagen and Ernst Hamer.
* **Organization:**         National Institute for Public Health and the Environment (RIVM)
* **Department:**           Centre for Research Infectious Diseases Diagnostics and Screening (IDS), Bacteriology (BPD)
* **Start date:**           01 - 04 - 2020
* **Commissioned by:**      Maaike van den Beld

## About this project

The goal of this pipeline is to generate assemblies from raw fastq files. The input of the pipeline is raw Illumina paired-end data (read length > 99 bp) in the form of two fastq files (with extension .fastq, .fastq.gz, .fq or .fq.gz), containing the forward and the reversed reads ('R1' and 'R2' must be part of the file name, respectively). On the basis of the generated genome assemblies, low quality and contaminated samples can be excluded for downstream analysis. __Note:__ The pipeline has been validated only in gastroenteric bacteria (_Salmonella_, _Shigella_, _Listeria_ and STEC) but it could theoretically be used in other genera/species.

The pipeline uses the following tools:
1. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (Andrews, 2010) is used to assess the quality of the raw Illumina reads
2. [FastP](https://github.com/OpenGene/fastp) (Chen, Zhou, Chen and Gu, 2018) is used to remove poor quality data and adapter sequences 
3. [Picard](https://broadinstitute.github.io/picard/) determines the library fragment lengths
4. The reads are assembled into scaffolds by [SPAdes](https://cab.spbu.ru/software/spades/) (Bankevich et al., 2012) by means of _de novo_ assembly of the genome. SPAdes uses k-mers for building an initial de Bruijn graph and on following stages it performs graph-theoretical operations to assemble the genome. Kmer sizes of 21, 33, 55, 77 and 99 were used. For _de novo_ assembly, SPAdes isolate mode is used 
5. [QUAST](http://quast.sourceforge.net/) (Gurevich, Saveliev, Vyahhi, & Tesler, 2013) is used to assess the quality of the filtered scaffolds. 
6. To assess the quality of the microbial genomes, [CheckM](https://ecogenomics.github.io/CheckM/) (Parks, Imelfort, Skennerton, Hugenholtz, & Tyson, 2015) is used. CheckM calculates scores for completeness, contamination and strain heterogeneity 
7. [Bbtools](https://jgi.doe.gov/data-and-tools/bbtools/) (Bushnell, 2014) is used to generate scaffold alignment metrics 
8. [MultiQC](https://multiqc.info/) (Ewels, Magnusson, Lundin, & Käller, 2016) is used to summarize analysis results and quality assessments in a single report for dynamic visualization.
9. [Kraken2](https://ccb.jhu.edu/software/kraken2/) and [Bracken](http://ccb.jhu.edu/software/bracken/) for identification of bacterial species.  

![Image of pipeline](https://github.com/RIVM-bioinformatics/Juno_pipeline/blob/master/files/DAG.svg)

## Prerequisities

* **Linux + conda** A Linux-like environment with at least 'miniconda' installed
* Preferentially **[Singularity](https://sylabs.io/guides/3.8/user-guide/)** installed and working. If you are not using singularity, you have to run the pipeline using the argument `--no-containers`.


## Installation

**IMPORTANT NOTE**: You need to have [Conda](https://docs.conda.io/en/latest/) installed and, preferentially, also [Singularity](https://sylabs.io/guides/3.8/user-guide/) so that every step is containerized and therefore more reproducible. There is an option to run without using singularity (option `--no-containers`) but Conda is mandatory. 

1. Clone the repository:

```
git clone https://github.com/RIVM-bioinformatics/Juno_pipeline.git
```
Alternatively, you can download it manually as a zip file (you will need to unzip it then). If you decide to download the zip only, the pipeline version will not be stored in the audit trail.

2. Enter the directory with the pipeline and install the pipeline:

```
bash install_juno_assembly.sh
```

## Parameters & Usage

### Command for help

* ```-h, --help``` Shows the help of the pipeline
* ```--help-genera``` Prints the genera accepted by this pipeline.

### Required parameters

* ```-i, --input``` Relative or absolute path to the input directory. It must contain all the raw reads (fastq) files for all samples to be processed (not in subfolders). The input files should have the extension '.fastq', '.fastq.gz', '.fq' or '.fq.gz'. Important, fastq files are expected to have at least 1000 lines (so not be empty), otherwise they will be ignored (with a warning) and not included in the pipeline results since they are probably failed samples.  

### Optional parameters

* `--genus` Genus of the samples to be analyzed. Only one name is allowed, so if multiple samples are included in the input directory, it will be assumed that they all have the same genus. If metadata is given, the genus in the metadata will overwrite the one given through this option. If no `--genus` and no `--metadata` are provided, then CheckM will not be run. Keep in mind that if that is the case, the coverage of the assembly and other stats calculated by CheckM will not be obtained.

* `-m --metadata` Relative or absolute path to a .csv file. If provided, it must contain at least one column with the 'sample' name (name of the file but removing \_S##\_R1.fastq.gz) and a column called 'genus' (mind the small letters). The genus provided will be used to choose the reference genome to analyze de QC of the _de novo_ assembly. Default is: None. An example metadata would be:

| __sample__ | __genus__ |
| :---: | :--- |
| sample1 | salmonella |

*Note:* The fastq files corresponding to this sample would probably be something like sample1_S1_R1_0001.fastq.gz and sample1_S1_R2_0001.fastq.gz.

* ```-o --output``` Relative or absolute path to the output directory. If non is given, an 'output' directory will be created in the current directory. 
* `-mpt --mean-quality-threshold` Phred score to be used as threshold for cleaning (filtering) fastq files. Default is 28.
* `-w --window-size` Window size to use for cleaning (filtering) fastq files. Default is 5.
* `-ml --minimum-length` Minimum length for fastq reads to be kept after trimming. Default is 50.
* `-k --kmer-size` Kmersizes to be used for the de novo assembly. Default is 21,33,55,77,99. Note that the kmer size must be given without commas (`-k 21 33 55 77 99`)
* `-cl --contig-length-threshold` Minimum length to filter the contigs generated by the de novo assembly. Default is 500. 
* `--no-containers` Use conda environments instead of containers. Default is to run in docker/singularity containers. Note that **YOU NEED TO HAVE SINGULARITY INSTALLED AND WORKING** to use the containers option.
* `-p --prefix`     Conda or singularity prefix. Basically a path to the place where you want to store the conda environments or the singularity images.
* ```-c --cores```  Number of cores to use. Default is 4 if running locally (--local) or 300 otherwise.
* ```-l --local```  If this flag is present, the pipeline will be run locally (not attempting to send the jobs to an HPC cluster**). The default is to assume that you are working on a cluster because the pipeline was developed in an environment where it is the case. **Note that currently only LSF clusters are supported.
* ```-q --queue```  Name of the queue that the job will be submitted to if working on a cluster. This argument will be ignored if working locally (`--local`). It defaults to 'bio'. 
* ```-n --dryrun```, ```-u --unlock``` and ```--rerunincomplete``` are all parameters passed to Snakemake. If you want the explanation of these parameters, please refer to the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/).
* `--snakemake-args` Extra arguments to be passed to snakemake API (https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html). Passed in arg=value form. Only use this if you are familiar with the pipeline and with which parameters cannot be accepted as extra argument. 

### The base command to run this program. 

Before running the pipeline you should always activate the conda environment.

```
conda activate juno_assembly
```

```
python juno_assembly.py -i [path/to/input/dir] -o [path/to/output/dir]
```

You can deactivate the environment once you are finished.

```
conda deactivate
```

### An example on how to run the pipeline.

```
conda activate juno_assembly
python juno_assembly.py -i my_data -o my_results --genus Salmonella --local
conda deactivate
```

## Explanation of the output

* **log:** Log files with output and error files from each Snakemake rule/step that is performed. 
* **audit_trail:** Information about the versions of software and databases used.
* **output per sample:** The pipeline will create one subfolder per each step performed. These subfolders will in turn contain the results per sample and/or one with results combined for all samples (if applicable). To understand the ouptut, please refer to the manuals of each individual tool. The generated subfolders are:
    - clean_fastq: contains the fastq files after filtering low quality reads and trimming ends and/or adapters. This is what you would use in downstream analyses that require fastq files as input.
    - de_novo_assembly: results from the assembly step without filtering out the small contigs.
    - de_novo_assembly_filtered: results from the assembly step containing only contigs larger than 500bp. This is what you would normally use in downstream analyses that require assembly (fasta files) as input.
    - multiqc: MultiQC report for all samples run.
    - qc_clean_fastq: results of the quality control for fastq files. Run after filtering and trimming.
    - qc_de_novo_assembly: results of the quality control of the assemblies. Includes results of different tools (refer to those tools for interpretation).
    - qc_raw_fastq: results of the quality control for fastq files. Run before filtering and trimming. 
    - identify_species: results of bracken for species identification. The main result is a file called top1_species_multireport.csv which summarizes the main species found in the sample.  
        
## Issues  

* The pipeline currently only supports LSF clusters.
* Any issue can be reported in the [Issues section](https://github.com/RIVM-bioinformatics/Juno_pipeline/issues) of this repository.

## Future ideas for this pipeline

* Add a tool to identify bacteria/species in the sample.

## License
This pipeline is licensed with an AGPL3 license. Detailed information can be found inside the 'LICENSE' file in this repository.

## Contact
* **Contact person:**       Alejandra Hernández Segura
* **Email**                 alejandra.hernandez.segura@rivm.nl
