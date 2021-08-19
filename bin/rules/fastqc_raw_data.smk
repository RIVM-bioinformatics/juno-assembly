#############################################################################
##### Data quality control and cleaning                                 #####
#############################################################################

rule qc_raw_fastq:
    input:
        lambda wildcards: SAMPLES[wildcards.sample][wildcards.read],
    output:
        html = OUT + "/qc_raw_fastq/{sample}_{read}_fastqc.html",
        zip = OUT + "/qc_raw_fastq/{sample}_{read}_fastqc.zip"
    conda:
        "../../envs/fastqc_trimmomatic.yaml"
    container:
        "biocontainers/fastqc:v0.11.9_cv8"
    threads: config["threads"]["fastqc"]
    resources: mem_gb=config["mem_gb"]["fastqc"]
    log:
        OUT + "/log/qc_raw_fastq/qc_raw_fastq_{sample}_{read}.log"
    params:
        output_dir = OUT + "/qc_raw_fastq"
    shell:
        """
        bash bin/fastqc_wrapper.sh {input} {params.output_dir} {output.html} {output.zip} {log} 
        """
