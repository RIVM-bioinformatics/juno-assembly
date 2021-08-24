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
        "docker://biocontainers/fastqc:v0.11.9_cv8"
    threads: config["threads"]["fastqc"]
    resources: mem_gb=config["mem_gb"]["fastqc"]
    log:
        OUT + "/log/qc_raw_fastq/qc_raw_fastq_{sample}_{read}.log"
    params:
        output_dir = OUT + "/qc_raw_fastq"
    shell:
        """
fastqc --quiet --outdir {params.output_dir} {input} &> {log} 

if [ ! -f {output.html} ]
then
    find {params.output_dir} -type f -name {wildcards.sample}*{wildcards.read}*.html -exec mv {{}} {output.html} \;
fi

if [ ! -f {output.zip} ]
then
    find {params.output_dir} -type f -name {wildcards.sample}*{wildcards.read}*.zip -exec mv {{}} {output.zip} \;
fi
        """
