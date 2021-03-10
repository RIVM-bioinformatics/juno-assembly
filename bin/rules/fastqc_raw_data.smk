#############################################################################
##### Data quality control and cleaning                                 #####
#############################################################################

rule qc_raw_fastq:
    input:
        lambda wildcards: SAMPLES[wildcards.sample][0][wildcards.read],
    output:
        html = OUT + "/qc_raw_fastq/{sample}_{read}_fastqc.html",
        zip = OUT + "/qc_raw_fastq/{sample}_{read}_fastqc.zip"
    conda:
        "../../envs/fastqc_trimmomatic.yaml"
    benchmark:
        OUT + "/log/benchmark/qc_raw_fastq_{sample}_{read}.txt"
    threads: config["threads"]["fastqc"]
    resources: mem_mb=config["mem_mb"]["fastqc"]
    log:
        OUT + "/log/qc_raw_fastq/qc_raw_fastq_{sample}_{read}.log"
    params:
        output_dir = OUT + "/qc_raw_fastq"
    shell:
        """
        bash bin/fastqc_wrapper.sh {input} {params.output_dir} {output.html} {output.zip} {log} 
        """
