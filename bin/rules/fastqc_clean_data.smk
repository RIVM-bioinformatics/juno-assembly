#############################################################################
##### Data quality control and cleaning                                 #####
#############################################################################

rule qc_clean_fastq:
    input:
        OUT + "/clean_fastq/{sample}_{read}.fastq.gz"
    output:
        html = OUT + "/qc_clean_fastq/{sample}_{read}_fastqc.html",
        zip = OUT + "/qc_clean_fastq/{sample}_{read}_fastqc.zip"
    conda:
        "../../envs/fastqc_trimmomatic.yaml"
    benchmark:
        OUT + "/log/benchmark/qc_clean_fastq_{sample}_{read}.txt"
    threads: config["threads"]["fastqc"]
    resources: mem_mb=config["mem_mb"]["fastqc"]
    log:
        OUT + "/log/qc_clean_fastq/qc_clean_fastq_{sample}_{read}.log"
    params:
        output_dir = OUT + "/qc_clean_fastq/"
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
