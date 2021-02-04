#############################################################################
##### Data quality control and cleaning                                 #####
#############################################################################

rule fastqc_clean_data:
    input:
        str(OUT / "trimmomatic/{sample}_{read}.fastq.gz")
    output:
        html=str(OUT / "FastQC_posttrim/{sample}_{read}_fastqc.html"),
        zip=str(OUT / "FastQC_posttrim/{sample}_{read}_fastqc.zip")
    conda:
        "../../envs/fastqc_trimmomatic.yaml"
    benchmark:
        str(OUT / "log/benchmark/fastqc_clean_data_{sample}_{read}.txt")
    threads: config["threads"]["fastqc"]
    resources: mem_mb=config["mem_mb"]["fastqc"]
    log:
        str(OUT / "log/fastqc/fastqc_clean_data_{sample}_{read}.log")
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
