#############################################################################
##### Data quality control and cleaning                                 #####
#############################################################################

rule QC_clean_data:
    input:
        str(OUT / "trimmomatic/{sample}_{read}.fastq")
    output:
        html=str(OUT / "FastQC_posttrim/{sample}_{read}_fastqc.html"),
        zip=str(OUT / "FastQC_posttrim/{sample}_{read}_fastqc.zip")
    conda:
        "../../environments/QC_and_clean.yaml"
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
