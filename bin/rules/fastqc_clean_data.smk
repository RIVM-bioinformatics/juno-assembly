#############################################################################
##### Data quality control and cleaning                                 #####
#############################################################################

rule qc_clean_fastq:
    input:
        pr1 = OUT + "/clean_fastq/{sample}_pR1.fastq.gz",
        pr2 = OUT + "/clean_fastq/{sample}_pR2.fastq.gz"
    output:
        html1 = OUT + "/qc_clean_fastq/{sample}_pR1_fastqc.html",
        zip1 = OUT + "/qc_clean_fastq/{sample}_pR1_fastqc.zip",
        html2 = OUT + "/qc_clean_fastq/{sample}_pR2_fastqc.html",
        zip2 = OUT + "/qc_clean_fastq/{sample}_pR2_fastqc.zip"
    conda:
        "../../envs/fastqc_trimmomatic.yaml"
    threads: config["threads"]["fastqc"]
    resources: mem_gb=config["mem_gb"]["fastqc"]
    log:
        OUT + "/log/qc_clean_fastq/qc_clean_fastq_{sample}.log"
    params:
        output_dir = OUT + "/qc_clean_fastq/"
    shell:
        """
if [ -s "{input.pr1}" ] # If file exists and is NOT empty
then
    fastqc --quiet --outdir {params.output_dir} {input.pr1} > {log}
else
    touch {output.html1} {output.zip1}
fi
if [ -s "{input.pr1}" ] # If file exists and is NOT empty
then
    fastqc --quiet --outdir {params.output_dir} {input.pr2} >> {log}
else
    touch {output.html2} {output.zip2}
fi
    """
