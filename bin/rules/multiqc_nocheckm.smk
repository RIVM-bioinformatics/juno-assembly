#############################################################################
##### Scaffold analyses: QUAST, picard, bbmap and QC-metrics    #####
#############################################################################

rule multiqc:
    input:
        expand(str(OUT / "FastQC_pretrim/{sample}_{read}_fastqc.zip"), sample = SAMPLES, read = "R1 R2".split()),
        expand(str(OUT / "FastQC_posttrim/{sample}_{read}_fastqc.zip"), sample = SAMPLES, read = "pR1 pR2 uR1 uR2".split()),
        str( OUT / "quast/report.tsv"),
        expand(str(OUT / "log/trimmomatic/trimmomatic_{sample}.log"), sample = SAMPLES),
        expand(str(OUT / "scaffolds_filtered/{sample}_insert_size_metrics.txt"), sample = SAMPLES),
    output:
        str(OUT / "multiqc/multiqc.html"),
    conda:
        "../../envs/multiqc.yaml"
    benchmark:
        str(OUT / "log/benchmark/multiqc.txt")
    threads: config["threads"]["multiqc"]
    resources: mem_mb=config["mem_mb"]["multiqc"]
    params:
        config_file="files/multiqc_config.yaml",
        output_dir=str(OUT / "multiqc")
    log:
        str(OUT / "log/multiqc/multiqc.log")
    shell:
        """
multiqc --force --config {params.config_file} \
-o {params.output_dir} -n multiqc.html {input} > {log} 2>&1
    """
