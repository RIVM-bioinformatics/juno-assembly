#############################################################################
##### Scaffold analyses: QUAST, picard, bbmap and QC-metrics    #####
#############################################################################

rule MultiQC_report:
    input:
        expand(str(OUT / "FastQC_pretrim/{sample}_{read}_fastqc.zip"), sample = SAMPLES, read = "R1 R2".split()),
        expand(str(OUT / "FastQC_posttrim/{sample}_{read}_fastqc.zip"), sample = SAMPLES, read = "pR1 pR2 uR1 uR2".split()),
        str( OUT / "QUAST/report.tsv"),
        expand(str(OUT / "log/trimmomatic/Clean_the_data_{sample}.log"), sample = SAMPLES),
        expand(str(OUT / "scaffolds_filtered/{sample}_insert_size_metrics.txt"), sample = SAMPLES),
    output:
        str(OUT / "MultiQC/multiqc.html"),
    conda:
        "../../environments/MultiQC_report.yaml"
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
