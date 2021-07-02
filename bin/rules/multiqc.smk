#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################

rule multiqc:
    input:
        expand(OUT + "/qc_raw_fastq/{sample}_{read}_fastqc.zip", sample = SAMPLES, read = "R1 R2".split()),
        expand(OUT + "/qc_clean_fastq/{sample}_{read}_fastqc.zip", sample = SAMPLES, read = "pR1 pR2".split()),
        OUT + "/qc_de_novo_assembly/quast/report.tsv",
        OUT + "/qc_de_novo_assembly/checkm/checkm_report.tsv",
        expand(OUT + "/log/clean_fastq/clean_fastq_{sample}.log", sample = SAMPLES),
        expand(OUT + "/qc_de_novo_assembly/insert_size/{sample}_insert_size_metrics.txt", sample = SAMPLES)
    output:
        OUT + "/multiqc/multiqc.html",
    conda:
        "../../envs/multiqc.yaml"
    benchmark:
        OUT + "/log/benchmark/multiqc.txt"
    threads: config["threads"]["multiqc"]
    resources: mem_mb=config["mem_mb"]["multiqc"]
    params:
        config_file="files/multiqc_config.yaml",
        output_dir=OUT + "/multiqc"
    log:
        OUT + "/log/multiqc/multiqc.log"
    shell:
        """
multiqc --force --config {params.config_file} \
-o {params.output_dir} -n multiqc.html {input} > {log} 2>&1
    """
