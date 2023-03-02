#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################


rule run_quast_combined:
    input:
        expand(OUT + "/de_novo_assembly_filtered/{sample}.fasta", sample=SAMPLES),
    output:
        report=OUT + "/qc_de_novo_assembly/quast/report.tsv",
        quast=OUT + "/qc_de_novo_assembly/quast/transposed_report.tsv",
    message:
        "Running QUAST in all samples and creating multireport."
    conda:
        "../../envs/quast.yaml"
    container:
        "docker://quay.io/biocontainers/quast:5.2.0--py310pl5321hc8f18ef_2"
    threads: config["threads"]["quast"]
    resources:
        mem_gb=config["mem_gb"]["quast"],
    params:
        output_dir=OUT + "/qc_de_novo_assembly/quast",
    log:
        OUT + "/log/qc_de_novo_assembly/quast.log",
    shell:
        """
        quast.py --threads {threads} {input} --output-dir {params.output_dir} > {log}
        """
