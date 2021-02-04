#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################


rule run_quast_combined:
    input:
        expand(str(OUT / "scaffolds_filtered/{sample}.fasta"), sample=SAMPLES)
    output:
        str(OUT / "quast/report.tsv")
    conda:
        "../../envs/quast.yaml"
    threads: config["threads"]["quast"],
    resources: mem_mb=config["mem_mb"]["quast"]
    params:
        output_dir = str(OUT / "quast"),
    log:
        str(OUT / "log/quast/quast_combined_quality.log")
    shell:
        """
        quast --threads {threads} {input:q} --output-dir {params.output_dir:q} > {log:q}
        """
