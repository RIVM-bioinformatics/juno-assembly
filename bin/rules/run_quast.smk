#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################


rule run_QUAST_combined:
    input:
        expand(str(OUT / "scaffolds_filtered/{sample}.fasta"), sample=SAMPLES)
    output:
        str(OUT / "QUAST/report.tsv")
    conda:
        "../../environments/QUAST.yaml"
    threads: 4
    params:
        output_dir = str(OUT / "QUAST"),
    log:
        str(OUT / "log/quast/quast_combined_quality.log")
    shell:
        """
        quast --threads {threads} {input:q} --output-dir {params.output_dir:q} > {log:q}
        """
