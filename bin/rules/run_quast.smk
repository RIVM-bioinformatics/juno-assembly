#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################


rule run_quast_combined:
    input:
        expand(OUT + "/de_novo_assembly_filtered/{sample}.fasta", sample=SAMPLES)
    output:
        OUT + "/qc_de_novo_assembly/quast/report.tsv"
    conda:
        "../../envs/quast.yaml"
    threads: config["threads"]["quast"],
    resources: mem_gb=config["mem_gb"]["quast"]
    params:
        output_dir = OUT + "/qc_de_novo_assembly/quast"
    log:
        OUT + "/log/qc_de_novo_assembly/quast.log"
    shell:
        """
        quast --threads {threads} {input:q} --output-dir {params.output_dir:q} > {log:q}
        """
