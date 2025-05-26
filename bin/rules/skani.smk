rule skani:
    input:
        expand(
            OUT + "/de_novo_assembly/{sample}/scaffolds.fasta",
            sample=SAMPLES,
        ),
    output:
        OUT + "/skani/skani_results.tsv",
    message:
        "Generating skani report."
    conda:
        "../../envs/skani.yaml"
    container:
        "docker://quay.io/biocontainers/skani:0.2.2--ha6fb395_2"
    threads: config["threads"]["skani"]
    resources:
        mem_gb=config["mem_gb"]["skani"],
    log:
        OUT + "/log/skani/skani_report.log",
    params:
        max_no_hits=config["skani_max_no_hits"],
        gtdb_db_dir=config["skani_gtdb_db_dir"],
    shell:
        """
skani {input} \
    --output {output} \
    -d {params.gtdb_db_dir} \
    --ci \
    -n {params.max_no_hits} 2>&1>{log}
        """
