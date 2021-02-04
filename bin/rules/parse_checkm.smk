#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################

rule parse_checkm:
    input:
        expand(str(OUT / "checkm/per_sample/{sample}/checkm_{sample}.tsv"), sample=SAMPLES)
    output:
        str(OUT / "checkm/checkm_combined/checkm_report.tsv")
    threads: config["threads"]["parsing"]
    resources: mem_mb=config["mem_mb"]["parsing"]
    log:
        str(OUT / "log/checkm/checkm_combined.log")
    script:
        "../parse_checkM.py"

