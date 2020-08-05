#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################

rule parse_CheckM:
    input:
        expand(str(OUT / "CheckM/per_sample/{sample}/CheckM_{sample}.tsv"), sample=SAMPLES)
    output:
        str(OUT / "CheckM/CheckM_combined/CheckM_report.tsv")
    threads: 1
    log:
        str(OUT / "log/checkm/CheckM_combined.log")
    script:
        "../parse_checkM.py"

