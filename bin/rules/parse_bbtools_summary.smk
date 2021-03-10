#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################

rule parse_bbtools_summary:
    input:
        expand(OUT + "/qc_de_novo_assembly/bbtools_scaffolds/per_sample/{sample}_MinLenFiltSummary.tsv", sample=SAMPLES)
    output:
        OUT + "/qc_de_novo_assembly/bbtools_scaffolds/bbtools_summary_report.tsv"
    threads: config["threads"]["parsing"]
    resources: mem_mb=config["mem_mb"]["parsing"]
    log:
        OUT + "/log/qc_de_novo_assembly/pileup_contig_metrics_combined.log"
    script:
        "../parse_bbtools_summary.py"
