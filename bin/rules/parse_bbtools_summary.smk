#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################

rule parse_bbtools_summary:
    input:
        expand(str(OUT / "bbtools_scaffolds/per_sample/{sample}_MinLenFiltSummary.tsv"), sample=SAMPLES)
    output:
        str(OUT / "bbtools_scaffolds/bbtools_combined/bbtools_summary_report.tsv")
    threads: config["threads"]["parsing"]
    resources: mem_mb=config["mem_mb"]["parsing"]
    log:
        str(OUT / "log/contigs_metrics/pileup_contig_metrics_combined.log")
    script:
        "../parse_bbtools_summary.py"
