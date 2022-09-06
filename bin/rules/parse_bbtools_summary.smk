#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################


rule parse_bbtools_summary:
    input:
        expand(
            OUT
            + "/qc_de_novo_assembly/bbtools_scaffolds/per_sample/{sample}_MinLenFiltSummary.tsv",
            sample=SAMPLES,
        ),
    output:
        OUT + "/qc_de_novo_assembly/bbtools_scaffolds/bbtools_summary_report.tsv",
    message:
        "Parsing the results of bbtools (pileup contig metrics) and making a multireport."
    threads: config["threads"]["parsing"]
    resources:
        mem_gb=config["mem_gb"]["parsing"],
    log:
        OUT + "/log/qc_de_novo_assembly/pileup_contig_metrics_combined.log",
    shell:
        "python bin/parse_bbtools_summary.py -i {input} -o {output} > {log}"
