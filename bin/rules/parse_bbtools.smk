#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################

rule parse_bbtools:
    input:
        expand(OUT + "/qc_de_novo_assembly/bbtools_scaffolds/per_sample/{sample}_perMinLenFiltScaffold.tsv", sample=SAMPLES),
    output:
        OUT + "/qc_de_novo_assembly/bbtools_scaffolds/bbtools_scaffolds.tsv",
    message: "Parsing the results of bbtools (pileup contig metrics)."
    threads: config["threads"]["parsing"]
    resources: mem_gb=config["mem_gb"]["parsing"]
    log:
        OUT + "/log/qc_de_novo_assembly/pileup_contig_metrics_combined.log"
    script:
        "../parse_bbtools.py"

