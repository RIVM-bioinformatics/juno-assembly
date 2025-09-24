rule top_species_multireport:
    input:
        bracken=expand(
            OUT + "/identify_species/contigs/{sample}/{sample}_species_content.txt",
            sample=SAMPLES,
        ),
        skani=OUT + "/identify_species/skani_results.tsv",
    output:
        OUT + "/identify_species/top1_species_multireport.csv",
    message:
        "Generating multireport for species identification."
    log:
        OUT + "/log/identify_species/multireport.log",
    threads: config["threads"]["parsing"]
    resources:
        mem_gb=config["mem_gb"]["parsing"],
    shell:
        """
        python bin/make_summary_main_species.py \
            --input-bracken-files {input.bracken} \
            --input-skani-file {input.skani} \
            --output-multireport {output} > {log} 2>&1
        """