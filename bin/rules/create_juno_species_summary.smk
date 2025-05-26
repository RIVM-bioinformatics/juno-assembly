rule top_species_multireport:
    input:
        expand(
            OUT + "/identify_species/contigs/{sample}/{sample}_species_content.txt",
            sample=SAMPLES,
        ),
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
        python bin/make_summary_main_species.py --input-files {input} \
                                                --output-multireport {output} > {log}
        """