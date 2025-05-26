rule identify_species_reads:
    input:
        r1=OUT + "/clean_fastq/{sample}_pR1.fastq.gz",
        r2=OUT + "/clean_fastq/{sample}_pR2.fastq.gz",
    output:
        kraken2_kreport=temp(OUT + "/identify_species/reads/{sample}/{sample}.kreport2"),
        bracken_s=OUT + "/identify_species/reads/{sample}/{sample}_species_content.txt",
        bracken_kreport=OUT
        + "/identify_species/reads/{sample}/{sample}_bracken_species.kreport2",
    message:
        "Running species identification for {wildcards.sample}."
    log:
        OUT + "/log/identify_species/reads/{sample}.log",
    threads: config["threads"]["kraken2"]
    conda:
        "../../envs/identify_species.yaml"
    container:
        "library://alesr13/default/kraken2_bracken:v2.1.2_v2.6.1"
    params:
        kraken_db=config["db_dir"],
    resources:
        mem_gb=config["mem_gb"]["kraken2"],
    shell:
        """
        # Adding --confidence 0.05 causes 0 kmers assigned to species in some samples
        # That breaks bracken
        kraken2 --db {params.kraken_db} \
            --threads {threads} \
            --report {output.kraken2_kreport} \
            --gzip-compressed \
            --paired \
            --output - \
            {input.r1} {input.r2}  &> {log} 

        bracken -d {params.kraken_db} \
            -i {output.kraken2_kreport} \
            -o {output.bracken_s} \
            -r 150 \
            -l S \
            -t 0  &>> {log} 

        """


rule identify_species:
    input:
        OUT + "/de_novo_assembly_filtered/{sample}.fasta",
    output:
        kraken2_kreport=temp(
            OUT + "/identify_species/contigs/{sample}/{sample}.kreport2"
        ),
        bracken_s=OUT
        + "/identify_species/contigs/{sample}/{sample}_species_content.txt",
        bracken_kreport=OUT
        + "/identify_species/contigs/{sample}/{sample}_bracken_species.kreport2",
    message:
        "Running species identification for {wildcards.sample}."
    log:
        OUT + "/log/identify_species/contigs/{sample}.log",
    threads: config["threads"]["kraken2"]
    conda:
        "../../envs/identify_species.yaml"
    container:
        "library://alesr13/default/kraken2_bracken:v2.1.2_v2.6.1"
    params:
        kraken_db=config["db_dir"],
    resources:
        mem_gb=config["mem_gb"]["kraken2"],
    shell:
        """
        # Adding --confidence 0.05 causes 0 kmers assigned to species in some samples
        # That breaks bracken
        kraken2 --db {params.kraken_db} \
            --threads {threads} \
            --report {output.kraken2_kreport} \
            {input}  &> {log} 

        bracken -d {params.kraken_db} \
            -i {output.kraken2_kreport} \
            -o {output.bracken_s} \
            -r 150 \
            -l S \
            -t 0  &>> {log} 

        """


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
