#############################################################################
##### Data quality control and cleaning                                 #####
#############################################################################

rule clean_fastq:
    input:
        lambda wildcards: (SAMPLES[wildcards.sample][i] for i in ("R1", "R2"))
    output:
        r1=OUT + "/clean_fastq/{sample}_pR1.fastq.gz",
        r2=OUT + "/clean_fastq/{sample}_pR2.fastq.gz",
        r1_unpaired=OUT + "/clean_fastq/{sample}_uR1.fastq.gz",
        r2_unpaired=OUT + "/clean_fastq/{sample}_uR2.fastq.gz",
        joined_unpaired=OUT + "/clean_fastq/{sample}_unpaired_joined.fastq.gz"
    conda:
        "../../envs/fastqc_trimmomatic.yaml"
    container:
        "biocontainers/trimmomatic:v0.38dfsg-1-deb_cv1"
    threads: config["threads"]["trimmomatic"]
    resources: mem_gb=config["mem_gb"]["trimmomatic"]
    log:
        OUT + "/log/clean_fastq/clean_fastq_{sample}.log"
    params:
        adapter_removal_config=config["trimmomatic"]["adapter_removal_config"],
        min_qual=config["trimmomatic"]["quality_trimming_config"],
        minlen=config["trimmomatic"]["minimum_length_config"],
    shell:
        """
trimmomatic PE -threads {threads} \
    {input[0]} {input[1]} \
    {output.r1} {output.r1_unpaired} \
    {output.r2} {output.r2_unpaired} \
    {params.adapter_removal_config} \
    {params.min_qual} \
    {params.minlen} > {log}

# touch the output in the case 100% of the reads are trimmed and no file is created
touch -r {output.r1} {output.r1_unpaired}
touch -r {output.r2} {output.r2_unpaired}

cat {output.r1_unpaired} {output.r2_unpaired} > {output.joined_unpaired}
        """
