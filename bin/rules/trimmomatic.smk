#############################################################################
##### Data quality control and cleaning                                 #####
#############################################################################

rule trimmomatic:
    input:
        lambda wildcards: (SAMPLES[wildcards.sample][0][i] for i in ("R1", "R2"))
    output:
        r1=str(OUT / "trimmomatic/{sample}_pR1.fastq.gz"),
        r2=str(OUT / "trimmomatic/{sample}_pR2.fastq.gz"),
        r1_unpaired=str(OUT / "trimmomatic/{sample}_uR1.fastq.gz"),
        r2_unpaired=str(OUT / "trimmomatic/{sample}_uR2.fastq.gz"),
    conda:
        "../../envs/fastqc_trimmomatic.yaml"
    benchmark:
        str(OUT / "log/benchmark/trimmomatic_{sample}.txt")
    threads: config["threads"]["trimmomatic"]
    resources: mem_mb=config["mem_mb"]["trimmomatic"]
    log:
        str(OUT / "log/trimmomatic/trimmomatic_{sample}.log")
    params:
        adapter_removal_config=config["Trimmomatic"]["adapter_removal_config"],
        quality_trimming_config=config["Trimmomatic"]["quality_trimming_config"],
        minimum_length_config=config["Trimmomatic"]["minimum_length_config"],
    shell:
        """
trimmomatic PE -threads {threads} \
{input[0]:q} {input[1]:q} \
{output.r1} {output.r1_unpaired} \
{output.r2} {output.r2_unpaired} \
{params.adapter_removal_config} \
{params.quality_trimming_config} \
{params.minimum_length_config} > {log} 2>&1
touch -r {output.r1} {output.r1_unpaired}
touch -r {output.r2} {output.r2_unpaired}
        """
# touch the output in the case 100% of the reads are trimmed and no file is created
