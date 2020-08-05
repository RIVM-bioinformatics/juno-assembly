#############################################################################
##### Data quality control and cleaning                                 #####
#############################################################################

rule Clean_the_data:
    input:
        lambda wildcards: (SAMPLES[wildcards.sample][0][i] for i in ("R1", "R2"))
    output:
        r1=str(OUT / "trimmomatic/{sample}_pR1.fastq"),
        r2=str(OUT / "trimmomatic/{sample}_pR2.fastq"),
        r1_unpaired=str(OUT / "trimmomatic/{sample}_uR1.fastq"),
        r2_unpaired=str(OUT / "trimmomatic/{sample}_uR2.fastq"),
    conda:
        "../../environments/QC_and_clean.yaml"
    benchmark:
        str(OUT / "log/benchmark/Clean_the_data_{sample}.txt")
    threads: config["threads"]["Clean_the_data"]
    log:
        str(OUT / "log/trimmomatic/Clean_the_data_{sample}.log")
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
