#############################################################################
##### Data quality control and cleaning                                 #####
#############################################################################

rule clean_fastq:
    input:
        lambda wildcards: (SAMPLES[wildcards.sample][0][i] for i in ("R1", "R2"))
    output:
        r1=OUT + "/clean_fastq/{sample}_pR1.fastq.gz",
        r2=OUT + "/clean_fastq/{sample}_pR2.fastq.gz",
        r1_unpaired=OUT + "/clean_fastq/{sample}_uR1.fastq.gz",
        r2_unpaired=OUT + "/clean_fastq/{sample}_uR2.fastq.gz",
        joined_unpaired=OUT + "/clean_fastq/{sample}_unpaired_joined.fastq.gz"
    conda:
        "../../envs/fastqc_trimmomatic.yaml"
    benchmark:
        OUT + "/log/benchmark/clean_fastq_{sample}.txt"
    threads: config["threads"]["trimmomatic"]
    resources: mem_mb=config["mem_mb"]["trimmomatic"]
    log:
        OUT + "/log/clean_fastq/clean_fastq_{sample}.log"
    params:
        adapter_removal_config=config["trimmomatic"]["adapter_removal_config"],
        quality_trimming_config=config["trimmomatic"]["quality_trimming_config"],
        minimum_length_config=config["trimmomatic"]["minimum_length_config"],
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

cat {output.r1_unpaired} {output.r2_unpaired} > {output.joined_unpaired}

        """
# touch the output in the case 100% of the reads are trimmed and no file is created
