#############################################################################
##### Data quality control and cleaning                                 #####
#############################################################################

# Create a concatenated unpaired fastq file to support SPAdes de novo assembly
rule cat_unpaired:
    input:
        r1_unpaired=str(OUT / "trimmomatic/{sample}_uR1.fastq.gz"),
        r2_unpaired=str(OUT / "trimmomatic/{sample}_uR2.fastq.gz")
    threads: config["threads"]["parsing"]
    resources: mem_mb=config["mem_mb"]["parsing"]
    output:
        str(OUT / "trimmomatic/{sample}_unpaired_joined.fastq.gz")
    shell:
        """
        cat {input.r1_unpaired} {input.r2_unpaired} > {output}
        """

