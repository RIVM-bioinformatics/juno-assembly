#############################################################################
##### Data quality control and cleaning                                 #####
#############################################################################

# Create a concatenated unpaired fastq file to support SPAdes de novo assembly
rule cat_unpaired:
    input:
        r1_unpaired=str(OUT / "trimmomatic/{sample}_uR1.fastq"),
        r2_unpaired=str(OUT / "trimmomatic/{sample}_uR2.fastq"),
    output:
        str(OUT / "trimmomatic/{sample}_unpaired_joined.fastq")
    shell:
        """
        cat {input.r1_unpaired} {input.r2_unpaired} > {output}
        """

