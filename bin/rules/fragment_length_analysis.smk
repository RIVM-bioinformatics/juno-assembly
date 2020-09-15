#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################

rule Fragment_length_analysis:
    input:
        fasta=str(OUT / "scaffolds_filtered/{sample}_scaffolds_ge500nt.fasta"),
        pR1=str(OUT / "trimmomatic/{sample}_pR1.fastq"),
        pR2=str(OUT / "trimmomatic/{sample}_pR2.fastq"),
    output:
        bam=str(OUT / "scaffolds_filtered/{sample}_sorted.bam"),
        bam_bai=str(OUT / "scaffolds_filtered/{sample}_sorted.bam.bai"),
        txt=str(OUT / "scaffolds_filtered/{sample}_insert_size_metrics.txt"),
        pdf=str(OUT / "scaffolds_filtered/{sample}_insert_size_histogram.pdf")
    conda:
        "../../environments/scaffold_analyses.yaml"
    log:
        str(OUT / "log/Fragment_length_analysis/Fragment_length_analysis_{sample}.log")
    benchmark:
        str(OUT / "log/benchmark/Fragment_length_analysis_{sample}.txt")
    threads: config["threads"]["Fragment_length_analysis"]
    shell:
        """
bwa index {input.fasta} > {log} 2>&1
bwa mem -t {threads} {input.fasta} \
{input.pR1} \
{input.pR2} 2>> {log} |\
samtools view -@ {threads} -uS - 2>> {log} |\
samtools sort -@ {threads} - -o {output.bam} >> {log} 2>&1
samtools index -@ {threads} {output.bam} >> {log} 2>&1
picard -Dpicard.useLegacyParser=false CollectInsertSizeMetrics \
-I {output.bam} \
-O {output.txt} \
-H {output.pdf} >> {log} 2>&1
        """
