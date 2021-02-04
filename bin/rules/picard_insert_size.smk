#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################

rule picard_insert_size:
    input:
        fasta=str(OUT / "scaffolds_filtered/{sample}.fasta"),
        pR1=str(OUT / "trimmomatic/{sample}_pR1.fastq.gz"),
        pR2=str(OUT / "trimmomatic/{sample}_pR2.fastq.gz"),
    output:
        bam=temp(str(OUT / "scaffolds_filtered/{sample}_sorted.bam")),
        bam_bai=temp(str(OUT / "scaffolds_filtered/{sample}_sorted.bam.bai")),
        txt=str(OUT / "scaffolds_filtered/{sample}_insert_size_metrics.txt"),
        ann=temp(str(OUT / "scaffolds_filtered/{sample}.fasta.ann")),
        amb=temp(str(OUT / "scaffolds_filtered/{sample}.fasta.amb")),
        bwt=temp(str(OUT / "scaffolds_filtered/{sample}.fasta.bwt")),
        pac=temp(str(OUT / "scaffolds_filtered/{sample}.fasta.pac")),
        sa=temp(str(OUT / "scaffolds_filtered/{sample}.fasta.sa")),
        pdf=str(OUT / "scaffolds_filtered/{sample}_insert_size_histogram.pdf")
    conda:
        "../../envs/scaffold_analyses.yaml"
    log:
        str(OUT / "log/picard_insert_size/picard_insert_size_{sample}.log")
    benchmark:
        str(OUT / "log/benchmark/picard_insert_size_{sample}.txt")
    threads: config["threads"]["picard"],
    resources: mem_mb=config["mem_mb"]["picard"]
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
