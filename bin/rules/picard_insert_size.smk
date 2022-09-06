#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################


rule picard_insert_size:
    input:
        fasta=OUT + "/de_novo_assembly_filtered/{sample}.fasta",
        pR1=OUT + "/clean_fastq/{sample}_pR1.fastq.gz",
        pR2=OUT + "/clean_fastq/{sample}_pR2.fastq.gz",
    output:
        temp(
            multiext(
                OUT + "/de_novo_assembly_filtered/{sample}.fasta",
                ".ann",
                ".amb",
                ".bwt.2bit.64",
                ".pac",
                ".0123",
            )
        ),
        bam=temp(OUT + "/qc_de_novo_assembly/insert_size/{sample}_sorted.bam"),
        bam_bai=temp(OUT + "/qc_de_novo_assembly/insert_size/{sample}_sorted.bam.bai"),
        txt=OUT + "/qc_de_novo_assembly/insert_size/{sample}_insert_size_metrics.txt",
        pdf=OUT + "/qc_de_novo_assembly/insert_size/{sample}_insert_size_histogram.pdf",
    message:
        "Calculating insert size for {wildcards.sample}."
    conda:
        "../../envs/scaffold_analyses.yaml"
    container:
        "library://alesr13/default/samtools_bwa_picard:v0.1"
    log:
        OUT
        + "/log/qc_de_novo_assembly/picard_insert_size/picard_insert_size_{sample}.log",
    threads: config["threads"]["picard"]
    resources:
        mem_gb=config["mem_gb"]["picard"],
    params:
        run_in_container=config["run_in_container"],
    shell:
        """
        set -euo pipefail
        bwa-mem2 index {input.fasta} > {log} 2>&1

        bwa-mem2 mem -t {threads} {input.fasta} \
        {input.pR1} \
        {input.pR2} 2>> {log} |\
        samtools view -@ {threads} -uS - 2>> {log} |\
        samtools sort -@ {threads} - -o {output.bam} >> {log} 2>&1

        samtools index -@ {threads} {output.bam} >> {log} 2>&1

        if [ {params.run_in_container} == True ]; then
            picard-tools CollectInsertSizeMetrics \
            I={output.bam} \
            O={output.txt} \
            H={output.pdf} >> {log} 2>&1
        else
            picard -Dpicard.useLegacyParser=false CollectInsertSizeMetrics \
            -I {output.bam} \
            -O {output.txt} \
            -H {output.pdf} >> {log} 2>&1
        fi
        """
