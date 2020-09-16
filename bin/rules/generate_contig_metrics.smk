#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################

rule Generate_contigs_metrics:
    input:
        bam=str(OUT / "scaffolds_filtered/{sample}_sorted.bam"),
        fasta=str(OUT / "scaffolds_filtered/{sample}.fasta"),
    output:
        summary=str(OUT / "bbtools_scaffolds/per_sample/{sample}_MinLenFiltSummary.tsv"),
        perScaffold=str(OUT / "bbtools_scaffolds/per_sample/{sample}_perMinLenFiltScaffold.tsv"),
    conda:
        "../../environments/scaffold_analyses.yaml"
    log:
        str(OUT / "log/contigs_metrics/Generate_contigs_metrics_{sample}.log")
    benchmark:
        str(OUT / "log/benchmark/Generate_contigs_metrics_{sample}.txt")
    threads: 1
    shell:
        """
pileup.sh in={input.bam} \
ref={input.fasta} \
out={output.perScaffold} \
secondary=f \
samstreamer=t > {log} 2>&1
cp {log} {output.summary}
        """

