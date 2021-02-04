#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################

rule pileup_contig_metrics:
    input:
        bam=str(OUT / "scaffolds_filtered/{sample}_sorted.bam"),
        fasta=str(OUT / "scaffolds_filtered/{sample}.fasta"),
    output:
        summary=str(OUT / "bbtools_scaffolds/per_sample/{sample}_MinLenFiltSummary.tsv"),
        perScaffold=str(OUT / "bbtools_scaffolds/per_sample/{sample}_perMinLenFiltScaffold.tsv"),
    conda:
        "../../envs/scaffold_analyses.yaml"
    log:
        str(OUT / "log/contigs_metrics/pileup_contig_metrics_{sample}.log")
    benchmark:
        str(OUT / "log/benchmark/pileup_contig_metrics_{sample}.txt")
    threads: config["threads"]["pileup"]
    resources: mem_mb=config["mem_mb"]["pileup"]
    shell:
        """
pileup.sh in={input.bam} \
ref={input.fasta} \
out={output.perScaffold} \
secondary=f \
samstreamer=t > {log} 2>&1
cp {log} {output.summary}
        """

