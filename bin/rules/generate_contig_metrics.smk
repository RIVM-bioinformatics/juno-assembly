#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################

rule pileup_contig_metrics:
    input:
        bam = OUT + "/qc_de_novo_assembly/insert_size/{sample}_sorted.bam",
        fasta = OUT + "/de_novo_assembly_filtered/{sample}.fasta"
    output:
        summary = OUT + "/qc_de_novo_assembly/bbtools_scaffolds/per_sample/{sample}_MinLenFiltSummary.tsv",
        perScaffold = OUT + "/qc_de_novo_assembly/bbtools_scaffolds/per_sample/{sample}_perMinLenFiltScaffold.tsv",
    conda:
        "../../envs/scaffold_analyses.yaml"
    log:
        OUT + "/log/qc_de_novo_assembly/pileup_contig_metrics_{sample}.log"
    benchmark:
        OUT + "/log/benchmark/pileup_contig_metrics_{sample}.txt"
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

