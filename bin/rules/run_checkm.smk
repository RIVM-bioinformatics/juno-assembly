#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################

rule run_CheckM:
    input:
        expand(str(OUT / "SPAdes/{sample}/scaffolds.fasta"), sample=SAMPLES)
    output:
        str(OUT / "CheckM/per_sample/{sample}/CheckM_{sample}.tsv"),
    conda:
        "../../environments/CheckM.yaml"
    threads: 4
    params:
        input_dir=str(OUT / "SPAdes/{sample}/"),
        output_dir=str(OUT / "CheckM/per_sample/{sample}"),
        genus = lambda wildcards: SAMPLES[wildcards.sample][1],
    log:
        str(OUT / "log/checkm/run_CheckM_{sample}.log")
    benchmark:
        str(OUT / "log/benchmark/CheckM_{sample}.txt")
    shell:
        """
        checkm taxonomy_wf genus "{params.genus}" {params.input_dir} {params.output_dir} -t {threads} -x scaffolds.fasta > {output}
        mv {params.output_dir}/checkm.log {log}
        """
