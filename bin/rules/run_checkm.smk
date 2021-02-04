#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################

rule checkm:
    input:
        expand(str(OUT / "spades/{sample}/scaffolds.fasta"), sample=SAMPLES)
    output:
        result=str(OUT / "checkm/per_sample/{sample}/checkm_{sample}.tsv"),
        tmp_dir1=temp(directory(str(OUT / "checkm/per_sample/{sample}/bins"))),
        tmp_dir2=temp(directory(str(OUT / "checkm/per_sample/{sample}/storage")))
    conda:
        "../../envs/checkm.yaml"
    threads: config["threads"]["checkm"],
    resources: mem_mb=config["mem_mb"]["checkm"]
    params:
        input_dir=str(OUT / "spades/{sample}/"),
        output_dir=str(OUT / "checkm/per_sample/{sample}"),
        genus = lambda wildcards: SAMPLES[wildcards.sample][1],
    log:
        str(OUT / "log/checkm/checkm_{sample}.log")
    benchmark:
        str(OUT / "log/benchmark/checkm_{sample}.txt")
    shell:
        """
        checkm taxonomy_wf genus "{params.genus}" {params.input_dir} {params.output_dir} -t {threads} -x scaffolds.fasta > {output.result}
        mv {params.output_dir}/checkm.log {log}
        """
