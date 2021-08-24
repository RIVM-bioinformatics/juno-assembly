#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################

rule checkm:
    input:
        OUT + "/de_novo_assembly/{sample}/contigs.fasta"
        # OUT + "/de_novo_assembly/{sample}/scaffolds.fasta"
    output:
        result = OUT + "/qc_de_novo_assembly/checkm/per_sample/{sample}/checkm_{sample}.tsv",
        tmp_dir1 = temp(directory(OUT + "/qc_de_novo_assembly/checkm/per_sample/{sample}/bins")),
        tmp_dir2 = temp(directory(OUT + "/qc_de_novo_assembly/checkm/per_sample/{sample}/storage"))
    conda:
        "../../envs/checkm.yaml"
    container:
        "docker://quay.io/biocontainers/checkm-genome:1.1.3--py_1"
    threads: config["threads"]["checkm"],
    resources: mem_gb=config["mem_gb"]["checkm"]
    params:
        input_dir=OUT + "/de_novo_assembly/{sample}/",
        output_dir=OUT + "/qc_de_novo_assembly/checkm/per_sample/{sample}",
        genus = lambda wildcards: SAMPLES[wildcards.sample]['genus']
    log:
        OUT + "/log/qc_de_novo_assembly/checkm_{sample}.log"
    shell:
        """
if [ "{params.genus}" == "None" ];then
    touch {output.result}
    mkdir -p {output.tmp_dir1}
    mkdir -p {output.tmp_dir2}
    echo "No genus was provided and therefore checkm was skipped" > {log}
else
    genus_capitalized={params.genus}
    genus_capitalized=${{genus_capitalized^}}
    checkm taxonomy_wf genus "${{genus_capitalized}}" \
        {params.input_dir} \
        {params.output_dir} \
        -t {threads} \
        -x contigs.fasta > {output.result} #TODO: Change to scaffolds.fasta when bug with empty unpaired is solved
    mv {params.output_dir}/checkm.log {log}
fi
        """
