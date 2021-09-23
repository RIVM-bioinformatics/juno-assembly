#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################

rule checkm:
    input:
        assembly = OUT + "/de_novo_assembly/{sample}/scaffolds.fasta",
        genus_bracken = OUT + '/identify_species/{sample}_bracken_species.kreport2'
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
    genus=$(grep "\sG\s" {input.genus_bracken} | head -n 1 | cut -f6 | xargs)
    genus_capitalized=${{genus^}}
    checkm taxonomy_wf genus "${{genus_capitalized}}" \
        {params.input_dir} \
        {params.output_dir} \
        -t {threads} \
        -x scaffolds.fasta > {output.result}
    mv {params.output_dir}/checkm.log {log}
    echo "\nNOTE: No genus was provided and therefore the top 1 genus (${{genus}}) detected by Kraken2 + Bracken was used for reference." >> {log}
else
    genus_capitalized={params.genus}
    genus_capitalized=${{genus_capitalized^}}
    checkm taxonomy_wf genus "${{genus_capitalized}}" \
        {params.input_dir} \
        {params.output_dir} \
        -t {threads} \
        -x scaffolds.fasta > {output.result}
    mv {params.output_dir}/checkm.log {log}
fi
        """
