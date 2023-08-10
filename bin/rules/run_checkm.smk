#############################################################################
##### Scaffold analyses: QUAST, CheckM, picard, bbmap and QC-metrics    #####
#############################################################################


rule select_genus_checkm:
    input:
        genus_bracken=OUT
        + "/identify_species/contigs/{sample}/{sample}_bracken_species.kreport2",
        list_accepted_genera="files/accepted_genera_checkm.txt",
    output:
        selected_genus=OUT
        + "/qc_de_novo_assembly/checkm/per_sample/{sample}/selected_genus.txt",
    message:
        "Selecting genus for CheckM for {wildcards.sample}."
    params:
        genus=lambda wildcards: SAMPLES[wildcards.sample]["genus"],
    log:
        OUT + "/log/qc_de_novo_assembly/select_genus_checkm_{sample}.log",
    shell:
        """
        python bin/select_genus_checkm.py \
        --genus {params.genus} \
        --bracken-output {input.genus_bracken} \
        --output {output.selected_genus} 2>&1>{log}
        """


rule checkm:
    input:
        assembly=OUT + "/de_novo_assembly/{sample}/scaffolds.fasta",
        selected_genus=OUT
        + "/qc_de_novo_assembly/checkm/per_sample/{sample}/selected_genus.txt",
    output:
        result=OUT
        + "/qc_de_novo_assembly/checkm/per_sample/{sample}/checkm_{sample}.tsv",
        tmp_dir1=temp(
            directory(OUT + "/qc_de_novo_assembly/checkm/per_sample/{sample}/bins")
        ),
        tmp_dir2=temp(
            directory(OUT + "/qc_de_novo_assembly/checkm/per_sample/{sample}/storage")
        ),
    message:
        "Running CheckM for {wildcards.sample}."
    conda:
        "../../envs/checkm.yaml"
    container:
        "docker://ghcr.io/boasvdp/checkm-genome:1.1.3"
    threads: config["threads"]["checkm"]
    resources:
        mem_gb=config["mem_gb"]["checkm"],
    params:
        input_dir=OUT + "/de_novo_assembly/{sample}/",
        output_dir=OUT + "/qc_de_novo_assembly/checkm/per_sample/{sample}",
    log:
        OUT + "/log/qc_de_novo_assembly/checkm_{sample}.log",
    shell:
        """
        if [ $(<{input.selected_genus}) == "NOT_SUPPORTED" ]
        then
            echo -e "This is a mock report, because this genus is not supported by CheckM.\nscaffolds 0 100 100" > {output.result}
            mkdir -p {output.tmp_dir1} {output.tmp_dir2}
        else
            checkm taxonomy_wf genus "$(<{input.selected_genus})" \
                {params.input_dir} \
                {params.output_dir} \
                -t {threads} \
                -x scaffolds.fasta > {output.result}
            mv {params.output_dir}/checkm.log {log}
        fi
        """
