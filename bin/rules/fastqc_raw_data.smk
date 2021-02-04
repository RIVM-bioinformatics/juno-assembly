#############################################################################
##### Data quality control and cleaning                                 #####
#############################################################################

rule fastqc_raw_data:
    input:
        lambda wildcards: SAMPLES[wildcards.sample][0][wildcards.read],
    output:
        html=str(OUT / "FastQC_pretrim/{sample}_{read}_fastqc.html"),
        zip=str(OUT / "FastQC_pretrim/{sample}_{read}_fastqc.zip")
    conda:
        "../../envs/fastqc_trimmomatic.yaml"
    benchmark:
        str(OUT / "log/benchmark/fastqc_raw_data_{sample}_{read}.txt")
    threads: config["threads"]["fastqc"]
    resources: mem_mb=config["mem_mb"]["fastqc"]
    log:
        str(OUT / "log/fastqc/fastqc_raw_data_{sample}_{read}.log")
    params:
        output_dir=str(OUT / "FastQC_pretrim/")
    shell:
        """
        bash bin/fastqc_wrapper.sh {input} {params.output_dir} {output.html} {output.zip} {log} 
        """
