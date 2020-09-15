#############################################################################
##### Data quality control and cleaning                                 #####
#############################################################################

rule QC_raw_data:
    input:
        lambda wildcards: SAMPLES[wildcards.sample][0][wildcards.read],
    output:
        html=str(OUT / "FastQC_pretrim/{sample}_{read}_fastqc.html"),
        zip=str(OUT / "FastQC_pretrim/{sample}_{read}_fastqc.zip")
    conda:
        "../../environments/QC_and_clean.yaml"
    benchmark:
        str(OUT / "log/benchmark/QC_raw_data_{sample}_{read}.txt")
    threads: 1
    log:
        str(OUT / "log/fastqc/QC_raw_data_{sample}_{read}.log")
    params:
        output_dir=str(OUT / "FastQC_pretrim/")
    shell:
        """
        bash bin/fastqc_wrapper.sh {input} {params.output_dir} {output.html} {output.zip} {log} 
        """
