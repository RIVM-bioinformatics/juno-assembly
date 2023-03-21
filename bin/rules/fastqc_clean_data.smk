#############################################################################
##### Data quality control and cleaning                                 #####
#############################################################################


rule qc_clean_fastq:
    input:
        OUT + "/clean_fastq/{sample}_p{read}.fastq.gz",
    output:
        html=OUT + "/qc_clean_fastq/{sample}_p{read}_fastqc.html",
        zip=OUT + "/qc_clean_fastq/{sample}_p{read}_fastqc.zip",
    message:
        "Running FastQC after filtering/trimming {wildcards.sample}."
    conda:
        "../../envs/qc_and_clean.yaml"
    container:
        "docker://biocontainers/fastqc:v0.11.9_cv8"
    threads: int(config["threads"]["fastqc"])
    resources:
        mem_gb=config["mem_gb"]["fastqc"],
    log:
        OUT + "/log/qc_clean_fastq/qc_clean_fastq_{sample}_{read}.log",
    params:
        output_dir=OUT + "/qc_clean_fastq/",
    shell:
        """
        if [ -s {input} ]
        then
            fastqc --quiet --outdir {params.output_dir} {input} >> {log} 
        else  
            touch {output}
        fi
        """
