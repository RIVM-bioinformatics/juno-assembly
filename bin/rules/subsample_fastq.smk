rule subsample_fastq:
    input:
        r1=OUT + "/clean_fastq/{sample}_pR1.fastq.gz",
        r2=OUT + "/clean_fastq/{sample}_pR2.fastq.gz",
    output:
        r1 = OUT + "/subsampled_fastq/{sample}_pR1.fastq.gz",
        r2 = OUT + "/subsampled_fastq/{sample}_pR2.fastq.gz",
    message:
        "Subsampling reads for {wildcards.sample}."
    conda:
        "../../envs/seqtk_mash.yaml"
    container:
        "docker://ghcr.io/boasvdp/seqtk_mash:0.0.1"
    threads: int(config["threads"]["trimmomatic"])
    resources:
        mem_gb=config["mem_gb"]["trimmomatic"],
    log:
        OUT + "/log/subsample_fastq/subsample_fastq_{sample}.log",
    params:
        target_depth=config["target_depth"],
        r1_tmp = OUT + "/subsampled_fastq/{sample}_pR1.fastq",
        r2_tmp = OUT + "/subsampled_fastq/{sample}_pR2.fastq",
    shell:
        """
python bin/subsample_reads.py --input {input.r1} {input.r2} \
    --output {params.r1_tmp} {params.r2_tmp} \
    --depth {params.target_depth} 2>&1>{log}
gzip {params.r1_tmp}
gzip {params.r2_tmp}
        """