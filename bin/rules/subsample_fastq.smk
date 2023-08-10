rule subsample_fastq:
    input:
        r1=OUT + "/clean_fastq/{sample}_pR1.fastq.gz",
        r2=OUT + "/clean_fastq/{sample}_pR2.fastq.gz",
    output:
        r1=OUT + "/subsampled_fastq/{sample}_pR1.fastq.gz",
        r2=OUT + "/subsampled_fastq/{sample}_pR2.fastq.gz",
        cov_cutoff_file=OUT + "/subsampling/{sample}.txt",
    message:
        "Subsampling reads for {wildcards.sample}."
    conda:
        "../../envs/seqtk_mash.yaml"
    container:
        "docker://ghcr.io/boasvdp/seqtk_mash:0.0.2"
    threads: int(config["threads"]["trimmomatic"])
    resources:
        mem_gb=config["mem_gb"]["trimmomatic"],
    log:
        OUT + "/log/subsample_fastq/subsample_fastq_{sample}.log",
    params:
        target_depth=config["target_depth"],
        cov_cutoff=config["cov_cutoff"],
    shell:
        """
python bin/subsample_reads.py --input {input.r1} {input.r2} \
    --output {output.r1} {output.r2} \
    --depth {params.target_depth} \
    --cov-cutoff-in {params.cov_cutoff} \
    --cov-cutoff-out {output.cov_cutoff_file} \
    --threads {threads} 2>&1>{log}
        """
