rule clean_fastq:
    input: 
        lambda wildcards: (SAMPLES[wildcards.sample][i] for i in ['R1', 'R2'])
    output: 
        r1=OUT + "/clean_fastq/{sample}_pR1.fastq.gz",
        r2=OUT + "/clean_fastq/{sample}_pR2.fastq.gz",
        unpaired = OUT + "/clean_fastq/{sample}_unpaired_joined.fastq.gz",
        html = OUT + "/clean_fastq/{sample}_fastp.html",
        json = OUT + "/clean_fastq/{sample}_fastp.json"
    message: "Filtering low quality reads for {wildcards.sample}."
    conda:
        "../../envs/qc_and_clean.yaml"
    container:
        "docker://biocontainers/fastp:v0.20.1_cv1"
    threads: config["threads"]["trimmomatic"]
    resources: mem_gb=config["mem_gb"]["trimmomatic"]
    log:
        OUT + "/log/clean_fastq/clean_fastq_{sample}.log"
    params:
        mean_quality = config['mean_quality_threshold'],
        window_size = config['window_size'],
        min_length = config['min_read_length']
    shell:
        """
fastp --in1 {input[0]} \
    --in2 {input[1]} \
    --out1 {output.r1} \
    --out2 {output.r2} \
    --unpaired1 {output.unpaired} \
    --unpaired2 {output.unpaired} \
    --html {output.html} \
    --json {output.json} \
    --report_title "FastP report for sample {wildcards.sample}" \
    --detect_adapter_for_pe \
    --thread {threads} \
    --cut_right \
    --cut_window_size {params.window_size} \
    --cut_mean_quality {params.mean_quality} \
    --correction \
    --length_required {params.min_length} > {log} 2>&1
        """