#############################################################################
##### De novo assembly                                                  #####
#############################################################################

rule run_SPAdes:
    input:
        r1=str(OUT / "trimmomatic/{sample}_pR1.fastq"),        
        r2=str(OUT / "trimmomatic/{sample}_pR2.fastq"),
        fastq_unpaired=str(OUT / "trimmomatic/{sample}_unpaired_joined.fastq")
    output:
        all_scaffolds=str(OUT / "SPAdes/{sample}/scaffolds.fasta"),
        filt_scaffolds=str(OUT / "scaffolds_filtered/{sample}_scaffolds_ge500nt.fasta")
    conda:
        "../../environments/de_novo_assembly.yaml"
    benchmark:
        str(OUT / "log/benchmark/De_novo_assembly_{sample}.txt")
    threads: config["threads"]["De_novo_assembly"]
    params:
        output_dir = str(OUT / "SPAdes/{sample}"),
        max_GB_RAM="100",
        kmersizes=config["SPAdes"]["kmersizes"],
        minlength=config["scaffold_minLen_filter"]["minlen"],
    log:
        str(OUT / "log/spades/{sample}_SPAdes_assembly.log")
    shell:
        """
        spades.py --isolate --only-assembler\
            -1 {input.r1:q} -2 {input.r2:q} \
            -s {input.fastq_unpaired} \
            -o {params.output_dir:q} \
            -k {params.kmersizes} \
            -m {params.max_GB_RAM} > {log:q}
            seqtk seq {output.all_scaffolds} 2>> {log} |\
            gawk -F "_" '/^>/ {{if ($4 >= {params.minlength}) {{print $0; getline; print $0}};}}' 2>> {log} 1> {output.filt_scaffolds} 
        """
