#############################################################################
##### De novo assembly                                                  #####
#############################################################################

rule run_SPAdes:
    input:
        r1=str(OUT / "trimmomatic/{sample}_pR1.fastq.gz"),        
        r2=str(OUT / "trimmomatic/{sample}_pR2.fastq.gz"),
        fastq_unpaired=str(OUT / "trimmomatic/{sample}_unpaired_joined.fastq.gz")
    output:
        all_scaffolds=str(OUT / "SPAdes/{sample}/scaffolds.fasta"),
        filt_scaffolds=str(OUT / "scaffolds_filtered/{sample}.fasta"),
        k21=temp(directory(str(OUT / "SPAdes/{sample}/K21"))),
        k33=temp(directory(str(OUT / "SPAdes/{sample}/K33"))),
        k55=temp(directory(str(OUT / "SPAdes/{sample}/K55"))),
        k77=temp(directory(str(OUT / "SPAdes/{sample}/K77"))),
        k99=temp(directory(str(OUT / "SPAdes/{sample}/K99"))),
        misc=temp(directory(str(OUT / "SPAdes/{sample}/misc"))),
        tmp=temp(directory(str(OUT / "SPAdes/{sample}/tmp"))),
        state=temp(directory(str(OUT / "SPAdes/{sample}/pipeline_state"))),
        sh=temp(str(OUT / "SPAdes/{sample}/run_spades.sh")),
        yaml=temp(str(OUT / "SPAdes/{sample}/run_spades.yaml")),
        fastg=temp(str(OUT / "SPAdes/{sample}/assembly_graph.fastg")),
        gfa=temp(str(OUT / "SPAdes/{sample}/assembly_graph_with_scaffolds.gfa")),
        before=temp(str(OUT / "SPAdes/{sample}/before_rr.fasta")),
        contigs=temp(str(OUT / "SPAdes/{sample}/contigs.fasta")),
        contpath=temp(str(OUT / "SPAdes/{sample}/contigs.paths")),
        ds=temp(str(OUT / "SPAdes/{sample}/dataset.info")),
        dsyaml=temp(str(OUT / "SPAdes/{sample}/input_dataset.yaml")),
        params=temp(str(OUT / "SPAdes/{sample}/params.txt")),
        scaffpath=temp(str(OUT / "SPAdes/{sample}/scaffolds.paths")),
        splog=temp(str(OUT / "SPAdes/{sample}/spades.log"))
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