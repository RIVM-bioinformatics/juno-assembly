#############################################################################
##### De novo assembly                                                  #####
#############################################################################

rule de_novo_assembly:
    input:
        r1 = OUT + "/clean_fastq/{sample}_pR1.fastq.gz",        
        r2 = OUT + "/clean_fastq/{sample}_pR2.fastq.gz",
        fastq_unpaired = OUT + "/clean_fastq/{sample}_unpaired_joined.fastq.gz"
    output:
        scaffolds = OUT + "/de_novo_assembly/{sample}/scaffolds.fasta",
        contigs = temp(OUT + "/de_novo_assembly/{sample}/contigs.fasta"),
        k21 = temp(directory(OUT + "/de_novo_assembly/{sample}/K21")),
        k33 = temp(directory(OUT + "/de_novo_assembly/{sample}/K33")),
        k55 = temp(directory(OUT + "/de_novo_assembly/{sample}/K55")),
        k77 = temp(directory(OUT + "/de_novo_assembly/{sample}/K77")),
        k99 = temp(directory(OUT + "/de_novo_assembly/{sample}/K99")),
        misc = temp(directory(OUT + "/de_novo_assembly/{sample}/misc")),
        tmp = temp(directory(OUT + "/de_novo_assembly/{sample}/tmp")),
        state = temp(directory(OUT + "/de_novo_assembly/{sample}/pipeline_state")),
        sh=temp(OUT + "/de_novo_assembly/{sample}/run_spades.sh"),
        yaml=temp(OUT + "/de_novo_assembly/{sample}/run_spades.yaml"),
        fastg=temp(OUT + "/de_novo_assembly/{sample}/assembly_graph.fastg"),
        gfa=temp(OUT + "/de_novo_assembly/{sample}/assembly_graph_with_scaffolds.gfa"),
        before=temp(OUT + "/de_novo_assembly/{sample}/before_rr.fasta"),
        contpath=temp(OUT + "/de_novo_assembly/{sample}/contigs.paths"),
        ds=temp(OUT + "/de_novo_assembly/{sample}/dataset.info"),
        dsyaml=temp(OUT + "/de_novo_assembly/{sample}/input_dataset.yaml"),
        params=temp(OUT + "/de_novo_assembly/{sample}/params.txt"),
        scaffpath=temp(OUT + "/de_novo_assembly/{sample}/scaffolds.paths"),
        splog=temp(OUT + "/de_novo_assembly/{sample}/spades.log"),
        gfa_simplified=temp(OUT + "/de_novo_assembly/{sample}/assembly_graph_after_simplification.gfa")
    message: "Making de novo assembly for {wildcards.sample}."
    conda:
        "../../envs/spades.yaml"
    container:
        "docker://quay.io/biocontainers/spades:3.15.3--h95f258a_0"
    threads: config["threads"]["spades"],
    resources: mem_gb=config["mem_gb"]["spades"]
    params:
        output_dir = OUT + "/de_novo_assembly/{sample}",
        kmersizes = config["kmer_size"]
    log:
        OUT + "/log/de_novo_assembly/{sample}_de_novo_assembly.log"
    shell:
        """
unpaired_file_size=$(wc {input.fastq_unpaired} | awk '{{print $1}}')

if [ ${{unpaired_file_size}} -gt 0 ];then
    echo "Running spades without using unpaired reads to make scaffolds (there were no unpaired reads).\n" > {log}
    spades.py --isolate \
        -1 {input.r1} \
        -2 {input.r2} \
        -s {input.fastq_unpaired} \
        -o {params.output_dir} \
        -k {params.kmersizes} \
        -m {resources.mem_gb} \
        -t {threads} >> {log}
else
    echo "Running spades using unpaired reads to make scaffolds.\n" > {log}
    spades.py --isolate \
        -1 {input.r1} \
        -2 {input.r2} \
        -o {params.output_dir} \
        -k {params.kmersizes} \
        -m {resources.mem_gb} \
        -t {threads} >> {log}
fi

if [ -f {output.scaffolds} ]; then
    cp {output.contigs} {output.scaffolds}
fi
        """




rule filter_de_novo_assembly:
    input:
        OUT + "/de_novo_assembly/{sample}/scaffolds.fasta"
    output:
        OUT + "/de_novo_assembly_filtered/{sample}.fasta"
    message: "Filtering out small contigs for the de novo assembly of {wildcards.sample}."
    conda:
        "../../envs/spades.yaml"
    container:
        "library://alesr13/default/seqtk_gawk:v1.3.1"
    threads: config["threads"]["parsing"],
    resources: mem_gb=config["mem_gb"]["parsing"]
    params:
        minlen=config["contig_length_threshold"]
    log:
        OUT + "/log/de_novo_assembly/{sample}_filter_assembly.log"
    shell:
        """
seqtk seq {input} 2>> {log} | gawk -F "_" '/^>/ {{if ($4 >= {params.minlen}) {{print $0; getline; print $0}};}}' 2>> {log} 1> {output} 
        """