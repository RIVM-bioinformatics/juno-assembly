#############################################################################
##### De novo assembly                                                  #####
#############################################################################

rule run_de_novo_assembly:
    input:
        r1=OUT + "/clean_fastq/{sample}_pR1.fastq.gz",        
        r2=OUT + "/clean_fastq/{sample}_pR2.fastq.gz",
        fastq_unpaired=OUT + "/clean_fastq/{sample}_unpaired_joined.fastq.gz"
    output:
        all_scaffolds = OUT + "/de_novo_assembly/{sample}/scaffolds.fasta",
        filt_scaffolds = OUT + "/de_novo_assembly_filtered/{sample}.fasta",
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
        contigs=temp(OUT + "/de_novo_assembly/{sample}/contigs.fasta"),
        contpath=temp(OUT + "/de_novo_assembly/{sample}/contigs.paths"),
        ds=temp(OUT + "/de_novo_assembly/{sample}/dataset.info"),
        dsyaml=temp(OUT + "/de_novo_assembly/{sample}/input_dataset.yaml"),
        params=temp(OUT + "/de_novo_assembly/{sample}/params.txt"),
        scaffpath=temp(OUT + "/de_novo_assembly/{sample}/scaffolds.paths"),
        splog=temp(OUT + "/de_novo_assembly/{sample}/spades.log")
    conda:
        "../../envs/spades.yaml"
    threads: config["threads"]["spades"],
    resources: mem_gb=config["mem_gb"]["spades"]
    params:
        output_dir = OUT + "/de_novo_assembly/{sample}",
        max_GB_RAM="100",
        kmersizes=config["spades"]["kmersizes"],
        minlen=config["scaffold_minLen_filter"]["minlen"],
    log:
        OUT + "/log/de_novo_assembly/{sample}_de_novo_assembly.log"
    shell:
        """
unpaired_file_size=$(gzip -l {input.fastq_unpaired} | awk 'NR==2 {{exit($2!=0)}}')

if [ ${{unpaired_file_size}} -gt 0 ];then
    spades.py --isolate --only-assembler\
        -1 {input.r1} \
        -2 {input.r2} \
        -s {input.fastq_unpaired} \
        -o {params.output_dir} \
        -k {params.kmersizes} \
        -m {resources.mem_gb} \
        -t {threads} > {log}
else
    spades.py --isolate --only-assembler\
        -1 {input.r1} \
        -2 {input.r2} \
        -o {params.output_dir} \
        -k {params.kmersizes} \
        -m {resources.mem_gb} \
        -t {threads} > {log}
fi
        
if [ $? -eq 0 ]; then
    seqtk seq {output.all_scaffolds} 2>> {log} | \
        gawk -F "_" '/^>/ {{if ($4 >= {params.minlen}) {{print $0; getline; print $0}};}}' 2>> {log} 1> {output.filt_scaffolds} 
else
    exit 1
fi
        """
