rule identify_species:
    input: 
        OUT + "/de_novo_assembly_filtered/{sample}.fasta"
    output:
        kraken2_kreport = temp(OUT + '/identify_species/{sample}.kreport2'),
        bracken_s = OUT + '/identify_species/{sample}_species_content.txt',
        bracken_kreport = OUT + '/identify_species/{sample}_bracken_species.kreport2'
    log:
        OUT + '/log/kraken2/{sample}_kraken2.log'
    threads: 
        config["threads"]["kraken2"]
    conda:
        '../../envs/identify_species.yaml'
    container:
        'library://alesr13/default/kraken2_bracken:v2.1.2_v2.6.1'
    params:
        kraken_db=config["db_dir"]
    resources:
        mem_gb=config["mem_gb"]["kraken2"]
    shell:
        """
kraken2 --db {params.kraken_db} \
    --threads {threads} \
    --report {output.kraken2_kreport} \
    {input}  &> {log} 

bracken -d {params.kraken_db} \
    -i {output.kraken2_kreport} \
    -o {output.bracken_s} \
    -r 150 \
    -l S \
    -t 0  &>> {log} 

        """
