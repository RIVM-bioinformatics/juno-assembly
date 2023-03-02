#############################################################################
##### Quality Control Report    #####
#############################################################################


rule create_juno_QC_report:
    input:
        species=OUT + "/identify_species/top1_species_multireport.csv",
        phred=OUT + "/multiqc/multiqc_data/multiqc_data.json",
        seq_len=OUT + "/multiqc/multiqc_data/multiqc_fastqc.txt",
        quast=OUT + "/qc_de_novo_assembly/quast/transposed_report.tsv",
        bbtools=OUT
        + "/qc_de_novo_assembly/bbtools_scaffolds/bbtools_summary_report.tsv",
        checkm=OUT + "/qc_de_novo_assembly/checkm/checkm_report.tsv",
    output:
        OUT + "/Juno_assembly_QC_report/QC_report.xlsx",
    message:
        "Creating Juno assembly QC report."
    threads: config["threads"]["parsing"]
    resources:
        mem_gb=config["mem_gb"]["parsing"],
    log:
        OUT + "/log/create_juno_QC_report.log",
    script:
        "../create_juno_qc_report.py"
