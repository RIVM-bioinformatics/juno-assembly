def parse_checkM(input_checkm, output_checkm):
    #extract completeness, contamination and Strain heterogeneity from each checkM output
    checkm_dict = {}
    for input_file in str(input_checkm).split():
        dict_key = input_file.split("_")[-1].split(".")[0]
        infile = open(input_file, "r")
        for line in infile:
            if "scaffolds " in line:
                completeness = line.split()[-3]
                contamination = line.split()[-2]
                strain_heterogeneity = line.split()[-1]
                checkm_dict[dict_key] = [completeness, contamination, strain_heterogeneity]
        infile.close()
    
    
    #write output to file used for multiQC input
    outfile = open(str(output_checkm),"w")

    outfile.write("sample\tcompleteness\tcontamination\tstrain_heterogeneity\n")
    for x, y in checkm_dict.items():
        outfile.write(x +"\t")
        for item in y:
            outfile.write(item + "\t")
        outfile.write("\n") 
    outfile.close()
    
parse_checkM(snakemake.input, snakemake.output)