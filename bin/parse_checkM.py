# extract completeness, contamination and Strain heterogeneity from each checkM output


def parse_checkM(input_checkm, output_checkm):
    checkm_dict = {}
    for input_file in str(input_checkm).split():
        dict_key = "_".join(
            input_file.split("/")[-1].split(".")[0].split("_")[1:]
        )  # get sample name

        # add a letter to every sample name so that multiqc will see the name as a string
        dict_key = dict_key + "L"

        infile = open(input_file, "r")
        for line in infile:
            if "scaffolds " in line:
                completeness = line.split()[-3]
                contamination = line.split()[-2]
                strain_heterogeneity = line.split()[-1]
                checkm_dict[dict_key] = [
                    completeness,
                    contamination,
                    strain_heterogeneity,
                ]
        infile.close()

    # write output to file used for multiQC input
    outfile = open(str(output_checkm), "w")

    # write headers to file
    outfile.write("sample\tcompleteness\tcontamination\tstrain_heterogeneity\n")

    # write CheckM stats to file
    for x, y in checkm_dict.items():
        outfile.write(x + "\t")
        for item in y:
            outfile.write(item + "\t")
        outfile.write("\n")
    outfile.close()


parse_checkM(snakemake.input, snakemake.output)
