# extract completeness, contamination and Strain heterogeneity from each checkM output
checkm_dict: dict[str, list[str]] = {}
with open(snakemake.output[0], "w") as outfile:  # type: ignore
    outfile.write("sample\tcompleteness\tcontamination\tstrain_heterogeneity\n")
    for input_file in snakemake.input:  # type: ignore
        # add a letter to every sample name so that multiqc will see the name as a string
        sample = snakemake.wildcards.sample + "L"  # type: ignore

        with open(input_file, "r") as infile:
            for line in infile:
                if "scaffolds " in line:
                    completeness = line.split()[-3]
                    contamination = line.split()[-2]
                    strain_heterogeneity = line.split()[-1]
                    outfile.write(
                        f"{sample}\t{completeness}\t{contamination}\t{strain_heterogeneity}\n"
                    )
