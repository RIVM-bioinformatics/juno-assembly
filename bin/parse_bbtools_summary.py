summary_dict = {}

variable_name_list: list[str] = []

for input_file in snakemake.input:  # type: ignore
    sample_name = snakemake.wildcards.sample  # type: ignore
    variable_name_list = []
    value_list = []
    with open(input_file, "r") as summary_file:
        for line in summary_file:
            if ":" in line:  # every line with : contains some information
                variable_name_list.append(line.split(":")[0].strip())
                value_list.append(line.split(":")[1].strip())
        summary_dict[sample_name] = value_list

# specify a file to write to
with open(snakemake.output[0], "w") as outfile:  # type: ignore
    # write headers to file. TODO: This does expect all files to have the same headers
    outfile.write("\t".join(["Sample"] + variable_name_list) + "\n")

    # write summary stats to file
    for sample, value_list in summary_dict.items():
        outfile.write("\t".join([sample] + value_list) + "\n")
