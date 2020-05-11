
def parse_bbtools_summary(input_bbtools, output_bbtools):
    summary_dict = {}

    for input_file in str(input_bbtools).split():
        # get the sample name from the file name
        sample_name = str(input_file).split("sample/")[1].split("_")[0]
        variable_name_list = []
        value_list = []

        # get the summary stats from each bbtools outputfile
        summary_file = open(input_file, "r")
        for line in summary_file:
            if ":" in line:  #every line with : contains some information
                variable_name_list.append(line.split(":")[0].strip())
                value_list.append(line.split(":")[1].strip())
        summary_dict[sample_name] = value_list

    # specify a file to write to
    outfile = open(str(output_bbtools),"w")

    # write headers to file
    outfile.write("Sample\t")
    for item in variable_name_list:
        outfile.write(item + "\t")
    outfile.write("\n")

    # write summary stats to file
    for key, value in summary_dict.items():
        outfile.write(key + "\t")
        for item in value:
            outfile.write(item + "\t")
        outfile.write("\n")
        
    outfile.close()

parse_bbtools_summary(snakemake.input, snakemake.output)