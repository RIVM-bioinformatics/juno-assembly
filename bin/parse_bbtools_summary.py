import argparse


def parse_bbtools_summary(input_bbtools, output_bbtools):
    summary_dict = {}

    for input_file in input_bbtools:
        # get the sample name from the file name
        sample_name = input_file.split("sample/")[1].split("_")[0]
        variable_name_list = []
        value_list = []

        # get the summary stats from each bbtools outputfile
        summary_file = open(input_file, "r")
        for line in summary_file:
            if ":" in line:  # every line with : contains some information
                variable_name_list.append(line.split(":")[0].strip())
                value_list.append(line.split(":")[1].strip())
        summary_dict[sample_name] = value_list

    # specify a file to write to
    outfile = open(str(output_bbtools), "w")

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


def main(input, output):
    parse_bbtools_summary(input, output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", type=str, required=True, metavar="PATH", nargs="+"
    )
    parser.add_argument("-o", "--output", type=str, required=True, metavar="FILE")
    args = parser.parse_args()
    main(args.input, args.output)
