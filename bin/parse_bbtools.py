import pandas


def parse_bbtools_perScaffold(input_bbtools, output_bbtools):
    # create an empty dataframe with the right headers
    bbtools_headers_file = open(input_bbtools[0], "r")
    bbtools_headers = bbtools_headers_file.readline().strip().split("\t")
    bbtools_headers.insert(0, "Sample")
    df_bbtools = pandas.DataFrame(columns=bbtools_headers)

    # loop over the bbtools files (perscaffold) and add them to the dataframe
    for input_file in str(input_bbtools).split():
        # get the sample name from the file name
        sample_name = str(input_file).split("sample/")[1].split("_")[0]

        # read the data into a pandas dataframe
        sample_dataframe = pandas.read_csv(input_file, sep="\t")

        # add sample name to the dataframe
        sample_dataframe.insert(0, "Sample", "")
        sample_dataframe["Sample"] = sample_name

        # append sample dataframe to master dataframe
        df_bbtools = df_bbtools.append(sample_dataframe, ignore_index=True)

    # write concat output to file
    df_bbtools.to_csv(str(output_bbtools), index=False)


parse_bbtools_perScaffold(snakemake.input, snakemake.output)
