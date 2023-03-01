import pandas as pd

sample_dataframes: list[pd.DataFrame] = []
# loop over the bbtools files (perscaffold) and add them to sample_dataframes
for input_file in snakemake.input:  # type: ignore
    sample_name = snakemake.wildcards.sample  # type:ignore

    sample_dataframe = pd.read_csv(input_file, sep="\t")
    sample_dataframe.insert(0, "Sample", sample_name)

    sample_dataframes.append(sample_dataframe)


# write concat output to file
pd.concat(sample_dataframes).to_csv(snakemake.output[0], index=False)  # type: ignore
