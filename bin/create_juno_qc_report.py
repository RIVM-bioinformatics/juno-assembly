"""
Parse several QC outputs from the Juno pipeline and
create a csv or excel file output
"""

import sys

print("something")
import os
import pandas as pd
import glob
from pathlib import Path
from functools import reduce
import openpyxl

import json
from pandas import json_normalize

import numpy as np

print("done with imports")

# input and output directories
input_dir = "/data/BioGrid/singhsp/Juno_assembly"  # path given as -i in pipeline
output_dir = (
    "/data/BioGrid/singhsp/Juno_assembly/output"  # path given as -o in pipeline
)


def get_metadata(species: str) -> pd.DataFrame:
    """
    Get the metadata file and create a pandas dataframe
    """
    filepath = Path(f"{species}")
    df = pd.read_csv(filepath, usecols=["sample", "genus"])
    df["sample"] = df["sample"].astype(str)
    print(df.dtypes)
    return df


def get_phred_score(phred: str) -> pd.DataFrame:
    """
    Parses multiqc_data json file to get phred quality scores for each sample into a pandas dataframe
    """
    with open(f"{phred}") as f:
        phred_data = json.load(f)

    phred = phred_data["report_plot_data"]["fastqc_per_sequence_quality_scores_plot"][
        "datasets"
    ]

    phred2 = []
    for val in phred:
        for i in val:
            phred2.append(i)

    df = pd.DataFrame(phred2)
    parsed_df = (
        df.groupby("name").data.apply(lambda x: pd.DataFrame(x.values[0])).reset_index()
    )

    parsed2_df = parsed_df.drop("level_1", axis=1)
    parsed2_df.rename(columns={0: "phred", 1: "Sequences"}, inplace=True)
    parsed3_df = parsed2_df[
        parsed2_df["name"].str.contains("_pR") == False
    ]  # Remove reads containing _pR

    parsed4_df = (
        parsed3_df.groupby("name")["Sequences"].max().reset_index()
    )  # get phred score by taking the score associated with the highest number of sequences per read
    parsed5_df = parsed3_df.merge(parsed4_df)
    parsed5_df["name"] = parsed5_df["name"].apply(
        lambda x: x.split("_")[0]
    )  # name to match metadata

    parsed6_df = (
        parsed5_df.groupby("name")["phred"].mean().reset_index()
    )  # mean phred score of R1 and R2
    parsed6_df.rename(
        columns={"name": "sample"}, inplace=True
    )  # rename name to sample to match rest of outputs
    parsed6_df["sample"] = parsed6_df["sample"].astype(str)

    print(parsed6_df.dtypes)
    return parsed6_df


def get_sequence_len(seq_len: str) -> pd.DataFrame:
    """
    Get the average sequence length from multiqc_fastqc.txt
    """

    seq_len_filepath = Path(f"{seq_len}")
    seq_len_df = pd.read_csv(
        seq_len_filepath, sep="\t", usecols=["Sample", "avg_sequence_length"]
    )
    seq_len_df.rename(columns={"Sample": "sample"}, inplace=True)
    seq_len_df = seq_len_df[
        seq_len_df["sample"].str.contains("_pR") == False
    ]  # remove reads containing _pR

    seq_len_df["sample"] = seq_len_df["sample"].apply(
        lambda x: x.split("_")[0]
    )  # sample name to match metadata

    seq_len_df = (
        seq_len_df.groupby("sample")["avg_sequence_length"].mean().reset_index()
    )  # average sequence length of R1 and R2
    seq_len_df["avg_sequence_length"] = seq_len_df["avg_sequence_length"].round(
        decimals=0
    )
    seq_len_df["sample"] = seq_len_df["sample"].astype(str)
    print(seq_len_df.dtypes)

    return seq_len_df


def get_transposed_report(quast: str) -> pd.DataFrame:
    """
    Creates a dataframe with necessary columns from the quast transposed report
    """
    trans_filepath = Path(f"{quast}")
    # print(trans_filepath)
    trans_report_df = pd.read_csv(
        trans_filepath,
        sep="\t",
        usecols=["Assembly", "Total length", "# contigs", "N50", "GC (%)"],
    )
    trans_report_df.rename(columns={"Assembly": "sample"}, inplace=True)
    trans_report_df["sample"] = trans_report_df["sample"].astype(str)
    print(trans_report_df.dtypes)

    return trans_report_df


def get_checkm_report(checkm: str) -> pd.DataFrame:
    """
    Creates a dataframe with necessary parameters from the checkm report
    """

    checkm_filepath = Path(f"{checkm}")
    # print(checkm_filepath)
    checkm_report_df = pd.read_csv(
        checkm_filepath,
        sep="\t",
        index_col=False,
        usecols=["sample", "completeness", "contamination"],
    )
    checkm_report_df["sample"] = checkm_report_df["sample"].str[:-1]
    checkm_report_df["sample"] = checkm_report_df["sample"].astype(str)

    print(checkm_report_df.dtypes)

    return checkm_report_df


def get_bbtools_report(bbtools: str) -> pd.DataFrame:
    """
    Creates a dataframe with necessary parameters from the bbtools report
    """

    bbtools_filepath = Path(f"{bbtools}")
    # print(bbtools_filepath)
    bbtools_report_df = pd.read_csv(
        bbtools_filepath, sep="\t", usecols=["Sample", "Reads", "Average coverage"]
    )
    bbtools_report_df.rename(columns={"Sample": "sample"}, inplace=True)
    bbtools_report_df["sample"] = bbtools_report_df["sample"].astype(str)
    print(bbtools_report_df.dtypes)

    return bbtools_report_df


def exceed_max(threshold: int, actual: int) -> str:
    """
    Highlight cell if actual value exceeds the threshold value
    """
    if actual > threshold:
        return "color:black; background-color:palegreen"
    return ""


def too_low(threshold: int, actual: int) -> str:
    """
    Highlight cell if actual value is less than the threshold
    """

    if actual < threshold:
        return "color:black; background-color:palegreen"
    return ""


def get_range(min: int, max: int, actual: int) -> str:
    """
    Highlight cell if value falls outside of a range
    """

    if actual < min or actual > max:
        return "color:black; background-color:palegreen"
    return ""


# def bacteria_colour(row: pd.Series) -> list(str):

#     """
#     Tolerance values for each parameter and bacterial species
#     """

#     genus = str(row['genus']).lower().strip()
#     # print(genus)
#     if genus == 'salmonella':
#         return ["", "", too_low(30, row['phred']), too_low(150, row['avg_sequence_length']), exceed_max(300, row['# contigs']),
#         get_range(4.4, 5.8, row['Total length']), get_range(51.6, 52.3, row['GC (%)']), too_low(30000, row['N50']),
#         "", too_low(30, row['Average coverage']), too_low(96, row['completeness']), exceed_max(4, row['contamination'])]

#     elif genus == 'escherichia':
#         return ["", "", too_low(30, row['phred']), too_low(150, row['avg_sequence_length']), exceed_max(300, row['# contigs']),
#         get_range(4.4, 5.8, row['Total length']), get_range(51.6, 52.3, row['GC (%)']), too_low(30000, row['N50']),
#         "", too_low(30, row['Average coverage']), too_low(96, row['completeness']), exceed_max(4, row['contamination'])]

#     elif genus == 'streptococcus':
#         return ["", "", too_low(30, row['phred']), too_low(150, row['avg_sequence_length']), exceed_max(300, row['# contigs']),
#         get_range(4.4, 5.8, row['Total length']), get_range(51.6, 52.3, row['GC (%)']), too_low(30000, row['N50']),
#         "", too_low(30, row['Average coverage']), too_low(96, row['completeness']), exceed_max(4, row['contamination'])]

#     elif genus == 'shigella':
#         return ["", "", too_low(30, row['phred']), too_low(150, row['avg_sequence_length']), exceed_max(762, row['# contigs']),
#         get_range(4.21, 5.03, row['Total length']), get_range(50.3, 51, row['GC (%)']), too_low(17500, row['N50']),
#         "", too_low(30, row['Average coverage']), too_low(96, row['completeness']), exceed_max(4, row['contamination'])]


#     elif genus == 'listeria':
#         return ["", "", too_low(30, row['phred']), too_low(150, row['avg_sequence_length']), exceed_max(300, row['# contigs']),
#         get_range(2.7, 3.23, row['Total length']), get_range(37.6, 38.2, row['GC (%)']), too_low(30000, row['N50']),
#         "", too_low(30, row['Average coverage']), too_low(96, row['completeness']), exceed_max(4, row['contamination'])]

#     elif genus == 'campylobacter':
#         return ["", "", too_low(30, row['phred']), too_low(150, row['avg_sequence_length']), exceed_max(100, row['# contigs']),
#         get_range(1.5, 1.9, row['Total length']), get_range(29.5, 31.5, row['GC (%)']), too_low(30000, row['N50']),
#         "", too_low(30, row['Average coverage']), too_low(96, row['completeness']), exceed_max(4, row['contamination'])]

#     elif genus == 'yersinia':
#         return ["", "", too_low(30, row['phred']), too_low(150, row['avg_sequence_length']), exceed_max(250, row['# contigs']),
#         get_range(3.9, 5.2, row['Total length']), get_range(46.2, 48.8, row['GC (%)']), too_low(30000, row['N50']),
#         "", too_low(30, row['Average coverage']), too_low(96, row['completeness']), exceed_max(4, row['contamination'])]

#     return ["", "", "", "", "", "", "", "", "", "", "", ""]


# def highlight_dataframe(dataframe: pd.DataFrame):
#     """
#     Apply conditional formatting to dataframe
#     """
#     styler = dataframe.style.apply(bacteria_colour, axis=1)

#     return styler


def bp_mbp(dataframe: pd.DataFrame, col: str) -> pd.DataFrame:
    """
    Convert base pairs to mega base pairs
    """
    dataframe[col] = dataframe[col] / 1000000
    return dataframe


def compile_report(species, phred, seq_len, quast, bbtools, checkm) -> pd.DataFrame:
    """
    Creates dataframe of required QC parameters from metadata, multiqc, quast, checkm and bbtools reports
    """

    # gather dataframes
    df_meta = get_metadata(species)
    df_phred = get_phred_score(phred)
    df_seq_length = get_sequence_len(seq_len)
    df_quast = get_transposed_report(quast)
    df_bbtools = get_bbtools_report(bbtools)
    df_checkm = get_checkm_report(checkm)

    # make list of dataframes
    dfs = [df_meta, df_phred, df_seq_length, df_quast, df_bbtools, df_checkm]
    # dfs = [df_meta, df_seq_length, df_quast, df_bbtools, df_checkm]

    # join dataframes
    final_df = reduce(
        lambda left, right: pd.merge(left, right, on=["sample"], how="outer"), dfs
    )

    # apply function to convert base pairs to mega base pairs
    final_df = bp_mbp(final_df, "Total length")
    # final_df2 = highlight_dataframe(final_df)

    return final_df


def get_csv_report(dataframe: pd.DataFrame) -> None:
    """
    Create csv format report from dataframe
    """

    dataframe.to_csv("juno_out.csv", index=False)


def get_excel_report(dataframe: pd.DataFrame) -> None:
    """
    Creates an excel format report with conditional formatting
    """
    with pd.ExcelWriter(sys.argv[-1], engine="openpyxl") as writer:
        for genus in dataframe["genus"].unique():
            newdf = dataframe[dataframe["genus"] == genus]
            # newdf = highlight_dataframe(newdf)
            newdf.to_excel(writer, sheet_name=genus, index=False)


def main(species, phred, seq_len, quast, bbtools, checkm):
    df_report = compile_report(species, phred, seq_len, quast, bbtools, checkm)
    print("at least compiling worked")
    return get_excel_report(df_report)


# def main(species, seq_len, quast, bbtools, checkm):
#     df_report = compile_report(species, seq_len, quast, bbtools, checkm)
#     print("at least compiling worked")
#     return get_excel_report(df_report)

print("running main")
main(
    *sys.argv[1:-1]
)  # snakemake.input.species, snakemake.input.phred, snakemake.input.seq_len, snakemake.input.quast, snakemake.input.bbtools, snakemake.input.checkm)
# (sys.argv[1:-1])
