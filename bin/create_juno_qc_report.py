"""
Parse several QC outputs from the Juno pipeline and
create a csv or excel file output
"""

import sys

# Redirect all printing and exceptions to the snakemake log of this rule
sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # type: ignore

import pandas as pd
from functools import reduce
import json


def get_genus(species_csv: str) -> pd.DataFrame:
    """
    Get the metadata file and create a pandas dataframe
    """
    df = pd.read_csv(species_csv, usecols=["sample", "genus"])
    df["sample"] = df["sample"].astype(str)
    return df


def get_phred_score(phred_json: str) -> pd.DataFrame:
    """
    Parses multiqc_data json file to get phred quality scores for each sample into a pandas dataframe
    """
    with open(phred_json) as f:
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
    # Remove reads containing _pR
    parsed3_df = parsed2_df[~parsed2_df["name"].str.contains("_pR")]

    # get phred score by taking the score associated with the highest number of sequences per read
    parsed4_df = parsed3_df.groupby("name")["Sequences"].max().reset_index()
    parsed5_df = parsed3_df.merge(parsed4_df)
    parsed5_df["name"] = parsed5_df["name"].astype(str).apply(lambda x: x.split("_")[0])  # type: ignore

    # set mean phred score of R1 and R2
    parsed6_df = parsed5_df.groupby("name")["phred"].mean().reset_index()
    # rename name to sample to match rest of outputs
    parsed6_df.rename(columns={"name": "sample"}, inplace=True)
    parsed6_df["sample"] = parsed6_df["sample"].astype(str)

    return parsed6_df


def get_sequence_len(seq_len_tsv: str) -> pd.DataFrame:
    """
    Get the average sequence length from multiqc_fastqc.txt
    """
    df = pd.read_csv(seq_len_tsv, sep="\t", usecols=["Sample", "avg_sequence_length"])
    df.rename(columns={"Sample": "sample"}, inplace=True)
    # remove reads containing _pR
    df = df[~df["sample"].str.contains("_pR")]

    # sample name to match metadata
    df["sample"] = df["sample"].apply(lambda x: x.split("_")[0])  # type: ignore

    # average sequence length of R1 and R2
    df = df.groupby("sample")["avg_sequence_length"].mean().reset_index()
    df["avg_sequence_length"] = df["avg_sequence_length"].round(decimals=0)
    df["sample"] = df["sample"].astype(str)

    return df


def get_transposed_report(quast_tsv: str) -> pd.DataFrame:
    """
    Creates a dataframe with necessary columns from the quast transposed report
    """
    df = pd.read_csv(
        quast_tsv,
        sep="\t",
        usecols=["Assembly", "Total length", "# contigs", "N50", "GC (%)"],
    )
    df.rename(columns={"Assembly": "sample"}, inplace=True)
    df["sample"] = df["sample"].astype(str)

    return df


def get_checkm_report(checkm_tsv: str) -> pd.DataFrame:
    """
    Creates a dataframe with necessary parameters from the checkm report
    """
    df = pd.read_csv(
        checkm_tsv,
        sep="\t",
        index_col=False,
        usecols=["sample", "completeness", "contamination"],
    )
    df["sample"] = df["sample"].str[:-1]
    df["sample"] = df["sample"].astype(str)

    return df


def get_bbtools_report(bbtools_tsv: str) -> pd.DataFrame:
    """
    Creates a dataframe with necessary parameters from the bbtools report
    """
    df = pd.read_csv(
        bbtools_tsv, sep="\t", usecols=["Sample", "Reads", "Average coverage"]
    )
    df.rename(columns={"Sample": "sample"}, inplace=True)
    df["sample"] = df["sample"].astype(str)

    return df


def compile_report(
    species_csv: str,
    phred_json: str,
    seq_len_tsv: str,
    quast_tsv: str,
    bbtools_tsv: str,
    checkm_tsv: str,
) -> pd.DataFrame:
    """
    Creates dataframe of required QC parameters from metadata, multiqc, quast, checkm and bbtools reports
    """

    # gather dataframes
    df_meta = get_genus(species_csv)
    df_phred = get_phred_score(phred_json)
    df_seq_length = get_sequence_len(seq_len_tsv)
    df_quast = get_transposed_report(quast_tsv)
    df_bbtools = get_bbtools_report(bbtools_tsv)
    df_checkm = get_checkm_report(checkm_tsv)

    # make list of dataframes
    dfs = [df_meta, df_phred, df_seq_length, df_quast, df_bbtools, df_checkm]
    # dfs = [df_meta, df_seq_length, df_quast, df_bbtools, df_checkm]

    # join dataframes
    final_df = reduce(
        lambda left, right: pd.merge(left, right, on=["sample"], how="outer"), dfs
    )

    # Convert base pairs to mega base pairs
    final_df["Total length"] = final_df["Total length"] / 1_000_000

    return final_df


def write_excel_report(df: pd.DataFrame, outfile: str) -> None:
    """
    Creates an excel format report with conditional formatting
    """
    with pd.ExcelWriter(outfile, engine="openpyxl") as writer:
        for genus, genus_df in df.groupby("genus"):
            genus_df["genus"] = str(genus)
            genus_df.to_excel(writer, sheet_name=str(genus), index=False)


report_df = compile_report(
    snakemake.input.species,  # type: ignore
    snakemake.input.phred,  # type: ignore
    snakemake.input.seq_len,  # type: ignore
    snakemake.input.quast,  # type: ignore
    snakemake.input.bbtools,  # type: ignore
    snakemake.input.checkm,  # type: ignore
)
write_excel_report(report_df, snakemake.output[0])  # type: ignore
