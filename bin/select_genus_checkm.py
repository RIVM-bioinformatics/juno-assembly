#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import logging


def read_bracken_report(path_to_report: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path_to_report,
        sep="\t",
        header=None,
        names=["pct", "count", "count_unique", "rank", "taxid", "name"],
    )
    df["name"] = df["name"].str.strip()
    return df


def get_top_microbial_hit(bracken_result: pd.DataFrame) -> str:
    df_genera = bracken_result[bracken_result["rank"] == "G"]
    top_hit = (
        df_genera.sort_values("count", ascending=False).reset_index().loc[0, "name"]
    )
    if top_hit == "Homo":
        logging.warning(
            f"The top species is the Homo genus, indicating contamination with an eukaryote."
        )
        top_hit_microbial = (
            df_genera.sort_values("count", ascending=False).reset_index().loc[1, "name"]
        )
    else:
        top_hit_microbial = top_hit
    return top_hit_microbial


def check_if_top_hit_is_supported(
    top_hit: str, path_to_list_accepted_genera: Path
) -> str:
    with open(path_to_list_accepted_genera, "r") as f:
        lines = f.readlines()
    list_accepted_genera = [accepted_genus.strip() for accepted_genus in lines]
    if top_hit in list_accepted_genera:
        selected_genus = top_hit
    else:
        logging.warning(
            f"The selected species is not supported by this version of CheckM."
        )
        selected_genus = "NOT_SUPPORTED"
    return selected_genus


def save_selected_genus(genus_name: str, output_path: Path) -> None:
    with open(output_path, "w+") as file:
        file.write(genus_name)


def main(args):
    # str "None" is provided on command line
    if args.genus == "None":
        logging.warning(
            f"No genus was provided, this will be guessed from Kraken2+bracken analysis."
        )
        bracken_result = read_bracken_report(args.bracken_output)
        top_hit = get_top_microbial_hit(bracken_result)
    else:
        top_hit = args.genus
    selected_genus = check_if_top_hit_is_supported(top_hit, args.list_accepted_genera)
    save_selected_genus(selected_genus, args.output)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--genus",
        type=str.capitalize,
        help="Genus supplied through metadata, overwriting Kraken2 analysis",
        default="None",
    )
    parser.add_argument(
        "--bracken-output", type=Path, help="Path to bracken output report"
    )
    parser.add_argument(
        "--list-accepted-genera",
        default=Path("files/accepted_genera_checkm.txt"),
        type=Path,
    )
    parser.add_argument(
        "--output", type=Path, help="Path to output file", required=True
    )

    args = parser.parse_args()

    main(args)
