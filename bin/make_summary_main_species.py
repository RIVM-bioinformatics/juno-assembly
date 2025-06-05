import argparse
import pathlib
import sys
from typing import Tuple, Self

import pandas as pd

argument_parser = argparse.ArgumentParser(description="Make bracken results multireport with top 1 species (species with higher score per sample).")
argument_parser.add_argument(
    "-i",
    "--input-bracken-dir",
    type=pathlib.Path,
    default=None,
    required="-f" not in sys.argv and "--input-bracken-files" not in sys.argv,
    help="Directory where to find the individual bracken results (<sample>_species_content.txt files).",
)
argument_parser.add_argument(
    "-f",
    "--input-bracken-files",
    type=pathlib.Path,
    default=[],
    nargs="+",
    required="-i" not in sys.argv and "--input-bracken-dir" not in sys.argv,
    help="List of bracken results files (<sample>_species_content.txt files).",
)
argument_parser.add_argument(
    "-s",
    "--input-skani-file",
    type=pathlib.Path,
    metavar="TSV",
    required=True,
    help="Skani species identification report.",
)
argument_parser.add_argument(
    "-o",
    "--output-multireport",
    type=pathlib.Path,
    metavar="CSV",
    default="top1_species_multireport.csv",
    help="Path and name of the output file (must have .csv extension) for the desired multireport. Default is '%(default)s'.",
)
argument_parser.add_argument(
    "--species-ani-cutoff",
    type=float,
    metavar="[0.0 - 100.0]",
    default=95.0,
    help="ANI cutoff value to consider a species as new or unclassified. Default is %(default)s%%.",
)
args = argument_parser.parse_args()


class BrackenResult:
    """
    Class that reads the bracken (species_content.txt) result file per sample
    and can extract the top species hit
    """

    def __init__(self: Self, filepath: pathlib.Path) -> None:
        self.filepath = pathlib.Path(filepath)
        assert self.filepath.exists(), f"The provided file {filepath} does not exist."

    def update_filepath(self: Self, filepath: pathlib.Path) -> None:
        self.filepath = filepath

    def read_bracken_result(self: Self) -> Tuple[str, pd.DataFrame]:
        """Read the <sample>_species_content.txt produced by Bracken"""
        bracken_result = pd.read_csv(self.filepath, sep="\t")
        sample_name = str(self.filepath.name)
        sample_name = sample_name.replace("_species_content.txt", "")
        return sample_name, bracken_result

    def find_top_hit(self: Self) -> pd.DataFrame:
        print(f"Finding top species for {self.filepath}")
        sample_name, bracken_result = self.read_bracken_result()
        idx_top_species = bracken_result["fraction_total_reads"].argmax()
        bracken_result = bracken_result.loc[[idx_top_species], ["name", "taxonomy_id", "fraction_total_reads"]]  # type: ignore
        bracken_result.rename(columns={"name": "full_species_name"}, inplace=True)
        bracken_result.insert(0, "sample", [sample_name])
        full_species_name = bracken_result.full_species_name.str.lower().str.split()
        bracken_result.insert(2, "genus", full_species_name.str.get(0))
        bracken_result.insert(3, "species", full_species_name.str.get(1))
        return bracken_result


def parse_skani_results(skani_output_path: pathlib.Path) -> pd.DataFrame:
    """
    Parse the skani output file and return a DataFrame with the top two hits per sample.
    Then it returns a wide format DataFrame with the top two hits for each sample (as
    requested by the Yersinia group). NB., the sample names are parsed from the query
    file as reported in the skani output file, so this assumes a pathing format as
    specified in the current smk: "*/[SAMPLE_NAME]/scaffolds.fasta".
    """

    column_format = {
        "ref_file": str,
        "query_file": str,
        "ani_value": float,
        "align_frac_ref": float,
        "align_frac_query": float,
        "ref_name": str,
        "query_name": str,
        "ani_lower_ci": float,
        "ani_upper_ci": float,
    }

    df = pd.read_csv(
        skani_output_path,
        sep="\t",
        header=0,
        names=column_format.keys(),  # type: ignore
        dtype=column_format,
    )

    #! parse the query_name from the query_file. NB., below assumes the suffix looks like */[SAMPLE_NAME]/scaffolds.fasta
    df["query_name"] = df["query_file"].apply(lambda x: x.split("/")[-2] if x.endswith("scaffolds.fasta") else None)  # type: ignore
    if df["query_name"].isnull().all():
        raise ValueError("The query_file does not have the expected format. It should end with 'scaffolds.fasta'.")

    df_sorted = df.sort_values(by=["query_name", "ani_value"], ascending=[True, False])
    top_hits = df_sorted.groupby("query_name").head(2).copy()  # ? get top two hits per sample only; as requested by Yersinia group
    top_hits["rank"] = top_hits.groupby("query_name").cumcount() + 1  # ? set rank 1 or 2 (0-based, so +1)

    top_hits_wide = top_hits.pivot(index="query_name", columns="rank")  # ? make wide format
    top_hits_wide.columns = [f"skani_{rank}_{col}" for col, rank in top_hits_wide.columns]  # ? flatten multi index columns

    column_order = [
        "skani_1_ref_file",
        "skani_1_query_file",
        "skani_1_ani_value",
        "skani_1_ani_lower_ci",
        "skani_1_ani_upper_ci",
        "skani_1_align_frac_ref",
        "skani_1_align_frac_query",
        "skani_1_ref_name",
        "skani_2_ref_file",
        "skani_2_query_file",
        "skani_2_ani_value",
        "skani_2_ani_lower_ci",
        "skani_2_ani_upper_ci",
        "skani_2_align_frac_ref",
        "skani_2_align_frac_query",
        "skani_2_ref_name",
    ]
    top_hits_wide = top_hits_wide[column_order].reset_index()  # ? reorder to match the expected output
    top_hits_wide = top_hits_wide.drop(columns=["skani_1_ref_file", "skani_1_query_file", "skani_2_ref_file", "skani_2_query_file"])

    return top_hits_wide


def throw_and_log_possible_new_species_warning(input_df: pd.DataFrame, species_cutoff: float = args.species_ani_cutoff) -> pd.DataFrame:
    """
    Takes as input the merged kraken/bracken and skani results DataFrame, it checks if the ANI value is below the species cutoff
    and if so, it throws a warning message and returns a DataFrame with the samples that have ANI values below the cutoff.
    """
    results = pd.DataFrame(columns=["sample", "skani_1_ani_value"])
    for sample, best_ani_value in zip(input_df["sample"], input_df["skani_1_ani_value"]):
        if best_ani_value < species_cutoff:
            print(
                f"Warning: The ANI value for sample {sample} is below the threshold of {species_cutoff}%. "
                "This may indicate a new or unclassified species."
            )
            new_row = pd.DataFrame([{"sample": sample, "skani_1_ani_value": best_ani_value}])
            results = pd.concat([results, new_row], ignore_index=True)
    return results


input_bracken_dir = args.input_bracken_dir
input_bracken_files = args.input_bracken_files
input_skani_file = args.input_skani_file

output_multireport = args.output_multireport

# ? Sanity checks
if input_bracken_dir is not None:
    input_bracken_dir = pathlib.Path(input_bracken_dir)
    assert input_bracken_dir.exists(), f"The provided input directory {input_bracken_dir} does not exist."
    if len(input_bracken_files) == 0:
        input_bracken_files = list(input_bracken_dir.glob("*_species_content.txt"))
else:
    assert len(input_bracken_files) > 0, "You need to provide either an input_dir or a list of input_files to make a Bracken multireport."
    input_bracken_dir = None

if input_skani_file.suffix != ".tsv":
    raise ValueError(f"The provided skani file {input_skani_file} must have a .tsv extension.")
input_skani_file = pathlib.Path(input_skani_file)
if not input_skani_file.exists():
    raise ValueError(f"The provided skani file {input_skani_file} must be a valid file path.")

if output_multireport.suffix != ".csv":
    raise ValueError(f"The output_multireport must have a csv extension (now provided {output_multireport}).")
output_multireport = pathlib.Path(output_multireport)


print("Reading kraken/bracken results...")
top1_per_sample = [BrackenResult(file_).find_top_hit() for file_ in input_bracken_files]
report = pd.concat(top1_per_sample)

print("Reading skani results...")
skani_results = parse_skani_results(input_skani_file)

# ? Merge the skani results with the bracken results
print("Merging skani results with bracken results...")
report = report.merge(
    skani_results,
    left_on="sample",
    right_on="query_name",
    how="outer",
)

# ? detect and flag possible new species based on ANI value
print(f"Checking for possible new species based on ANI value cutoff of {args.species_ani_cutoff}%...")
possible_new_species = throw_and_log_possible_new_species_warning(report, species_cutoff=args.species_ani_cutoff)

if not possible_new_species.empty:
    base_output_path = output_multireport.with_suffix("")
    possible_new_species_file = base_output_path.with_name(f"{base_output_path.name}_possible_new_species.csv")
    possible_new_species.to_csv(possible_new_species_file, index=False)
    print(f"\tPossible new species detected and saved to {possible_new_species_file}")
else:
    print("\tNo possible new species detected based on the provided ANI cutoff.")

print(f"Writing multireport to file {output_multireport}...")
report.to_csv(output_multireport, index=False)
