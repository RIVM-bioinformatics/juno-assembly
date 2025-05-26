import argparse
import pathlib
import sys
from typing import Tuple

import pandas as pd

argument_parser = argparse.ArgumentParser(
    description="Make bracken results multireport with top 1 species (species with higher score per sample)."
)
argument_parser.add_argument(
    "-i",
    "--input-dir",
    type=pathlib.Path,
    default=None,
    required="-f" not in sys.argv and "--input-files" not in sys.argv,
    help="Directory where to find the individual bracken results (<sample>_species_content.txt files).",
)
argument_parser.add_argument(
    "-f",
    "--input-files",
    type=pathlib.Path,
    default=[],
    nargs="+",
    required="-i" not in sys.argv and "--input-dir" not in sys.argv,
    help="List of bracken results files (<sample>_species_content.txt files).",
)
argument_parser.add_argument(
    "-o",
    "--output-multireport",
    type=str,
    default="top1_species_multireport.csv",
    help="Path and name of the output file (must have .csv extension) for the desired multireport.",
)
args = argument_parser.parse_args()


class BrackenResult:
    """
    Class that reads the bracken (species_content.txt) result file per sample
    and can extract the top species hit
    """

    def __init__(self, filepath: pathlib.Path) -> None:
        self.filepath = pathlib.Path(filepath)
        assert self.filepath.exists(), f"The provided file {filepath} does not exist."

    def update_filepath(self, filepath: pathlib.Path) -> None:
        self.filepath = filepath

    def read_bracken_result(self) -> Tuple[str, pd.DataFrame]:
        """Read the <sample>_species_content.txt produced by Bracken"""
        bracken_result = pd.read_csv(self.filepath, sep="\t")
        sample_name = str(self.filepath.name)
        sample_name = sample_name.replace("_species_content.txt", "")
        return sample_name, bracken_result

    def find_top_hit(self) -> pd.DataFrame:
        print(f"Finding top species for {self.filepath}")
        sample_name, bracken_result = self.read_bracken_result()
        idx_top_species = bracken_result["fraction_total_reads"].argmax()
        bracken_result = bracken_result.loc[
            [idx_top_species], ["name", "taxonomy_id", "fraction_total_reads"]
        ]  # type: ignore
        bracken_result.rename(columns={"name": "full_species_name"}, inplace=True)
        bracken_result.insert(0, "sample", [sample_name])
        full_species_name = bracken_result.full_species_name.str.lower().str.split()
        bracken_result.insert(2, "genus", full_species_name.str.get(0))
        bracken_result.insert(3, "species", full_species_name.str.get(1))
        return bracken_result


input_dir = args.input_dir
input_files = args.input_files
output_multireport = args.output_multireport

if input_dir is not None:
    input_dir = pathlib.Path(input_dir)
    assert (
        input_dir.exists()
    ), f"The provided input directory {input_dir} does not exist."
    if len(input_files) == 0:
        input_files = list(input_dir.glob("*_species_content.txt"))
else:
    assert (
        len(input_files) > 0
    ), "You need to provide either an input_dir or a list of input_files to make a Bracken multireport."
    input_dir = None

assert output_multireport.endswith(
    ".csv"
), "The output_multireport must have a csv extension (now provided {output_multireport})."
output_multireport = pathlib.Path(output_multireport)


"""Read and concatenate the top 'n' result for multiple multireports"""
print("Creating multireport...")
top1_per_sample = [BrackenResult(file_).find_top_hit() for file_ in input_files]
report = pd.concat(top1_per_sample)

print(f"Writing multireport to file {output_multireport}...")
report.to_csv(output_multireport, index=False)
