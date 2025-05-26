"""
Juno_assembly pipeline
Authors: Alejandra Hernandez-Segura
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infectieziekteonderzoek, Diagnostiek en Laboratorium
            Surveillance (IDS), Bacteriologie (BPD)     
Date: 18-08-2021   

Documentation: https://rivm-bioinformatics.github.io/ids_bacteriology_man/juno-assembly.html 
"""

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml
from juno_library import Pipeline

from version import __description__, __package_name__, __version__


def main() -> None:
    juno_assembly = JunoAssembly()
    juno_assembly.run()


def get_suppported_checkm_genera() -> list[str]:
    with open(
        Path(__file__).parent.joinpath("files", "accepted_genera_checkm.txt"), mode="r"
    ) as f:
        return [g.strip().lower() for g in f.readlines()]


@dataclass
class JunoAssembly(Pipeline):
    pipeline_name: str = __package_name__
    pipeline_version: str = __version__
    input_type: str = "fastq"
    supported_genera: list[str] = field(default_factory=get_suppported_checkm_genera)

    def _add_args_to_parser(self) -> None:
        super()._add_args_to_parser()
        supported_genera = self.supported_genera

        class HelpGeneraAction(argparse.BooleanOptionalAction):
            def __call__(self, *args, **kwargs) -> None:  # type: ignore
                print("\n".join(["The accepted genera are:"] + supported_genera))
                exit(0)

        self.parser.description = "Juno_assembly pipeline. Automated pipeline for pre-processing, QC and assembly of bacterial NGS sequencing data."

        self.add_argument(
            "--help-genera",
            action=HelpGeneraAction,
            help="Prints the genera accepted by this pipeline.",
        )
        self.add_argument(
            "-g",
            "--genus",
            type=str.lower,
            choices=self.supported_genera,
            default=None,
            metavar="GENUS",
            help="Genus of the samples to be analyzed. If metadata is given, the genus in the metadata will overwrite the one given through this option.",
        )
        self.add_argument(
            "-m",
            "--metadata",
            type=Path,
            default=None,
            metavar="FILE",
            dest="metadata_file",
            help="Relative or absolute path to a .csv file. If provided, it must contain at least one column with the 'Sample' name (name of the file but removing _R1.fastq.gz) and a column called 'Genus' (mind the capital in the first letter). The genus provided will be used to choose the reference genome to analyze de QC of the de novo assembly.",
        )
        self.add_argument(
            "-d",
            "--db-dir",
            type=Path,
            metavar="DIR",
            default="/mnt/db/juno/kraken2_db",
            help="Relative or absolute path to the Kraken2 database. Default: /mnt/db/juno/kraken2_db.",
        )
        self.add_argument(
            "-sdb",
            "--skani-gtdb-db-dir",
            type=Path,
            metavar="DIR",
            default="/mnt/db/skani_gtdb-r226/gtdb_skani_database_ani-version-r226",
            help="Relative or absolute path to the Skani GTDB database. Default: /mnt/db/skani_gtdb-r226/gtdb_skani_database_ani-version-r226.",
        )
        self.add_argument(
            "-mpt",
            "--mean-quality-threshold",
            type=int,
            metavar="INT",
            default=28,
            help="Phred score to be used as threshold for cleaning (filtering) fastq files.",
        )
        self.add_argument(
            "-ws",
            "--window-size",
            type=int,
            metavar="INT",
            default=5,
            help="Window size to use for cleaning (filtering) fastq files.",
        )
        self.add_argument(
            "-ml",
            "--minimum-length",
            type=int,
            metavar="INT",
            default=50,
            dest="min_read_length",
            help="Minimum length for fastq reads to be kept after trimming.",
        )
        self.add_argument(
            "-k",
            "--kmer-size",
            nargs="+",
            type=int,
            metavar="INT INT...",
            default=[21, 33, 55, 77, 99],
            help="Kmersizes to be used for the de novo assembly.",
        )
        self.add_argument(
            "-cc",
            "--cov-cutoff",
            type=str,
            metavar="STR/INT",
            default="calculate",
            help="SPAdes k-mer coverage cut-off to use. Can be calculate, off, or a specified integer. "
            '"Calculate" lets the script calculate a sample-specific value that works for most use cases.',
        )
        self.add_argument(
            "-cl",
            "--contig-length-threshold",
            type=str,
            metavar="INT",
            default=500,
            help="Minimum length to filter the contigs generated by the de novo assembly.",
        )
        self.add_argument(
            "-td",
            "--target-depth",
            type=int,
            metavar="INT",
            default=150,
            help="Target depth for subsampling prior to de novo assembly",
        )
        self.add_argument(
            "-sm",
            "--skani-max-no-hits",
            type=int,
            metavar="INT",
            default=2,
            dest="skani_max_no_hits",
            help="Maximum number of hits to report for each contig in the Skani step. Default is 2, change value for debugging or development only.",
        )

    def _parse_args(self) -> argparse.Namespace:
        args = super()._parse_args()

        self.db_dir: Path = args.db_dir.resolve()
        self.metadata_file: Optional[Path] = args.metadata_file
        self.genus: Optional[str] = args.genus

        self.mean_quality_threshold = args.mean_quality_threshold
        self.window_size = args.window_size
        self.min_read_length = args.min_read_length
        self.kmer_size = args.kmer_size
        self.cov_cutoff = args.cov_cutoff
        self.contig_length_threshold = args.contig_length_threshold
        self.target_depth = args.target_depth
        self.skani_max_no_hits = args.skani_max_no_hits
        self.skani_gtdb_db_dir = args.skani_gtdb_db_dir.resolve()

        return args

    def __validate_kraken2_db_dir(self) -> bool:
        hash_file_exists = self.db_dir.joinpath("hash.k2d").exists()
        opts_file_exists = self.db_dir.joinpath("opts.k2d").exists()
        taxo_file_exists = self.db_dir.joinpath("taxo.k2d").exists()
        if hash_file_exists and opts_file_exists and taxo_file_exists:
            return True
        else:
            raise ValueError(
                f"The provided path to the database for Kraken2 ({str(self.db_dir)}) does not contain the expected files. Please download it again!"
            )

    def __validate_skani_gtdb_db_dir(self) -> bool:
        skani_db_dir = self.skani_gtdb_db_dir
        if not skani_db_dir.exists() or not skani_db_dir.is_dir():
            raise ValueError(
                f"The provided path for skani GTDB database ({str(skani_db_dir)}) does not exist or is not a directory. Please check the path and try again!"
            )
        else:
            return True

    def update_sample_dict_with_metadata(self) -> None:
        self.get_metadata_from_csv_file(
            filepath=self.metadata_file, expected_colnames=["sample", "genus"]
        )
        for sample, properties in self.sample_dict.items():
            try:
                properties["genus"] = (
                    self.juno_metadata[sample]["genus"].strip().lower()
                )
            except (KeyError, TypeError, AttributeError):
                properties["genus"] = self.genus  # type: ignore

    def setup(self) -> None:
        super().setup()
        if self.snakemake_args["use_singularity"]:
            self.snakemake_args["singularity_args"] = " ".join(
                [
                    self.snakemake_args["singularity_args"],
                    f"--bind {self.db_dir}:{self.db_dir}",
                ]
            )

        self.update_sample_dict_with_metadata()
        with open(
            Path(__file__).parent.joinpath("config/pipeline_parameters.yaml")
        ) as f:
            parameters_dict = yaml.safe_load(f)
        self.snakemake_config.update(parameters_dict)

        self.user_parameters = {
            "input_dir": str(self.input_dir),
            "out": str(self.output_dir),
            "exclusion_file": str(self.exclusion_file),
            "genus": self.genus,
            "mean_quality_threshold": self.mean_quality_threshold,
            "window_size": self.window_size,
            "min_read_length": self.min_read_length,
            "kmer_size": self.kmer_size,
            "cov_cutoff": self.cov_cutoff,
            "contig_length_threshold": self.contig_length_threshold,
            "run_in_container": self.snakemake_args["use_singularity"],
            "db_dir": str(self.db_dir),
            "gtdb_db_dir": str(self.skani_gtdb_db_dir),
            "target_depth": self.target_depth,
            "max_no_hits": int(self.skani_max_no_hits),
        }

        if not self.dryrun or self.unlock:
            self.__validate_kraken2_db_dir()
            self.__validate_skani_gtdb_db_dir()


if __name__ == "__main__":
    main()
