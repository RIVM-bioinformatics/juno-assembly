"""
Juno_assembly pipeline
Authors: Alejandra Hernandez-Segura
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium
            Surveillance (IDS), Bacteriologie (BPD)     
Date: 18-08-2021   

Documentation: https://rivm-bioinformatics.github.io/ids_bacteriology_man/juno-assembly.html 
"""

from juno_library.helper_functions import error_formatter
import argparse
import os
from pathlib import Path
import sys
import yaml
from dataclasses import dataclass, field
from version import __package_name__, __version__, __description__
from typing import Any

from juno_library import PipelineStartup, RunSnakemake, helper_functions


@dataclass(kw_only=True)
class JunoAssemblyRun(PipelineStartup, RunSnakemake):
    """Class with the arguments and specifications that are only for the
    Juno_assembly pipeline but inherit from PipelineStartup and RunSnakemake
    """

    input_type: str = "fastq"
    pipeline_name: str = __package_name__
    pipeline_version: str = __version__
    workdir: Path = Path(__file__).parent.resolve()
    db_dir: Path = Path("/mnt/db/juno/kraken2_db")
    min_num_lines: int = 1000  # TODO: Find ideal min num of reads/lines needed
    name_snakemake_report: str = "juno_assembly_report.html"

    exclusion_file: None | Path = None
    genus: None | str = None
    metadata_file: None | Path = None
    mean_quality_threshold: int = 28
    window_size: int = 5
    min_read_length: int = 50
    kmer_size: list[int] = field(default_factory=lambda: [21, 33, 55, 77, 99])
    contig_length_threshold: int = 500
    help_genera: bool = False
    cores: int = 300
    time_limit: int = 60
    local: bool = False
    queue: str = "bio"
    unlock: bool = False
    rerunincomplete: bool = False
    dryrun: bool = False
    run_in_container: bool = False
    singularity_prefix: None | str = None

    def __post_init__(
        self,
        **kwargs: dict[str, Any],
    ) -> None:
        """Initiating Juno_assembly pipeline"""

        if self.help_genera:
            print("The accepted genera are:")
            os.system("cat files/accepted_genera_checkm.txt")
            sys.exit(0)

        self.singularityargs = f"--bind {self.input_dir}:{self.input_dir} --bind {self.output_dir}:{self.output_dir} --bind {self.db_dir}:{self.db_dir}"
        # Build class
        PipelineStartup.__post_init__(self)
        RunSnakemake.__post_init__(
            self,
            **kwargs,
        )

        # Specific for Juno assembly
        self.supported_genera = []
        with open(
            self.workdir.joinpath("files", "accepted_genera_checkm.txt")
        ) as file_:
            for line in file_:
                genus_name = line.replace("\n", "").lower()
                self.supported_genera.append(genus_name)

    def __check_genus_is_supported(self) -> bool:
        if self.genus is not None:
            if self.genus.strip().lower() in self.supported_genera:
                return True
            else:
                raise ValueError(
                    error_formatter(
                        f'The genus {self.genus} is not supported. You can leave the "genus" empty for samples with unsupported genera.'
                    )
                )
        else:
            return True

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

    def update_sample_dict_with_metadata(self) -> None:
        self.__check_genus_is_supported()
        self.get_metadata_from_csv_file(
            filepath=self.metadata_file, expected_colnames=["sample", "genus"]
        )
        for sample in self.sample_dict:
            try:
                self.sample_dict[sample]["genus"] = (
                    self.juno_metadata[sample]["genus"].strip().lower()
                )
            except (KeyError, TypeError):
                self.sample_dict[sample]["genus"] = self.genus

    def write_userparameters(self) -> dict[str, Any]:
        config_params = {
            "input_dir": str(self.input_dir),
            "out": str(self.output_dir),
            "exclusion_file": str(self.exclusion_file),
            "genus": self.genus,
            "mean_quality_threshold": self.mean_quality_threshold,
            "window_size": self.window_size,
            "min_read_length": self.min_read_length,
            "kmer_size": self.kmer_size,
            "contig_length_threshold": self.contig_length_threshold,
            "run_in_container": self.usesingularity,
            "db_dir": str(self.db_dir),
        }

        with open(self.user_parameters, "w") as file:
            yaml.dump(config_params, file, default_flow_style=False)

        return config_params

    def run_juno_assembly_pipeline(self) -> None:
        self.update_sample_dict_with_metadata()
        with open(self.sample_sheet, "w") as file:
            yaml.dump(self.sample_dict, file, default_flow_style=False)

        self.user_params = self.write_userparameters()
        if not self.dryrun or self.unlock:
            self.__validate_kraken2_db_dir()
        self.successful_run = self.run_snakemake()
        assert self.successful_run, f"Please check the log files"
        if not self.dryrun or self.unlock:
            self.make_snakemake_report()


def main() -> None:
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument(
        "--help-genera",
        action="store_true",
        help="Prints the genera accepted by this pipeline.",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=not "--help-genera" in sys.argv,
        metavar="DIR",
        help="Relative or absolute path to the input directory. It must contain all the raw reads (fastq) files for all samples to be processed (not in subfolders).",
    )
    parser.add_argument(
        "-g",
        "--genus",
        type=str.lower,
        default=None,
        metavar="GENUS",
        help="Genus of the samples to be analyzed. If metadata is given, the genus in the metadata will overwrite the one given through this option.",
    )
    parser.add_argument(
        "-m",
        "--metadata",
        type=Path,
        default=None,
        metavar="FILE",
        help="Relative or absolute path to a .csv file. If provided, it must contain at least one column with the 'Sample' name (name of the file but removing _R1.fastq.gz) and a column called 'Genus' (mind the capital in the first letter). The genus provided will be used to choose the reference genome to analyze de QC of the de novo assembly.",
    )
    parser.add_argument(
        "-ex",
        "--exclusionfile",
        type=Path,
        metavar="FILE",
        dest="exclusion_file",
        help="Path to the file that contains samplenames to be excluded.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        metavar="DIR",
        default="output",
        help="Relative or absolute path to the output directory. If non is given, an 'output' directory will be created in the current directory.",
    )
    parser.add_argument(
        "-d",
        "--db-dir",
        type=Path,
        metavar="DIR",
        default="/mnt/db/juno/kraken2_db",
        help="Relative or absolute path to the Kraken2 database. Default: /mnt/db/juno/kraken2_db.",
    )
    parser.add_argument(
        "-mpt",
        "--mean-quality-threshold",
        type=int,
        metavar="INT",
        default=28,
        help="Phred score to be used as threshold for cleaning (filtering) fastq files.",
    )
    parser.add_argument(
        "-w",
        "--window-size",
        type=int,
        metavar="INT",
        default=5,
        help="Window size to use for cleaning (filtering) fastq files.",
    )
    parser.add_argument(
        "-ml",
        "--minimum-length",
        type=int,
        metavar="INT",
        default=50,
        help="Minimum length for fastq reads to be kept after trimming.",
    )
    parser.add_argument(
        "-k",
        "--kmer-size",
        nargs="+",
        type=int,
        metavar="INT INT...",
        default=[21, 33, 55, 77, 99],
        help="Kmersizes to be used for the de novo assembly.",
    )
    parser.add_argument(
        "-cl",
        "--contig-length-threshold",
        type=str,
        metavar="INT",
        default=500,
        help="Minimum length to filter the contigs generated by the de novo assembly.",
    )
    parser.add_argument(
        "--no-containers",
        action="store_false",
        help="Use conda environments instead of containers.",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=Path,
        metavar="PATH",
        default=None,
        help="Conda or singularity prefix. Basically a path to the place where you want to store the conda environments or the singularity images.",
    )
    parser.add_argument(
        "-c",
        "--cores",
        type=int,
        metavar="INT",
        default=300 if not "--local" in sys.argv else 4,
        help="Number of cores to use. Default is 4 if running locally (--local) or 300 otherwise.",
    )
    parser.add_argument(
        "-q",
        "--queue",
        type=str,
        metavar="STR",
        default="bio",
        help="Name of the queue that the job will be submitted to if working on a cluster.",
    )
    parser.add_argument(
        "-l",
        "--local",
        action="store_true",
        help="If this flag is present, the pipeline will be run locally (not attempting to send the jobs to an HPC cluster**). The default is to assume that you are working on a cluster. **Note that currently only LSF clusters are supported.",
    )
    parser.add_argument(
        "-tl",
        "--time-limit",
        type=int,
        metavar="INT",
        default=60,
        help="Time limit per job in minutes (passed as -W argument to bsub). Jobs will be killed if not finished in this time.",
    )
    parser.add_argument(
        "-u",
        "--unlock",
        action="store_true",
        help="Unlock output directory (passed to snakemake).",
    )
    parser.add_argument(
        "-n",
        "--dryrun",
        action="store_true",
        help="Dry run printing steps to be taken in the pipeline without actually running it (passed to snakemake).",
    )
    parser.add_argument(
        "--rerunincomplete",
        action="store_true",
        help="Re-run jobs if they are marked as incomplete (passed to snakemake).",
    )
    parser.add_argument(
        "--snakemake-args",
        nargs="*",
        default={},
        action=helper_functions.SnakemakeKwargsAction,
        help="Extra arguments to be passed to snakemake API (https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html).",
    )
    args = parser.parse_args()
    juno_assembly_run = JunoAssemblyRun(
        input_dir=args.input,
        output_dir=args.output,
        genus=args.genus,
        exclusion_file=args.exclusion_file,
        db_dir=args.db_dir,
        help_genera=args.help_genera,
        metadata_file=args.metadata,
        cores=args.cores,
        time_limit=args.time_limit,
        local=args.local,
        queue=args.queue,
        unlock=args.unlock,
        rerunincomplete=args.rerunincomplete,
        dryrun=args.dryrun,
        run_in_container=args.no_containers,
        mean_quality_threshold=args.mean_quality_threshold,
        window_size=args.window_size,
        min_read_length=args.minimum_length,
        kmer_size=args.kmer_size,
        contig_length_threshold=args.contig_length_threshold,
        singularity_prefix=args.prefix,
        restarttimes=2,
        latency_wait=60,
        **args.snakemake_args,
    )
    juno_assembly_run.setup_and_validate()
    juno_assembly_run.run_juno_assembly_pipeline()


if __name__ == "__main__":
    main()
