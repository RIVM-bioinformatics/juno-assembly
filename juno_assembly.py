"""
Juno_assembly pipeline
Authors: Alejandra Hernandez-Segura
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium
            Surveillance (IDS), Bacteriologie (BPD)     
Date: 18-08-2021   

Documentation: https://rivm-bioinformatics.github.io/ids_bacteriology_man/juno-assembly.html 
"""

from base_juno_pipeline import *
import argparse
import os
import pathlib
import sys
import yaml

class JunoAssemblyRun(base_juno_pipeline.PipelineStartup,
                        base_juno_pipeline.RunSnakemake):
    """Class with the arguments and specifications that are only for the 
    Juno_assembly pipeline but inherit from PipelineStartup and RunSnakemake
    """

    def __init__(self, 
                input_dir, 
                output_dir, 
                db_dir='/mnt/db/juno/kraken2_db',
                genus=None,
                metadata=None,
                mean_quality_threshold=28,
                window_size=5,
                min_read_length=50,
                kmer_size=[21,33,55,77,99],
                contig_length_threshold=500,
                help_genera=False,
                cores=300,
                local=False,
                queue='bio',
                unlock=False,
                rerunincomplete=False,
                dryrun=False,
                run_in_container=False,
                prefix=None,
                **kwargs):
        """Initiating Juno_assembly pipeline"""
        
        if help_genera:
            print("The accepted genera are:")
            os.system("cat files/accepted_genera_checkm.txt")
            sys.exit(0)

        # Build class
        output_dir = pathlib.Path(output_dir).resolve()
        workdir = pathlib.Path(__file__).parent.resolve()
        self.db_dir = pathlib.Path(db_dir).resolve()
        self.path_to_audit = output_dir.joinpath('audit_trail')
        base_juno_pipeline.PipelineStartup.__init__(self,
            input_dir=pathlib.Path(input_dir).resolve(), 
            input_type='fastq',
            min_num_lines=1000) # TODO: Find ideal min num of reads/lines needed
        base_juno_pipeline.RunSnakemake.__init__(self,
            pipeline_name='Juno_assembly',
            pipeline_version='v2.0.2',
            output_dir=output_dir,
            workdir=workdir,
            cores=cores,
            local=local,
            queue=queue,
            unlock=unlock,
            rerunincomplete=rerunincomplete,
            dryrun=dryrun,
            useconda=not run_in_container,
            conda_prefix=prefix,
            usesingularity=run_in_container,
            singularityargs=f"--bind {self.input_dir}:{self.input_dir} --bind {output_dir}:{output_dir} --bind {db_dir}:{db_dir}",
            singularity_prefix=prefix,
            restarttimes=2,
            latency_wait=60,
            name_snakemake_report=str(self.path_to_audit.joinpath('juno_assembly_report.html')),
            **kwargs)

        # Specific for Juno assembly
        self.mean_quality_threshold = int(mean_quality_threshold)
        self.window_size = int(window_size)
        self.min_read_length = int(min_read_length)
        self.kmer_size = ','.join([str(ks) for ks in kmer_size])
        self.contig_length_threshold = int(contig_length_threshold)
        self.supported_genera=[]
        with open(self.workdir.joinpath('files', 'accepted_genera_checkm.txt')) as file_:
            for line in file_:
                genus_name = line.replace('\n', '').lower()
                self.supported_genera.append(genus_name)
        self.genus = genus
        self.metadata_file = metadata
        
        # Start pipeline
        self.run_juno_assembly_pipeline()

    def __get_supported_genera(self):
        supported_genera=[]
        with open(self.workdir.joinpath('files', 'accepted_genera_checkm.txt')) as file_:
            for line in file_:
                genus_name = line.replace('\n', '').lower()
                supported_genera.append(genus_name)
        self.supported_genera = supported_genera

    def __check_genus_is_supported(self):
        self.__get_supported_genera()
        if self.genus is not None:
            if self.genus.strip().lower() in self.supported_genera:
                return True
            else:
                raise ValueError(
                    self.error_formatter(
                        f'The genus {self.genus} is not supported. You can leave the "genus" empty for samples with unsupported genera.'
                        )
                    )
        else:
            return True

    def __validate_kraken2_db_dir(self):
        hash_file_exists = self.db_dir.joinpath('hash.k2d').exists()
        opts_file_exists =  self.db_dir.joinpath('opts.k2d').exists()
        taxo_file_exists = self.db_dir.joinpath('taxo.k2d').exists()
        if hash_file_exists and opts_file_exists and taxo_file_exists:
            return True
        else:
            raise ValueError(f'The provided path to the database for Kraken2 ({str(self.db_dir)}) does not contain the expected files. Please download it again!')

    def update_sample_dict_with_metadata(self):
        self.__check_genus_is_supported()
        self.get_metadata_from_csv_file(filepath=self.metadata_file, expected_colnames=['sample', 'genus'])
        for sample in self.sample_dict:
            try:
                self.sample_dict[sample]['genus'] = self.juno_metadata[sample]['genus'].strip().lower()
            except (KeyError, TypeError):
                self.sample_dict[sample]['genus'] = self.genus

    def start_juno_assembly_pipeline(self):
        """
        Function to start the pipeline (some steps from PipelineStartup need to
        be modified for the Juno_assembly pipeline to accept metadata
        """
        self.start_juno_pipeline()
        self.update_sample_dict_with_metadata()
        with open(self.sample_sheet, 'w') as file:
            yaml.dump(self.sample_dict, file, default_flow_style=False)
    
    def write_userparameters(self):

        config_params = {'input_dir': str(self.input_dir),
                        'out': str(self.output_dir),
                        'genus': self.genus,
                        'mean_quality_threshold': self.mean_quality_threshold,
                        'window_size': self.window_size,
                        'min_read_length': self.min_read_length,
                        'kmer_size': self.kmer_size, 
                        'contig_length_threshold': self.contig_length_threshold,
                        'run_in_container': self.usesingularity,
                        'db_dir': str(self.db_dir)}
        
        with open(self.user_parameters, 'w') as file:
            yaml.dump(config_params, file, default_flow_style=False)

        return config_params

    def run_juno_assembly_pipeline(self):
        self.start_juno_assembly_pipeline()
        self.user_params = self.write_userparameters()
        self.get_run_info()
        if not self.dryrun or self.unlock:
            self.__validate_kraken2_db_dir()
        self.successful_run = self.run_snakemake()
        assert self.successful_run, f'Please check the log files'
        if not self.dryrun or self.unlock:
            self.make_snakemake_report()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "Juno_assembly pipeline. Automated pipeline for pre-processing, QC and assembly of bacterial NGS sequencing data."
    )
    parser.add_argument(
        "--help-genera",
        action = 'store_true',
        help = "Prints the genera accepted by this pipeline."
    )
    parser.add_argument(
        "-i",
        "--input",
        type = pathlib.Path,
        required = not '--help-genera' in sys.argv,
        metavar = "DIR",
        help = "Relative or absolute path to the input directory. It must contain all the raw reads (fastq) files for all samples to be processed (not in subfolders)."
    )
    parser.add_argument(
        "-g",
        "--genus",
        type = str.lower,
        default = None,
        metavar = "GENUS",
        help = "Genus of the samples to be analyzed. If metadata is given, the genus in the metadata will overwrite the one given through this option."
    )
    parser.add_argument(
        "-m",
        "--metadata",
        type = pathlib.Path,
        default = None,
        metavar = "FILE",
        help = "Relative or absolute path to a .csv file. If provided, it must contain at least one column with the 'Sample' name (name of the file but removing _R1.fastq.gz) and a column called 'Genus' (mind the capital in the first letter). The genus provided will be used to choose the reference genome to analyze de QC of the de novo assembly."
    )
    parser.add_argument(
        "-o",
        "--output",
        type = pathlib.Path,
        metavar = "DIR",
        default = "output",
        help = "Relative or absolute path to the output directory. If non is given, an 'output' directory will be created in the current directory."
    )
    parser.add_argument(
        "-d",
        "--db-dir",
        type = pathlib.Path,
        metavar = "DIR",
        default = "/mnt/db/juno/kraken2_db",
        help = "Relative or absolute path to the Kraken2 database. Default: /mnt/db/juno/kraken2_db."
    )
    parser.add_argument(
        "-mpt",
        "--mean-quality-threshold",
        type = int,
        metavar = "INT",
        default = 28,
        help = "Phred score to be used as threshold for cleaning (filtering) fastq files."
    )
    parser.add_argument(
        "-w",
        "--window-size",
        type = int,
        metavar = "INT",
        default = 5,
        help = "Window size to use for cleaning (filtering) fastq files."
    )
    parser.add_argument(
        "-ml",
        "--minimum-length",
        type = int,
        metavar = "INT",
        default = 50,
        help = "Minimum length for fastq reads to be kept after trimming."
    )
    parser.add_argument(
        "-k",
        "--kmer-size",
        type = str,
        nargs='+',
        metavar = "INT INT...",
        default = [21,33,55,77,99],
        help = "Kmersizes to be used for the de novo assembly."
    )
    parser.add_argument(
        "-cl",
        "--contig-length-threshold",
        type = str,
        metavar = "INT",
        default = 500,
        help = "Minimum length to filter the contigs generated by the de novo assembly."
    )
    parser.add_argument(
        "--no-containers",
        action = 'store_false',
        help = "Use conda environments instead of containers."
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type = pathlib.Path,
        metavar="PATH",
        default=None,
        help = "Conda or singularity prefix. Basically a path to the place where you want to store the conda environments or the singularity images."
    )
    parser.add_argument(
        "-c",
        "--cores",
        type = int,
        metavar = "INT",
        default = 300 if not '--local' in sys.argv else 4,
        help="Number of cores to use. Default is 4 if running locally (--local) or 300 otherwise."
    )
    parser.add_argument(
        "-q",
        "--queue",
        type = str,
        metavar = "STR",
        default = 'bio',
        help = 'Name of the queue that the job will be submitted to if working on a cluster.'
    )
    parser.add_argument(
        "-l",
        "--local",
        action='store_true',
        help="If this flag is present, the pipeline will be run locally (not attempting to send the jobs to an HPC cluster**). The default is to assume that you are working on a cluster. **Note that currently only LSF clusters are supported."
    )
    parser.add_argument(
        "-u",
        "--unlock",
        action = 'store_true',
        help = "Unlock output directory (passed to snakemake)."
    )
    parser.add_argument(
        "-n",
        "--dryrun",
        action='store_true',
        help="Dry run printing steps to be taken in the pipeline without actually running it (passed to snakemake)."
    )
    parser.add_argument(
        "--rerunincomplete",
        action='store_true',
        help="Re-run jobs if they are marked as incomplete (passed to snakemake)."
    )
    parser.add_argument(
        "--snakemake-args",
        nargs='*',
        default={},
        action=helper_functions.SnakemakeKwargsAction,
        help="Extra arguments to be passed to snakemake API (https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html)."
    )
    args = parser.parse_args()
    juno_assembly_run = JunoAssemblyRun(input_dir=args.input, 
                                            genus=args.genus,
                                            output_dir=args.output, 
                                            db_dir=args.db_dir,
                                            help_genera=args.help_genera,
                                            metadata=args.metadata,
                                            cores=args.cores,
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
                                            prefix=args.prefix,
                                            **args.snakemake_args)