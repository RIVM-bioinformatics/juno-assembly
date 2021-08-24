"""
Juno_assembly pipeline
Authors: Alejandra Hernandez-Segura
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium
            Surveillance (IDS), Bacteriologie (BPD)     
Date: 18-08-2021   

Documentation: https://rivm-bioinformatics.github.io/ids_bacteriology_man/juno-assembly.html 
"""

import atexit
import base_juno_pipeline
import argparse
import os
import pandas as pd
import pathlib
import sys
import warnings
import yaml

class JunoAssemblyRun(base_juno_pipeline.helper_functions.JunoHelpers):
    """Class with the arguments and specifications that are only for the 
    Juno_assembly pipeline but inherit from PipelineStartup and RunSnakemake
    """

    def __init__(self, 
                input_dir, 
                output_dir, 
                genus=None,
                metadata=None,
                help_genera=False,
                cores=300,
                local=False,
                queue='bio',
                unlock=False,
                rerunincomplete=False,
                dryrun=False,
                update_dbs=False):
        """Initiating Juno_assembly pipeline"""
        
        if help_genera:
            print("The accepted genera are:")
            os.system("cat files/accepted_genera_checkm.txt")
            sys.exit(0)

        # Pipeline attributes
        self.pipeline_info = {'pipeline_name': "Juno_assembly",
                                'pipeline_version': "2.0"}
        self.snakefile = "Snakefile"
        self.sample_sheet = "config/sample_sheet.yaml"
        self.workdir = pathlib.Path(__file__).parent.absolute()
        self.useconda = False
        self.usesingularity = True
        self.user_parameters = pathlib.Path("config/user_parameters.yaml")
        self.restarttimes = 0 
        self.latency_wait = 90
        self.supported_genera=[]
        with open(self.workdir.joinpath('files', 'accepted_genera_checkm.txt')) as file_:
            for line in file_:
                genus_name = line.replace('\n', '').lower()
                self.supported_genera.append(genus_name)

        self.input_dir = pathlib.Path(input_dir).resolve()
        self.output_dir = pathlib.Path(output_dir).resolve()
        self.singularityargs = f"--bind {self.input_dir}:{self.input_dir} --bind {self.output_dir}:{self.output_dir}"
        
        if genus is not None:
            self.genus = genus.strip().lower()
            self.__check_genus_is_supported(self.genus)
        else:
            self.genus = None
        if metadata is not None:
            self.metadata = pathlib.Path(metadata)
        else:
            self.metadata = None
        # Start pipeline
        self.startup = self.start_pipeline()
        self.user_params = self.write_userparameters()
        snakemake_run = base_juno_pipeline.RunSnakemake(pipeline_name = self.pipeline_info['pipeline_name'],
                                            pipeline_version = self.pipeline_info['pipeline_version'],
                                            sample_sheet = self.sample_sheet,
                                            output_dir = self.output_dir,
                                            workdir = self.workdir,
                                            snakefile = self.snakefile,
                                            cores = cores,
                                            local = local,
                                            queue = queue,
                                            unlock = unlock,
                                            rerunincomplete = rerunincomplete,
                                            dryrun = dryrun,
                                            useconda = self.useconda,
                                            usesingularity = self.usesingularity,
                                            singularityargs = self.singularityargs,
                                            restarttimes = self.restarttimes,
                                            latency_wait = self.latency_wait)
        
        self.successful_run = snakemake_run.run_snakemake()
        assert self.successful_run, f'Please check the log files'

    def __check_genus_is_supported(self, genus):
        if genus.lower() in self.supported_genera:
            return True
        else:
            raise ValueError(
                self.error_formatter(
                    f'The genus {genus} is not supported. You can leave the "genus" empty for samples with unsupported genera.'
                    )
                )

    def add_metadata(self, samples_dic):
        assert self.metadata.is_file(), f"Provided metadata file ({self.metadata}) does not exist"
        # Load species file
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            species_dic = pd.read_csv(self.metadata, usecols=['Sample', 'Genus'], index_col=0)
        species_dic.index = species_dic.index.map(str)
        species_dic['Genus'] = species_dic['Genus'].apply(lambda x: x.strip().lower())
        species_dic = species_dic.transpose().to_dict()
        # Update dictionary with species
        for sample_name in samples_dic :
            try:
                self.__check_genus_is_supported(species_dic[sample_name]['Genus'])
                samples_dic[sample_name]['genus'] =  species_dic[sample_name]['Genus']
            except KeyError :
                pass

    def start_pipeline(self):
        """Function to start the pipeline (some steps from PipelineStartup 
        need to be modified for the Juno_assembly pipeline to accept metadata
        """
        
        startup = base_juno_pipeline.PipelineStartup(self.input_dir, input_type = 'fastq')
        
        # Add genus metadata if existing
        for sample in startup.sample_dict:
            startup.sample_dict[sample]['genus'] = self.genus
        if self.metadata is not None:
            print('\nAdding genus information from metadata file...\n')
            self.add_metadata(startup.sample_dict)
        # Write sample_sheet
        with open(self.sample_sheet, 'w') as file:
            yaml.dump(startup.sample_dict, file, default_flow_style=False)
        return startup
    
    def write_userparameters(self):

        config_params = {'input_dir': str(self.input_dir),
                        'out': str(self.output_dir),
                        'genus': self.genus,
                        'using_containers': self.usesingularity}
        
        with open(self.user_parameters, 'w') as file:
            yaml.dump(config_params, file, default_flow_style=False)

        return config_params
    
    

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
        help = "Relative or absolute path to the metadata .csv file. If provided, it must contain at least one column with the 'Sample' name (name of the file but removing _R1.fastq.gz) and a column called 'Genus' (mind the capital in the first letter). The genus provided will be used to choose the reference genome to analyze de QC of the de novo assembly."
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
        "-c",
        "--cores",
        type = int,
        metavar = "INT",
        default = 300,
        help="Number of cores to use. Default is 300"
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
        help="Running pipeline locally (instead of in a computer cluster). Default is running it in a cluster."
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
    args = parser.parse_args()
    juno_assembly_run = JunoAssemblyRun(input_dir=args.input, 
                                            genus=args.genus,
                                            output_dir=args.output, 
                                            help_genera=args.help_genera,
                                            metadata=args.metadata,
                                            cores=args.cores,
                                            local=args.local,
                                            queue=args.queue,
                                            unlock=args.unlock,
                                            rerunincomplete=args.rerunincomplete,
                                            dryrun=args.dryrun)