PIPELINE_NAME="$1"
VERSION="$2"

echo """
${PIPELINE_NAME}
Version used: 'https://github.com/AleSR13/Juno-typing/tree/${VERSION}'
Built with Snakemake
  Usage: 
bash $(echo "${PIPELINE_NAME}" | awk '{print tolower($0)}') -i <INPUT_DIR> --genus <genus> <parameters>

  N.B. it is designed for Illumina paired-end data only

Input:
  -i, --input [DIR]                 This is the folder containing your input 
                                    fastq files. Default 
                                    is raw_data/
  -o, --output [DIR]                This is the folder containing your output 
                                    files. Default is output/
                                    
  --genus [STR]                     Provide name of an approved genus to be used
                                    for CheckM. Please check --help-genera for a 
                                    list of approved genera). Only one name is 
                                    allowed, so if multiple samples are included 
                                    in the input dir, it will be assumed that 
                                    they all have the same genus. If different 
                                    samples have different genus, then do not 
                                    provide a --genus and use --metadata instead

  --metadata  [FILE]                Excel file (.xlsx) containing the information 
                                    of the samples and their corresponding genus. 
                                    It should contain at least a column called 
                                    'Monsternummer' containing the sample ID (same 
                                    name than fastq files but removing the suffix 
                                    _R1|2.fastq.gz) and another one called 'genus'
                                    containing the name of the genus. Mind the 
                                    capital in the name 'Monsternummer'. Default
                                    is: /data/BioGrid/NGSlab/BAC_in_house_NGS/In-
                                    house_NGS_selectie_2021.xlsx

  --queue, -q [STR]                 If using a cluster, this is the name of the queue
                                    to which the jobs should be sent. Default
                                    is 'bio'.  
  --cores [INT]                     Number of cores to use to run the pipeline. If
                                    running in a cluster, the default is 300. If 
                                    running locally, the default is 4.
                                  
  --local, -l                       If this flag is present, the pipeline is run 
                                    locally instead of in a cluster. The default 
                                    is to run in a cluster ('bio' queue)

Output:
  out/                              Contains dir contains the results of every 
                                    step of the pipeline.

  out/log/                          Contains the log files for every step of the 
                                    pipeline

  out/log/drmaa                     Contains the .out and .err files of every job 
                                    sent to the grid/cluster.

  out/log/results                  Contains the log files and parameters that 
                                    the pipeline used for the current run


Parameters:
  -h, --help                        Print the help document.

  --help-genera                     Prints list of accepted genera for this 
                                    pipeline (based on CheckM taxon list).

  -sh, --snakemake-help             Print the snakemake help document.

  --clean (-y)                      Removes output (-y forces 'Yes' on all prompts).
  
  --no-checkm                       Not run CheckM or update the genus database 
                                    from CheckM

  --no-genus-update                 Not update the genus database from CheckM

  -n, --dry-run                     Useful snakemake command that displays 
                                    the steps to be performed without actually 
                                    executing them. Useful to spot any potential 
                                    issues while running the pipeline.

  -u, --unlock                      Unlocks the working directory. A directory 
                                    is locked when a run ends abruptly and it 
                                    prevents you from doing subsequent analyses 
                                    on that directory until it gets unlocked.

  Other snakemake parameters        Any other parameters will be passed to 
                                    snakemake. Read snakemake help (-sh) to see
                                    the options.
                                    
"""