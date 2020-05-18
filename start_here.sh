#!/bin/bash

#load in functions
source bin/functions.sh
eval "$(parse_yaml profile/pipeline_parameters.yaml "params_")"
eval "$(parse_yaml profile/config.yaml "configuration_")"

UNIQUE_ID=$(bin/generate_id.sh)
SET_HOSTNAME=$(bin/gethostname.sh)


### conda environment
PATH_MASTER_YAML="environments/master_env.yaml"
MASTER_NAME=$(head -n 1 ${PATH_MASTER_YAML} | cut -f2 -d ' ') # Extract Conda environment name as specified in yaml file
PATH_CHECKM_YAML="environments/CheckM.yaml"
CHECKM_NAME=$(head -n 1 ${PATH_MASTER_YAML} | cut -f2 -d ' ') # Extract Conda environment name as specified in yaml file

### Default values for CLI parameters
INPUT_DIR="raw_data/"
SKIP_CONFIRMATION="FALSE"
SNAKEMAKE_UNLOCK="FALSE"
CLEAN="FALSE"
HELP="FALSE"
MAKE_SAMPLE_SHEET="FALSE"
SHEET_SUCCESS="FALSE"
UPDATE_GENUS="TRUE"

### Parse the commandline arguments, if they are not part of the pipeline, they get send to Snakemake
POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        -i|--input)
        INPUT_DIR="$2"
        shift # Next
        shift # Next
        ;;
        -h|--help)
        HELP="TRUE"
        shift # Next
        ;;
        -sh|--snakemake-help)
        SNAKEMAKE_HELP="TRUE"
        shift # Next
        ;;
        --clean)
        CLEAN="TRUE"
        shift # Next
        ;;
        --make-sample-sheet)
        MAKE_SAMPLE_SHEET="TRUE"
        shift # Next
        ;;
        -y)
        SKIP_CONFIRMATION="TRUE"
        shift # Next
        ;;
        -u|--unlock)
        SNAKEMAKE_UNLOCK="TRUE"
        shift # Next
        ;;
        *) # Any other option
        POSITIONAL+=("$1") # save in array
        shift # Next
        ;;
    esac
done
set -- "${POSITIONAL[@]:-}" # Restores the positional arguments (i.e. without the case arguments above) which then can be called via `$@` or `$[0-9]` etc. These parameters are send to Snakemake.


### Remove all output
if [ "${CLEAN:-}" == "TRUE" ]; then
    line
    spacer
    echo -e "The following files and folders will be deleted:\ndata/\nlogs/\nresults/\nprofile/variables.yaml\nsample_sheet.yaml\n\n"
    if [ "${SKIP_CONFIRMATION}" == "TRUE" ]; then
        echo -e "Removing output: data/ logs/ results/ profile/variables.yaml sample_sheet.yaml"
            rm -rf data/
            rm -rf logs/
            rm -rf out/
            #rm -f sample_sheet.yaml
            rm -f profile/variables.yaml

            ## clean the yaml files
            sed -i '\|databases|d' profile/pipeline_parameters.yaml
            sed -i '\|background_ref|d' profile/pipeline_parameters.yaml
            sed -i '\|Krona_taxonomy|d' profile/pipeline_parameters.yaml
            sed -i '\|virusHostDB|d' profile/pipeline_parameters.yaml
            sed -i '\|NCBI_new_taxdump_rankedlineage|d' profile/pipeline_parameters.yaml
            sed -i '\|NCBI_new_taxdump_host|d' profile/pipeline_parameters.yaml
            #sed -i '\|drmaa|d' profile/config.yaml
    fi
    exit 0
fi


### Print bac_gastro help message
if [ "${HELP:-}" == "TRUE" ]; then
    line
    cat <<HELP_USAGE
Bac_gastro pipeline, version $VERSION, built with Snakemake
  Usage: bash $0 -i <INPUT_DIR> <parameters>
  N.B. it is designed for Illumina paired-end data only

Input:
  -i, --input [DIR]                 This is the folder containing your input fastq files.
                                    Default iSNAKEMAKE_UNLOCKaw_data/' and only relative paths are accepted.
Output (automatically generated):
  out/                             Contains dSNAKEMAKE_UNLOCKled intermediate files.
  logs/                             Contains SNAKEMAKE_UNLOCKlog files.

Parameters:
  -h, --help                        Print theSNAKEMAKE_UNLOCKp document.
  -sh, --snakemake-help             Print theSNAKEMAKE_UNLOCKkemake help document.
  --clean (-y)                      Removes oSNAKEMAKE_UNLOCKt. (-y forces "Yes" on all prompts)
  -n, --dry-run                     Useful snSNAKEMAKE_UNLOCKake command: Do not execute anything, and
                                    display wSNAKEMAKE_UNLOCKwould be done.
  -u, --unlock                      Removes tSNAKEMAKE_UNLOCKock on the working directory. This happens when
                                    a run endSNAKEMAKE_UNLOCKruptly and prevents you from doing subsequent
                                    analyses.SNAKEMAKE_UNLOCK
  -q, --quiet                       Useful snakemake command: Do not output any progress or
                                    rule information.

HELP_USAGE
    exit 0
fi


###############################################################################################################
##### Installation block                                                                                  #####
###############################################################################################################

### Pre-flight check: Assess availability of required files, conda and master environment
if [ ! -e "${PATH_MASTER_YAML}" ]; then # If this yaml file does not exist, give error.
    line
    spacer
    echo -e "ERROR: Missing file \"${PATH_MASTER_YAML}\""
    exit 1
fi

if [[ $PATH != *${MASTER_NAME}* ]]; then # If the master environment is not in your path (i.e. it is not currently active), do...
    line
    spacer
    set +ue # Turn bash strict mode off because that breaks conda
    source activate "${MASTER_NAME}" # Try to activate this env
    if [ ! $? -eq 0 ]; then # If exit statement is not 0, i.e. master conda env hasn't been installed yet, do...
        installer_intro
        if [ "${SKIP_CONFIRMATION}" = "TRUE" ]; then
            echo -e "\tInstalling master environment..." 
            conda env create -f ${PATH_MASTER_YAML} 
            source activate "${MASTER_NAME}"
            echo -e "DONE"
        else
            while read -r -p "The master environment hasn't been installed yet, do you want to install this environment now? [y/N] " envanswer
            do
                envanswer=${envanswer,,}
                if [[ "${envanswer}" =~ ^(yes|y)$ ]]; then
                    echo -e "\tInstalling master environment..." 
                    conda env create -f ${PATH_MASTER_YAML}
                    source activate "${MASTER_NAME}"
                    echo -e "DONE"
                    break
                elif [[ "${envanswer}" =~ ^(no|n)$ ]]; then
                    echo -e "The master environment is a requirement. Exiting because bac_gastro cannot continue without this environment"
                    exit 1
                else
                    echo -e "Please answer with 'yes' or 'no'"
                fi
            done
        fi
    fi
    set -ue # Turn bash strict mode on again
    echo -e "Succesfully activated master environment"
fi


if [ "${SNAKEMAKE_UNLOCK}" == "TRUE" ]; then
    printf "\nUnlocking working directory...\n"
    snakemake -s Snakefile --profile profile --unlock
    printf "\nDone.\n"
    exit 0
fi


### Print Snakemake help
if [ "${SNAKEMAKE_HELP:-}" == "TRUE" ]; then
    line
    snakemake --help
    exit 0
fi


### Pass other CLI arguments along to Snakemake
if [ ! -d "${INPUT_DIR}" ]; then
    minispacer
    echo -e "The input directory specified (${INPUT_DIR}) does not exist"
    echo -e "Please specify an existing input directory"
    exit 1
fi

### Generate sample sheet
if [ -n "$(ls -A "${INPUT_DIR}")" ]; then
    minispacer
    echo -e "Files in input directory (${INPUT_DIR}) are present"
    echo -e "Generating sample sheet..."
    bin/generate_sample_sheet.py "${INPUT_DIR}" > sample_sheet.yaml
    SHEET_SUCCESS="TRUE"
else
    minispacer
    echo -e "The input directory you specified (${INPUT_DIR}) exists but is empty...\nPlease specify a directory with input-data."
    exit 0
fi

### Checker for succesfull creation of sample_sheet
if [ "${SHEET_SUCCESS}" == "TRUE" ]; then
    echo -e "Succesfully generated the sample sheet"
    echo -e "ready_for_start"
else
    echo -e "Couldn't find files in the input directory that ended up being in a .FASTQ, .FQ or .GZ format"
    echo -e "Please inspect the input directory (${INPUT_DIR}) and make sure the files are in one of the formats listed below"
    echo -e ".fastq.gz (Zipped Fastq)"
    echo -e ".fq.gz (Zipped Fq)"
    echo -e ".fastq (Unzipped Fastq)"
    echo -e ".fq (unzipped Fq)"
    exit 1
fi


if [ "${MAKE_SAMPLE_SHEET}" == "TRUE" ]; then
    echo -e "bac_gastro_run:\n    identifier: ${UNIQUE_ID}" > profile/variables.yaml
    echo -e "Server_host:\n    hostname: http://${SET_HOSTNAME}" >> profile/variables.yaml
    echo -e "The sample sheet and variables file has now been created, you can now run the snakefile manually"
    exit 0
fi


if [ "${UPDATE_GENUS}" == "TRUE" ]; then
    printf "\ncollecting available genera from CheckM...\n"
    set +ue # Turn bash strict mode off because that breaks conda
    source activate "${CHECKM_NAME}" # Try to activate checkM env
    if [ ! $? -eq 0 ]; then # If exit statement is not 0, i.e. checkM conda env hasn't been installed yet, do...
        echo -e "\tInstalling master environment..." 
        conda env create -f ${PATH_CHECKM_YAML} 
        source activate "${CHECKM_NAME}"
        echo -e "DONE"
    fi
    checkm taxon_list > checkm_taxon_list.txt  
fi

### Actual snakemake command with checkers for required files. N.B. here the UNIQUE_ID and SET_HOSTNAME variables are set!
if [ -e sample_sheet.yaml ]; then
    echo -e "Starting snakemake"
    set +ue #turn off bash strict mode because snakemake and conda can't work with it properly
    echo -e "pipeline_run:\n    identifier: ${UNIQUE_ID}" > profile/variables.yaml
    echo -e "Server_host:\n    hostname: http://${SET_HOSTNAME}" >> profile/variables.yaml
    eval $(parse_yaml profile/variables.yaml "config_")
    snakemake -s Snakefile --profile profile ${@}
    #echo -e "\nUnique identifier for this run is: $config_run_identifier "
    echo -e "bac gastro run complete"
    set -ue #turn bash strict mode back on
else
    echo -e "Sample_sheet.yaml could not be found"
    echo -e "This also means that the pipeline was unable to generate a new sample sheet for you"
    echo -e "Please inspect the input directory (${INPUT_DIR}) and make sure the right files are present"
    exit 1
fi

exit 0 