#!/bin/bash
set -x
# use set -x for extra logging

# Wrapper for juno assembly pipeline

#----------------------------------------------#
# User parameters
input_dir="${1%/}"
output_dir="${2%/}"
PROJECT_NAME="${irods_input_projectID}" # This should be an environment variable

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" > /dev/null 2>&1 && pwd )"
cd ${DIR}

# Sanity checks
if [ ! -z "${1}" ] || [ ! -z "${2}" ] || [ ! -z "${irods_input_projectID}" ]
then
  INPUTDIR="${1}"
  OUTPUTDIR="${2}"
  PROJECT_NAME="${irods_input_projectID}"
else
  echo "No inputdir, outputdir or project name (param 1, 2, irods_input_projectID)"
  exit 1
fi

#check if there is an exclusion file, if so change the parameter
if [ ! -z "${irods_input_sequencing__run_id}" ] && [ -f "/data/BioGrid/NGSlab/sample_sheets/${irods_input_sequencing__run_id}.exclude" ]
then
  EXCLUSION_FILE_COMMAND="-ex /data/BioGrid/NGSlab/sample_sheets/${irods_input_sequencing__run_id}.exclude"
else
  EXCLUSION_FILE_COMMAND=""
fi

if [ ! -d "${INPUTDIR}" ] || [ ! -d "${OUTPUTDIR}" ]
then
  echo "inputdir $INPUTDIR or output dir $OUTPUTDIR does not exist!"
  exit 1
fi

set -euo pipefail

#----------------------------------------------#
## make sure conda works

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/mnt/miniconda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/mnt/miniconda/etc/profile.d/conda.sh" ]; then
        . "/mnt/miniconda/etc/profile.d/conda.sh"
    else
        export PATH="/mnt/miniconda/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<export -f conda
export -f __conda_activate
export -f __conda_reactivate
export -f __conda_hashr


#----------------------------------------------#
# we can use the base installation of mamba to create the environment. 
# Swapping to a parent env is not necessary anymore.
mamba env create -f envs/juno_assembly.yaml --name pipeline_env
conda activate pipeline_env


#----------------------------------------------#
# Run the pipeline

case $PROJECT_NAME in
  adhoc|rogas|svgasuit|bsr_rvp)
    GENUS_COMMAND=""
    ;;
  dsshig|svshig)
    GENUS_COMMAND="--genus Shigella"
    ;;
  salm|svsalent|svsaltyp|vdl_salm)
    GENUS_COMMAND="--genus Salmonella"
    ;;
  svlismon|vdl_list)
    GENUS_COMMAND="--genus Listeria"
    ;;
  svstec|vdl_ecoli|vdl_stec)
    GENUS_COMMAND="--genus Escherichia"
    ;;
  campy|vdl_campy)
    GENUS_COMMAND="--genus Campylobacter"
    ;;
  refsamp)
    GENUS_FILE=`realpath $(find ../ -type f -name genus_sheet_refsamp.csv)`
    GENUS_COMMAND="--metadata $GENUS_FILE"
    ;;
  *)
    GENUS_COMMAND=""
    ;;
esac

echo -e "\nRun pipeline..."

if [ ! -z ${irods_runsheet_sys__runsheet__lsf_queue} ]; then
    QUEUE="${irods_runsheet_sys__runsheet__lsf_queue}"
else
    QUEUE="bio"
fi

set -euo pipefail
set -x

# Setting up the tmpdir for singularity as the current directory (default is /tmp but it gets full easily)
# Containers will use it for storing tmp files when building a container
export SINGULARITY_TMPDIR="$(pwd)"

python juno_assembly.py \
  --queue $QUEUE \
  -i $input_dir \
  -o $output_dir \
  $EXCLUSION_FILE_COMMAND \
  $GENUS_COMMAND \
  --prefix /mnt/db/juno/sing_containers

result=$?

# Propagate metadata

set +euo pipefail

SEQ_KEYS=
SEQ_ENV=`env | grep irods_input_sequencing`
for SEQ_AVU in ${SEQ_ENV}
do
    SEQ_KEYS="${SEQ_KEYS} ${SEQ_AVU%%=*}"
done

for key in $SEQ_KEYS irods_input_illumina__Flowcell irods_input_illumina__Instrument \
    irods_input_illumina__Date irods_input_illumina__Run_number irods_input_illumina__Run_Id
do
    if [ ! -z ${!key} ] ; then
        attrname=${key:12}
        attrname=${attrname/__/::}
        echo "${attrname}: '${!key}'" >> ${OUTPUTDIR}/metadata.yml
    fi
done

exit ${result}
