#!/bin/bash

# Wrapper for juno assembly pipeline

set -euo pipefail

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
  echo "No inputdir, outputdir or project name (param 1, 2 or irods_input_projectID)"
  exit 1
fi

if [ ! -d "${INPUTDIR}" ] || [ ! -d "${OUTPUTDIR}" ]
then
  echo "inputdir $INPUTDIR or output dir $OUTPUTDIR does not exist!"
  exit 1
fi

#----------------------------------------------#
# Create/update necessary environments
PATH_MAMBA_YAML="envs/mamba.yaml"
PATH_MASTER_YAML="envs/master_env.yaml"
MAMBA_NAME=$(head -n 1 ${PATH_MAMBA_YAML} | cut -f2 -d ' ')
MASTER_NAME=$(head -n 1 ${PATH_MASTER_YAML} | cut -f2 -d ' ')

echo -e "\nUpdating necessary environments to run the pipeline..."

# Removing strict mode because it sometimes breaks the code for 
# activating an environment and for testing whether some variables
# are set or not
set +euo pipefail 
bash install_juno_assembly.sh
source activate "${MAMBA_NAME}"
source activate "${MASTER_NAME}"

#----------------------------------------------#
# Run the pipeline

case $PROJECT_NAME in
  adhoc|rogas|svgasuit|bsr_rvp)
    GENUS_ALL="NotProvided"
    ;;
  dsshig|svshig)
    GENUS_ALL="Shigella"
    ;;
  salm|svsalent|svsaltyp|vdl_salm)
    GENUS_ALL="Salmonella"
    ;;
  svlismon|vdl_list)
    GENUS_ALL="Listeria"
    ;;
  svstec|vdl_ecoli|vdl_stec)
    GENUS_ALL="Escherichia"
    ;;
  campy|vdl_campy)
    GENUS_ALL="Campylobacter"
    ;;
  *)
    GENUS_ALL="NotProvided"
    ;;
esac


echo -e "\nRun pipeline..."

if [ ! -z ${irods_runsheet_sys__runsheet__lsf_queue} ]; then
    QUEUE="${irods_runsheet_sys__runsheet__lsf_queue}"
else
    QUEUE="bio"
fi

set -euo pipefail

# Setting up the tmpdir for singularity as the current directory (default is /tmp but it gets full easily)
# Containers will use it for storing tmp files when building a container
export SINGULARITY_TMPDIR="$(pwd)"

if [ "${irods_input_projectID}" == "refsamp" ]; then
    
    GENUS_FILE=`realpath $(find ../ -type f -name genus_sheet_refsamp.csv)`
    
    python juno_assembly.py --queue "${QUEUE}" \
      -i "${input_dir}" \
      -o "${output_dir}" \
      --metadata "${GENUS_FILE}" \
      --prefix "/mnt/db/juno/sing_containers"
    
    result=$?

elif [ "${GENUS_ALL}" == "NotProvided" ]; then

    python juno_assembly.py --queue "${QUEUE}" \
      -i "${input_dir}" \
      -o "${output_dir}" \
      --prefix "/mnt/db/juno/sing_containers"

    result=$?

else

    python juno_assembly.py --queue "${QUEUE}" \
      -i "${input_dir}" \
      -o "${output_dir}" \
      --genus "${GENUS_ALL}" \
      --prefix "/mnt/db/juno/sing_containers"
    
    result=$?

fi

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

# Produce svg with rules
# snakemake --config sample_sheet=config/sample_sheet.yaml \
#             --configfiles config/pipeline_parameters.yaml config/user_parameters.yaml \
#             -j 1 --use-conda \
#             --rulegraph | dot -Tsvg > files/DAG.svg
