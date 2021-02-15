#!/bin/bash
set -x

INPUT_DIR=$1
OUTPUT_DIR=$2

if [ -z "$3" ]; then
    GENUS_ALL="NotProvided"
else
    GENUS_ALL=$3
fi

if [ $GENUS_ALL == "NotProvided" ]; then
    CHECKM="FALSE"
else 
    CHECKM="TRUE"
fi

if [ ! -z "$4" ]; then
    METADATA=$4
    CHECKM="TRUE"
fi

conda env update -f envs/mamba.yaml -q -v
source activate mamba

PATH_MASTER_YAML="envs/master_env.yaml"
MASTER_NAME="$(head -n 1 ${PATH_MASTER_YAML} | cut -f2 -d ' ')" # Extract Conda environment name as specified in yaml file

mamba env update -f $PATH_MASTER_YAML -q -v

# if [ $CHECKM == "TRUE" ]; then
#     mamba env update -f envs/checkm.yaml -q -v
#     source activate checkM
#     checkm taxon_list > files/checkm_taxon_list.txt
#     source activate mamba # back to mamba again to start juno_master
# fi

source activate $MASTER_NAME

python bin/generate_sample_sheet.py "${INPUT_DIR}" > sample_sheet.yaml

UNIQUE_ID=$(bin/include/generate_id.sh)
SET_HOSTNAME=$(bin/include/gethostname.sh)

echo -e "pipeline_run:\n    identifier: ${UNIQUE_ID}" > profile/variables.yaml
echo -e "Server_host:\n    hostname: http://${SET_HOSTNAME}" >> profile/variables.yaml

if [ -z $METADATA ]; then
    echo -e "Juno call (this settings would overwrite any others in the configuration files): \n" > profile/juno_call.txt
    echo -e "snakemake --config checkm=$CHECKM out=$OUTPUT_DIR genus=$GENUS_ALL \
--profile profile \
--drmaa ' -q ${irods_runsheet_sys__runsheet__lsf_queue} -n {threads} -R \'span[hosts=1]\'' \
--drmaa-log-dir ${OUTPUT_DIR}/log/drmaa \
${@}\n" >> profile/juno_call.txt
    snakemake --config checkm=$CHECKM out=$OUTPUT_DIR genus=$GENUS_ALL --profile profile --drmaa " -q ${irods_runsheet_sys__runsheet__lsf_queue} -n {threads} -R \"span[hosts=1]\"" --drmaa-log-dir ${OUTPUT_DIR}/log/drmaa
else
    echo "This is the genus file: $METADATA"
    echo -e "Juno call (this settings would overwrite any others in the configuration files): \n" > profile/juno_call.txt
    echo -e "snakemake --config checkm=$CHECKM out=$OUTPUT_DIR genus=$GENUS_ALL metadata=$METADATA --profile profile \
--drmaa  -q ${irods_runsheet_sys__runsheet__lsf_queue} -n {threads} -R \'span[hosts=1]\'' \
--drmaa-log-dir ${OUTPUT_DIR}/log/drmaa ${@} \n" >> profile/juno_call.txt
    snakemake --config checkm=$CHECKM out=$OUTPUT_DIR genus=$GENUS_ALL metadata=$METADATA --profile profile --drmaa " -q ${irods_runsheet_sys__runsheet__lsf_queue} -n {threads} -R \"span[hosts=1]\"" --drmaa-log-dir ${OUTPUT_DIR}/log/drmaa
fi 
