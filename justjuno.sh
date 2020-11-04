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
    GENUS_FILE=$4
    CHECKM="TRUE"
fi

conda env update -f environments/mamba.yaml -q -v
source activate mamba

mamba env update -f environments/master_env.yaml -q -v

if [ $CHECKM == "TRUE" ]; then
    mamba env update -f environments/CheckM.yaml -q -v
    source activate checkM
    checkm taxon_list > checkm_taxon_list.txt
fi

source activate juno_master

python bin/generate_sample_sheet.py "${INPUT_DIR}" > sample_sheet.yaml

UNIQUE_ID=$(bin/generate_id.sh)
SET_HOSTNAME=$(bin/gethostname.sh)

echo -e "pipeline_run:\n    identifier: ${UNIQUE_ID}" > profile/variables.yaml
echo -e "Server_host:\n    hostname: http://${SET_HOSTNAME}" >> profile/variables.yaml
#eval $(parse_yaml profile/variables.yaml "config_")

if [ -z $GENUS_FILE ]; then
    snakemake --config checkm=$CHECKM out=$OUTPUT_DIR genus=$GENUS_ALL --profile profile --drmaa " -q bio -n {threads} -R \"span[hosts=1]\"" --drmaa-log-dir ${OUTPUT_DIR}/log/drmaa
else
    snakemake --config checkm=$CHECKM out=$OUTPUT_DIR genus=$GENUS_ALL genus_file=$GENUS_FILE --profile profile --drmaa " -q bio -n {threads} -R \"span[hosts=1]\"" --drmaa-log-dir ${OUTPUT_DIR}/log/drmaa
fi