#!/bin/bash
set -x

INPUT_DIR=$1
OUTPUT_DIR=$2
CHECKM="FALSE"

conda env update -f environments/mamba.yaml -q -v
source activate mamba

mamba env update -f environments/master_env.yaml -q -v
source activate juno_master

python bin/generate_sample_sheet.py "${INPUT_DIR}" > sample_sheet.yaml

UNIQUE_ID=$(bin/generate_id.sh)
SET_HOSTNAME=$(bin/gethostname.sh)

echo -e "pipeline_run:\n    identifier: ${UNIQUE_ID}" > profile/variables.yaml
echo -e "Server_host:\n    hostname: http://${SET_HOSTNAME}" >> profile/variables.yaml
#eval $(parse_yaml profile/variables.yaml "config_")

snakemake -s Snakefile --config checkm=$CHECKM out=$OUTPUT_DIR --profile profile --drmaa " -q bio -n {threads} -R \"span[hosts=1]\"" --drmaa-log-dir ${OUTPUT_DIR}/log/drmaa
