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
else
    METADATA="/data/BioGrid/NGSlab/BAC_in_house_NGS/In-house_NGS_selectie_2021.xlsx"
fi

conda env update -f envs/mamba.yaml -q -v
source activate mamba

PATH_MASTER_YAML="envs/master_env.yaml"
MASTER_NAME="$(head -n 1 ${PATH_MASTER_YAML} | cut -f2 -d ' ')" 

mamba env update -f $PATH_MASTER_YAML -q -v

source activate $MASTER_NAME

python bin/generate_sample_sheet.py "${INPUT_DIR}" > sample_sheet.yaml

UNIQUE_ID=$(bin/include/generate_id.sh)
SET_HOSTNAME=$(bin/include/gethostname.sh)
VERSION=$(git log -n 1 --pretty=format:"%H")

echo -e "pipeline_run:\n    identifier: ${UNIQUE_ID}" > profile/variables.yaml
echo -e "Server_host:\n    hostname: http://${SET_HOSTNAME}" >> profile/variables.yaml

__USERPARAMETERS="
pipeline_version: "${VERSION}"
out: "${OUTPUT_DIR}"
metadata: "$METADATA"
checkm: "$CHECKM"
genus: "$GENUS_ALL"
cores: 300
"

echo "$__USERPARAMETERS" > profile/user_parameters.yaml

snakemake --profile profile --cores 300 \
    --drmaa " -n {threads} \
    -o ${OUTPUT_DIR}/log/drmaa/{name}_{wildcards}_{jobid}.out \
    -e ${OUTPUT_DIR}/log/drmaa/{name}_{wildcards}_{jobid}.err \
    -R \"span[hosts=1] rusage[mem={resources.mem_mb}]\" "  \
    --drmaa-log-dir ${OUTPUT_DIR}/log/drmaa
