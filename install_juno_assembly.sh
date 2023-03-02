# Make sure to run from the main directory of the pipeline
set -euo pipefail
path_to_main=$(realpath $0)
workdir=$(dirname "$path_to_main")
cd "$workdir"

# Install conda
set +euo pipefail
conda env update -f envs/mamba.yaml
source activate mamba
set -euo pipefail

# Delete previous installations
rm -rf envs/src/
mamba env update -f envs/juno_assembly.yaml
