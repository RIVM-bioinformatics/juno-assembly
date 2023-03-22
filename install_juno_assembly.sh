# Make sure to run from the main directory of the pipeline
set -euo pipefail
path_to_main=$(realpath $0)
workdir=$(dirname "$path_to_main")
cd "$workdir"

# Delete previous installations
rm -rf envs/src/
mamba env update -f envs/juno_assembly.yaml
