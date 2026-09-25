#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --output="logs/%j_Node.log"
#SBATCH --partition=cpuq
#SBATCH --mem=10GB
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --job-name=Coloc_net
#SBATCH --time=00:30:00

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_number> <R_script>" >&2
    exit 1
fi

input=$1
r_script=$(realpath "$2")

if [[ ! -f "$r_script" ]]; then
    echo "ERROR: R script not found: $r_script" >&2
    exit 1
fi

conda_env_dir=$(dirname "$CONDA_EXE")
conda_env_name="beyond_secretome"
source "$conda_env_dir/activate" "$conda_env_name"

results_dir=$(realpath "../results")
cd "$results_dir"

Rscript "$r_script" "$input"
