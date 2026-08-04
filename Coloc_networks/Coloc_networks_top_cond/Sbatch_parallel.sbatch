#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@fht.org
#SBATCH --output="logs/%j_Node.log"
#SBATCH --partition=cpuq
#SBATCH --mem=10GB
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --job-name=Coloc_net
#SBATCH --time=00:30:00

source /center/healthds/singularity_functions
cd $HOME

input=$1
r_script=$2

Rscript "${r_script}" "${input}"
