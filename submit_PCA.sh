#!/bin/bash 
#SBATCH --job-name=EWAS
#SBATCH --partition=main
#SBATCH --time=48:00:00
#SBATCH --mail-user=rxxwang@umich.edu
#SBATCH --mail-type=end,fail
#SBATCH --nodes=1
#SBATCH --mem=7g
#SBATCH --array=1-100
#SBATCH --output=/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/slurm_output/%x_%A_%a.out

adbio=$1
model=$2
echo $adbio
echo "Model${model}"
if [[ $model =~ ^[1-6]$ ]]; then
  Rscript --vanilla "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/5.${model}-EWAS_model${model}.R" $SLURM_ARRAY_TASK_ID $adbio
else
  echo "Error: input must be a number between 1 and 6"
  exit 1
fi
