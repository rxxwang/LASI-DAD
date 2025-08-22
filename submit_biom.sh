#!/bin/bash 
#SBATCH --job-name=plotbiom
#SBATCH --partition=main
#SBATCH --time=48:00:00
#SBATCH --mail-user=rxxwang@umich.edu
#SBATCH --mail-type=end,fail
#SBATCH --nodes=1
#SBATCH --mem=15g
#SBATCH --array=1-5
#SBATCH --output=/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/slurm_output/%x_%A_%a.out

model=$1
echo "Model${model}"
Rscript --vanilla /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/6.0-qqplot.R $SLURM_ARRAY_TASK_ID $model
