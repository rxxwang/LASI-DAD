#!/bin/bash 
#SBATCH --job-name=cor_filt
#SBATCH --partition=main
#SBATCH --time=48:00:00
#SBATCH --mail-user=rxxwang@umich.edu
#SBATCH --mail-type=end,fail
#SBATCH --nodes=1
#SBATCH --mem=15g
#SBATCH --array=1-22
#SBATCH --output=/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/slurm_output/%x_%A_%a.out

Rscript --vanilla /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/3.0-cor_filter.R $SLURM_ARRAY_TASK_ID

