#!/bin/bash 
#SBATCH --job-name=preprocs
#SBATCH --partition=main
#SBATCH --time=48:00:00
#SBATCH --mail-user=rxxwang@umich.edu
#SBATCH --mail-type=end,fail
#SBATCH --nodes=1
#SBATCH --mem=7g
#SBATCH --output=/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/slurm_output/%x_%A_%a.out

Rscript --vanilla /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/1.1-Separate_data.R
