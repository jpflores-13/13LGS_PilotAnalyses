#!/bin/bash

#SBATCH --job-name=13lgs_doublets
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-8                       # one job per sample
#SBATCH --output=logs/doublets_%A_%a.out
#SBATCH --error=logs/doublets_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jflores@unc.edu

module purge
module load r/4.5.0

cd /work/users/j/p/jpflores/projects/13LGS_PilotAnalyses

Rscript scripts/processing/doublet_removal.R