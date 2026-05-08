#!/bin/bash

#SBATCH --job-name=13lgs_marker_plots
#SBATCH --time=02:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/marker_plots_%j.out
#SBATCH --error=logs/marker_plots_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jflores@unc.edu

module purge
module load r/4.5.0

cd /work/users/j/p/jpflores/projects/13LGS_PilotAnalyses

Rscript scripts/plots/marker_plots.R
