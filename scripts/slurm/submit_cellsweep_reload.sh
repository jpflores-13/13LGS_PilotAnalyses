#!/bin/bash

#SBATCH --job-name=13lgs_cellsweep_reload
#SBATCH --time=06:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/cellsweep_reload_%j.out
#SBATCH --error=logs/cellsweep_reload_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jflores@unc.edu

module purge
module load r/4.5.0

cd /work/users/j/p/jpflores/projects/13LGS_PilotAnalyses

Rscript scripts/processing/cellsweep_reload.R