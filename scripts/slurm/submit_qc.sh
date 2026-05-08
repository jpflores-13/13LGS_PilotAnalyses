#!/bin/bash

#SBATCH --job-name=13lgs_qc
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/qc_filtering_%j.out
#SBATCH --error=logs/qc_filtering_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jflores@unc.edu

module purge
module load r/4.5.0

cd /work/users/j/p/jpflores/projects/13LGS_PilotAnalyses

Rscript scripts/processing/qc_filtering.R
