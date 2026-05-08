#!/bin/bash

#SBATCH --job-name=13lgs_cluster
#SBATCH --time=08:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/normalize_cluster_%j.out
#SBATCH --error=logs/normalize_cluster_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jflores@unc.edu

module purge
module load r/4.5.0

cd /work/users/j/p/jpflores/projects/13LGS_PilotAnalyses

Rscript scripts/processing/normalize_cluster.R