#!/bin/bash

#SBATCH --job-name=13lgs_cellsweep
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-8
#SBATCH --output=logs/cellsweep_%A_%a.out
#SBATCH --error=logs/cellsweep_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jflores@unc.edu

module purge
module load python/3.12.4   # adjust to your cluster's Python module

cd /work/users/j/p/jpflores/projects/13LGS_PilotAnalyses

python scripts/processing/cellsweep_convert.py
python scripts/processing/cellsweep_run.py