#!/bin/bash

#SBATCH --job-name=sctransform_per_sample
#SBATCH --time=04:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-8                    # one job per sample (8 total)
#SBATCH --output=logs/sct_%A_%a.out   # %A = job ID, %a = array index
#SBATCH --error=logs/sct_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jflores@unc.edu

## Load R module
module purge
module load r/4.5.0

## Change into the project directory first so renv activates automatically
cd /work/users/j/p/jpflores/projects/13LGS_PilotAnalyses

## Then run the script
Rscript scripts/processing/sctransform_per_sample.R