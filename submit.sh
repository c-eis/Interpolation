#!/bin/bash
#
#SBATCH --array=0-123%100
#SBATCH --job-name=kriging
#SBATCH --output=slurm_%a.out
#SBATCH -p smp
#SBATCH --time=12:00:00

#SBATCH --mem=5000M

ulimit -s unlimited
module purge
module load gcc
module load python3
module load gdal/2.2.3
module load R/3.3.0
module load centoslibs

srun /global/AWIsoft/R/3.3.0/lib64/R/bin/Rscript --vanilla slurm_run.R

