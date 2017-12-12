#!/bin/bash
##
## Time-stamp: <2016-08-05 11:11:17 (cluettig)>
##    
#SBATCH --job-name=interpolation
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --mail-user=cluettig@awi.de
#SBATCH --mail-type=END

ulimit -s unlimited
module purge
module load gcc
module load python
module load gdal

srun python int.py /work/ollie/cluettig/Greenland/data/Kriging_Test/vx.tif /work/ollie/cluettig/Greenland/data/Kriging_Test/vx_int.tif


