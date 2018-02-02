#!/bin/bash
##
## Time-stamp: <2016-08-05 11:11:17 (cluettig)>
##
#SBATCH --job-name=interpolation
#SBATCH -N 1
#SBATCH --time=00:10:00
#SBATCH --mail-user=cluettig@awi.de
#SBATCH --mail-type=END

ulimit -s unlimited
module purge
module load gcc
module load python3
module load gdal/2.2.3
module load R/3.3.0
module load centoslibs

srun python3 int_onlyR.py /work/ollie/cluettig/Recovery/gefiltert/vx.tif /work/ollie/\
cluettig/Recovery/gefiltert/vx_int.tif
