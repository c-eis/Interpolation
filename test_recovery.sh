#!/bin/bash
#SBATCH --job-name=int
#SBATCH -N 1
#SBATCH --time=02:00:00
#SBATCH --mail-user=cluettig@awi.de
#SBATCH --mail-type=ALL
module load gcc
module load R/3.3.0.gcc
module load gdal
module load GMT
module load python

#cp -r Recovery/outline/ .
cp Recovery/merge_shift_ref_vx.tif .
rm -r Recovery/
mkdir Recovery
#mv outline/ Recovery/
mv merge_shift_ref_vx.tif Recovery/
srun ./kombi.sh Recovery/merge_shift_ref_vx Recovery/int_merge_shift_ref_vx S
