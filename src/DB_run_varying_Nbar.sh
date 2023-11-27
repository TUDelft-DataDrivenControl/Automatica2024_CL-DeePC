#!/bin/sh
#
#SBATCH --job-name="d_Nbar"
#SBATCH --partition=compute
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=41
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --account=research-3me-dcsc
#SBATCH --output=d_Nbar.%j.out
#SBATCH --error=d_Nbar.%j.err

module load matlab

matlab -r varying_Nbar
echo "Finished MATLAB calculations, varying Nbar."