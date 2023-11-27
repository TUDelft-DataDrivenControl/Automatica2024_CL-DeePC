#!/bin/sh
#
#SBATCH --job-name="d_Re"
#SBATCH --partition=compute
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=41
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --account=research-3me-dcsc
#SBATCH --output=d_Re.%j.out
#SBATCH --error=d_Re.%j.err

cd ${HOME}/../../scratch/${USER}/"Automatica 2023"/src/
module load matlab

matlab -r varying_noise_level
echo "Finished MATLAB calculations, varying Re."