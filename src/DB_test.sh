#!/bin/sh
#
#SBATCH --job-name="d_Nbar"
#SBATCH --partition=compute
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --account=research-3me-dcsc
#SBATCH --output=d_Nbar.%j.out
#SBATCH --error=d_Nbar.%j.err

module load matlab

matlab -r test_script
echo "Finished MATLAB calculations."
