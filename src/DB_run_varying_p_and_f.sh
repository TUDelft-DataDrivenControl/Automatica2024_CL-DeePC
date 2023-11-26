#!/bin/sh
#
#SBATCH --job-name="d_pf"
#SBATCH --partition=compute
#SBATCH --time=20:00:00
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --account=research-3me-dcsc
#SBATCH --output=d_pf.%j.out
#SBATCH --error=d_pf.%j.err

module load matlab

matlab -r varying_p_and_f
echo "Finished MATLAB calculations, varying p & f."