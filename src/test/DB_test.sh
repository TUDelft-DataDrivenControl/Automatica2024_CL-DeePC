#!/bin/sh
#
#SBATCH --job-name="DB_test"
#SBATCH --partition=compute
#SBATCH --time=00:30:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --account=research-3me-dcsc
#SBATCH --output=DB_test.%j.out
#SBATCH --error=DB_test.%j.err

module load matlab
module list

dir
cd ${HOME}/../../scratch/rogierdinkla/"Automatica 2023"/src/

matlab -r test_script
echo "Finished MATLAB calculations."
