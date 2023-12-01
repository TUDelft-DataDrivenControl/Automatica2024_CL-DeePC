#!/bin/sh
#
#SBATCH --job-name="d_pf"
#SBATCH --partition=compute
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=41
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --account=research-3me-dcsc
#SBATCH --output=d_pf.%j.out
#SBATCH --error=d_pf.%j.err
#SBATCH --mail-user=r.t.o.dinkla@tudelft.nl
#SBATCH --mail-type=ALL

cd ${HOME}/../../scratch/${USER}/Automatica_2023/src/d_pf
module load matlab

matlab -r varying_p_and_f
echo "Finished MATLAB calculations, varying p & f."