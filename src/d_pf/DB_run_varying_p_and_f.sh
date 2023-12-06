#!/bin/sh
#
#SBATCH --job-name="sched_d_pf"
#SBATCH --partition=compute
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --account=research-3me-dcsc
#SBATCH --output=scheduler_d_pf.%j.out
#SBATCH --error=scheduler_d_pf.%j.err
#SBATCH --mail-user=r.t.o.dinkla@tudelft.nl
#SBATCH --mail-type=ALL

cd ${HOME}/../../scratch/${USER}/Automatica_2023/src/d_pf
module load matlab

matlab -r varying_p_and_f #>> scheduler_d_pf.${SLURM_JOB_ID}.out
echo "Finished MATLAB calculations, varying p & f."