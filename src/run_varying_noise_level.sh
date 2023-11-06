#!/bin/bash
#PBS -l nodes=1:ppn=10
#PBS -N log.varying_noise_level
module load matlab
cd "$HOME/Rogier/Automatica_2023/src"

dir

matlab  -r varying_noise_level
echo "Finished MATLAB  calculations. Duration (s):" 