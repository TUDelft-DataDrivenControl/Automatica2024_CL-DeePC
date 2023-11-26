#!/bin/sh
#
#SBATCH --job-name="test_transfer"
#SBATCH --partition=trans
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=research-3me-dcsc
#SBATCH --output=DB_test_transfer.%j.out
#SBATCH --error=DB_test_transfer.%j.err

source="/scratch/${USER}/'Automatica 2023'/data/raw/"
destination="/tudelft.net/staff-umbrella/'RD Backup Drive'/'Contributed Papers'/journal/'Automatica 2023'/data/raw/"

rsync -av --no-perms "${source}" "${destination}"