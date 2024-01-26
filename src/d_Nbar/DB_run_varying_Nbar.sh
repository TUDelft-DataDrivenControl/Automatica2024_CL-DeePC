#!/bin/bash

# Function to submit SLURM job array
submit_job_array() {
    local index=$1

    # Generate a unique SLURM submission script for each set of parameters
    script_name="./Job${index}/submit_job_array_${index}.sh"

    # Create the SLURM submission script
    cat > "$script_name" <<EOL
#!/bin/bash
#SBATCH --job-name=d_Nbar${index}
#SBATCH --partition=compute
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --account=research-3me-dcsc
#SBATCH --output=./Job${index}/d_Nbar.%A_%a.out
#SBATCH --error=./Job${index}/d_Nbar.%A_%a.err
#SBATCH --array=1-120

# Load any necessary modules or set environment variables
module load matlab

# Your commands or script for each job array
matlab -nosplash -nodesktop -r "varying_Nbar(${index},\$SLURM_ARRAY_TASK_ID),quit"
echo "Finished MATLAB calculations."

EOL

    # Make the script executable
    chmod +x "$script_name"

    # Submit the SLURM job
    sbatch "$script_name"

    echo "Submitted job array ${index}"
}

# Specify the range of job arrays
start_index=1
end_index=50

cd ${HOME}/../../scratch/${USER}/Automatica_2023/src/d_Nbar

# Iterate over the range of job arrays
for index in $(seq $start_index $end_index); do
    # Call the function to submit the SLURM job array
    mkdir "./Job${index}"
    submit_job_array "$index"
done
echo "Submitted all Slurm job arrays"