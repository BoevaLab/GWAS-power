#!/bin/bash

#SBATCH -p gpu
#SBATCH -c 16
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=6G
#SBATCH --gres=gpu:rtx4090:1
#SBATCH --job-name=aitl_loocv_%j
#SBATCH --output=aitl_loocv_%j.out

source ~/.bashrc
conda activate /cluster/work/boeva/lrabuzin/conda_envs/spatial

# Define the input file
INPUT_FILE="method_parameters_independent_set.txt"

# Define the Python script
PYTHON_SCRIPT="gwas_snp_selection_fdr_bonf_method.py"

# Check if the input file exists
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file '$INPUT_FILE' not found!"
    exit 1
fi

# Loop through each line in the input file and call the Python script
while IFS= read -r input; do
    echo "Processing input: $input"
    python3 "$PYTHON_SCRIPT" "$input"
done < "$INPUT_FILE"

echo "All inputs processed successfully."
