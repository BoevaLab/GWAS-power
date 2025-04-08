#!/bin/bash
#SBATCH -c 16
#SBATCH --time=0:45:00
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=gwas_power_%A_%a
#SBATCH --output=gwas_power_%A_%a.out
#SBATCH --array=1-758

source ~/.bashrc
conda activate deepcast_gwas

# Define the input file and Python script
INPUT_FILE="method_icd_phenotypes.txt"
PYTHON_SCRIPT="gwas_snp_selection_fdr_bonf_method_updated.py"

# Check if the input file exists
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file '$INPUT_FILE' not found!"
    exit 1
fi

# Get total number of lines in the input file
TOTAL_LINES=$(wc -l < "$INPUT_FILE")
if [[ $SLURM_ARRAY_TASK_ID -gt $TOTAL_LINES ]]; then
    echo "Error: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) is greater than the number of lines in $INPUT_FILE ($TOTAL_LINES)."
    exit 1
fi

# Read the line corresponding to the current array task
input=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INPUT_FILE")
echo "Processing input: $input"
python3 "$PYTHON_SCRIPT" $input

echo "Input processed successfully."