import random
import subprocess
import os
import pandas as pd 
import matplotlib.pyplot as plt
import shutil
import warnings
import sys
warnings.filterwarnings("ignore", category=FutureWarning)

def main():
    
    window_sizes = [250000, 375000, 500000, 750000]

    brain_tracks_files = "../gwas_1_matching"

    for window in window_sizes:
        # run gwas_1_leading_snps_procedue with 0.0 standard deviation and window size given
        command= [
            'python', 
            './gwas_1_leading_snps_window_size.py',
            f"{brain_tracks_files}",  
            f'snp_lists_results/id=reference/window={window}/filtered_snps_gwas_1_sd=0.0.csv', 
            f"0.0",
            f"{window}",
            f"combined_results/id=reference/window={window}/combined_result_sd=0.0.csv",
            f'snp_lists_results/id=reference/window={window}/filtered_snps_gwas_1_sd=0.0.csv'
            ]
        
        subprocess.run(command, check=True)

        # using the new reference result as the reference set, run gwas_1_leading_snps_procedure with the brain tracks 
        reference_result = f'snp_lists_results/reference/filtered_snps_gwas_1_threshold=0.0.csv'

        command=[
            'python', 
            './gwas_1_leading_snps_window_size.py',
            f"{brain_tracks_files}",  
            f'{reference_result}', 
            f"1.0",
            f"{window}",
            f"combined_results/id=reference/window={window}/combined_result_sd=1.0.csv",
            f'snp_lists_results/id=reference/window={window}/filtered_snps_gwas_1_sd=1.0.csv'
            ]

        subprocess.run(command, check=True)
        # using the new reference result as the ref set, run gwas_1_leading_snps_procedure with the random sets
        random_sets = "./gwas_1_matching_random_sets/id="
        keys = [chr(i) for i in range(ord('a'), ord('a') + 19)]
        for key in keys:
            matched_files_path = os.path.join(random_sets,key)
            print(matched_files_path)
            command=[
            'python', 
            './gwas_1_leading_snps_window_size.py',
            f"{matched_files_path}",  
            f'{reference_result}', 
            f"1.0",
            f"{window}",
            f"combined_results/id={key}/window={window}/combined_result_sd=1.0.csv",
            f'snp_lists_results/id={key}/window={window}/filtered_snps_gwas_1_sd=1.0.csv'
            ]
            subprocess.run(command, check=True)

    return


if __name__ == "__main__":
    main()