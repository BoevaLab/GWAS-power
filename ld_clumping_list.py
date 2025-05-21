import os
import glob
import pandas as pd
from scipy.stats import ttest_1samp
from statsmodels.stats.multitest import multipletests
import numpy as np
import sys
import warnings
warnings.filterwarnings("ignore")
import gzip
import re
from intervaltree import IntervalTree
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import subprocess
import tempfile
import shutil
from datetime import datetime

# Author; Sophie Sigfstead
# Purpose: Latest version of our filtering method to improve GWAS power. 

def overlaps(df1, df2):
    """
    Input: Two dataframes with snp pos, left border and right border. 
           - df1 should be the original/reference list
           - df2 should be the new list
    Output: The number of overlapping loci between df1 and df2
    """
    overlap_count = 0
    # Iterate through each interval in df1
    for i, row1 in df1.iterrows():
        left1, right1 = row1['left_border'], row1['right_border']
        chr1  = row1['chr']
        for j, row2 in df2.iterrows():
            left2, right2 = row2['left_border'], row2['right_border']
            chr2  = row2['chr']
            if (left1 <= right2 and right1 >= left2) and (chr1 == chr2):
                overlap_count += 1
                break  # Exit loop once an overlap is found
    return overlap_count



# takes everything above r^2 0.2 that is within 500kb
#  
def ld_based_clumping(df, plink_path, ld_reference_path, p_value_column='p_value', snp_column='snp', chr_column='chr', bp_column='pos', r2_threshold=0.2, kb_radius=500, p_threshold=1):
    # Create temporary directory for intermediate files
    temp_dir = tempfile.mkdtemp()

    # Prepare input file for PLINK
    target_file = os.path.join(temp_dir, "targets.txt")
    clumped_output = os.path.join(temp_dir, "plink_output")
    
    # Format the association file for PLINK
    # Extract the first column (assumed to be the SNP ID) and drop duplicates
    target_ids = df[snp_column]

    # Write out to targets.txt (one SNP ID per line, no header)
    target_ids.to_csv(target_file, index=False, header=False)
    
    # Run PLINK clumping
    print(f"Running PLINK clumping", flush=True)
    cmd_str = f"{plink_path} --bfile {ld_reference_path} --show-tags {target_file}" \
                f"    --tag-r2 {r2_threshold} --tag-kb {kb_radius} --list-all"
    
    try:
        subprocess.run(cmd_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except Exception as e:
        print(f"Error executing PLINK: {e}", flush=True)
        return pd.DataFrame()
        
    #     # Process PLINK results
    #     print(f"Processing PLINK results", flush=True)
    #     clumped_file = f"{clumped_output}.clumped"
    #     if not os.path.exists(clumped_file) or os.path.getsize(clumped_file) == 0:
    #         print(f"Warning: PLINK clumping produced no results", flush=True)
    #         return pd.DataFrame()
        
    #     clumped_results = pd.read_csv(clumped_file, delim_whitespace=True)
        
    #     # Set default borders based on kb_radius
    #     lead_snps = clumped_results['SNP'].tolist()
    #     result_df = df[df[snp_column].isin(lead_snps)].copy()
    #     result_df['left_border'] = result_df[bp_column] - (kb_radius * 1000)
    #     result_df['right_border'] = result_df[bp_column] + (kb_radius * 1000)
        
    #     # Get precise LD borders for each lead SNP
    #     print(f"Getting LD-based borders for each lead SNP ({len(result_df)} lead SNPs)...", flush=True)
    #     ld_count = 0
        
    #     for i, row in result_df.iterrows():
    #         lead_snp = row[snp_column]
    #         chr_val = row[chr_column]
    #         pos_val = row[bp_column]
            
    #         # Run PLINK to get LD information
    #         tags_output = os.path.join(temp_dir, f"tags_{lead_snp}")
    #         tags_cmd = f"{plink_path} --bfile {ld_reference_path} --r2 --ld-snp {lead_snp} " \
    #                   f"--ld-window-kb {kb_radius} --ld-window 999999 --ld-window-r2 {r2_threshold} " \
    #                   f"--out {tags_output}"
            
            
    #         try:
    #             subprocess.run(tags_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    #             ld_file = f"{tags_output}.ld"
                
    #             if os.path.exists(ld_file) and os.path.getsize(ld_file) > 0:
    #                 # Process LD results
    #                 ld_data = pd.read_csv(ld_file, delim_whitespace=True)
    #                 if not ld_data.empty:
    #                     # Find min and max positions of SNPs in LD
    #                     min_pos = float('inf')
    #                     max_pos = 0
    #                     for _, ld_row in ld_data.iterrows():
    #                         if ld_row['R2'] >= r2_threshold:
    #                             bp = ld_row['BP_B']
    #                             min_pos = min(min_pos, bp)
    #                             max_pos = max(max_pos, bp)
                        
    #                     # Only update if we found LD SNPs
    #                     if min_pos < float('inf') and max_pos > 0:
    #                         result_df.at[i, 'left_border'] = min_pos
    #                         result_df.at[i, 'right_border'] = max_pos
    #                         ld_count += 1
    #         except Exception as e:
    #             print(f"Error processing LD for SNP {lead_snp}: {e}", flush=True)
    #             # Continue with the next SNP - we'll use default borders for this one
            
    #     print(f"Updated LD borders for {ld_count} lead SNPs", flush=True)
    #     return result_df
            
    # except Exception as e:
    #     print(f"Error during LD-based clumping: {e}", flush=True)
    #     return pd.DataFrame()
    # finally:
    #     # Clean up temporary directory
    #     try:
    #         shutil.rmtree(temp_dir)
    #     except:
    #         pass  # Ignore cleanup errors

def main(df_summary_stats, directory_1000_genomes, track_list, coding_snp_list_path, alpha, ld_based=True, plink_path=None, ld_reference_path=None):
    ### CONSTANTS 
    R2_THRESHOLD = 0.2    # LD r² threshold for clumping
    KB_RADIUS = 500       # Distance in kb to look for LD (similar to WINDOW_SIZE but in kb)
    N_CPUS = os.cpu_count() or 1

    # 4 roughly equal chunks
    chunks = np.array_split(df_summary_stats, N_CPUS)

    with ProcessPoolExecutor(max_workers=N_CPUS) as exec:
        func = partial(ld_based_clumping, 
            plink_path = plink_path, 
            ld_reference_path = ld_reference_path,
            r2_threshold=R2_THRESHOLD,
            kb_radius=KB_RADIUS
        )
        processed = list(exec.map(func, chunks))
    df_result = pd.concat(processed, ignore_index=True)
    
    # Process data according to the selected method
    if ld_based:
        print(f"   - Using LD-based clumping with p < {p_value_threshold:.2e}", flush=True)
        result_df = ld_based_clumping(
            method_relevant_snps, 
            plink_path, 
            ld_reference_path,
            p_threshold=p_value_threshold,
            r2_threshold=R2_THRESHOLD,
            kb_radius=KB_RADIUS
        )

if __name__ == "__main__":
    """
    Inputs to main: 
        
        file_path_summary_stats: path to summary statistics file (.tsv or .tsv.bgz)
        directory_1000_genomes: directory housing the SAD score CSV files for all SNPs
        coding_snp_list_path: path to coding SNP dataframe (must have a 'snp' column)
        track_list: (optional) CSV file (no header) with one column of track numbers. If not provided, a default list is used.
        ld_based: (optional) boolean flag to use LD-based clumping (default: True)
        plink_path: (optional) path to PLINK executable (required if ld_based=True)
        ld_reference_path: (optional) path to LD reference panel (required if ld_based=True)
    
    Outputs: 
        result_df: final list of SNPs (without SAD values or FDR t-test info)
        metadata_df: metadata about the filtering process
    
    How to run:
        python gwas_snp_selection_fdr_bonf_method.py <file_path_summary_stats> <directory_1000_genomes> <coding_snp_list_path> [<track_list>] [<ld_based>] [<plink_path>] [<ld_reference_path>]
    """
    if len(sys.argv) < 4:
        print("Usage: python gwas_snp_selection_fdr_bonf_method.py <file_path_summary_stats> <directory_1000_genomes> <coding_snp_list_path> [<track_list>] [<ld_based>] [<plink_path>] [<ld_reference_path>]")
        sys.exit(1)
        
    file_path_summary_stats = sys.argv[1]
    plink_path = sys.argv[2]
    ld_reference_path = sys.argv[3]

    print("Opening sumstats", flush=True)
    file = gzip.open(file_path_summary_stats, "rt") if file_path_summary_stats.endswith(".tsv.bgz") else open(file_path_summary_stats, "r")
    
    print("Create df for sumstats", flush=True)
    df_summary_stats = pd.concat(pd.read_csv(file, sep="\t", usecols=selected_columns, chunksize=100000), ignore_index=True)

    main(
        df_summary_stats,
        plink_path=plink_path,
        ld_reference_path=ld_reference_path
    )
    