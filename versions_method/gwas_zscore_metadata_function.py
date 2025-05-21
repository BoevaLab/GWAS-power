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
import concurrent.futures
from functools import partial
import subprocess
import tempfile
import shutil
from datetime import datetime
from pathlib import Path

# Author; Sophie Sigfstead
# Purpose: Latest version of our filtering method to improve GWAS power. 


### CONSTANTS 
R2_THRESHOLD = 0.2    # LD r² threshold for clumping
KB_RADIUS = 500       # Distance in kb to look for LD (similar to WINDOW_SIZE but in kb)

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

def window_elimination(df, WINDOW_SIZE):
    """
    Input: Dataframe of summary statistics under the p-value threshold
    Output: Dataframe containing only leading SNPs (SNPs within leading SNP position +/- WINDOW_SIZE are eliminated)
    """
    result_df = pd.DataFrame()
    for chr_num in range(1, 23):
        chr_df = df[df['chr'] == chr_num].sort_values(by='p_value')
        while not chr_df.empty:
            top_snp = chr_df.iloc[0]
            result_df = result_df._append(top_snp, ignore_index=True)
            chr_df = chr_df[~((chr_df['pos'] - top_snp['pos']).abs() < WINDOW_SIZE)].reset_index(drop=True)
    return result_df 

def ld_based_clumping(df, ld_lookup_path, p_value_column='p_value', snp_column='snp', p_threshold=1):
    """
    Perform LD-based clumping using PLINK.
    
    Parameters:
        df (pd.DataFrame): Dataframe containing SNP information and p-values
        p_value_column (str): Column name for p-values
        snp_column (str): Column name for SNP IDs
        p_threshold (float): Significance threshold for index SNPs
        
    Returns:
        pd.DataFrame: Clumped SNP list with borders based on actual LD relationships
    """

    global R2_THRESHOLD
    global KB_RADIUS

    print(f"Performing LD-based clumping with lookup table (r²={R2_THRESHOLD}, window={KB_RADIUS})", flush=True)
    
    # Check if input DataFrame is empty or has no SNPs below threshold
    filtered_df = df[df[p_value_column] <= p_threshold].copy()
    if filtered_df.empty:
        print(f"Warning: No SNPs pass the p-value threshold of {p_threshold}", flush=True)
        return pd.DataFrame()
    
    print(f"Input has {len(filtered_df)} SNPs passing p-value threshold", flush=True)
    print("Performing ld-based border lookup", flush=True)

    df_windows = pd.read_pickle(ld_lookup_path)

    # 3a) Merge by index
    df_merged = filtered_df.rename(columns={snp_column:'SNP'}).merge(
        df_windows,
        how='left',           # keep all rows in df_main
        left_on='SNP',        # column in df_main
        right_index=True      # index in df_windows
    )

    print(f"Finished merge with ld-based border lookup table.", flush=True)

    return df_merged

def process_csv_file(csv_file, directory_1000_genomes, sad_columns, pattern, summary_stats_by_chr):
    """
    Process a single CSV file:
      - Reads the CSV and extracts the chromosome using the precompiled regex pattern.
      - Retrieves the corresponding summary stats from summary_stats_by_chr.
      - Sets the index on the merge keys and performs a join.
      - Returns the merged DataFrame for that CSV.
    """
    csv_file_path = os.path.join(directory_1000_genomes, csv_file)
    try:
        chunk = pd.read_csv(csv_file_path, usecols=['chr', 'pos', 'ref', 'alt', 'snp'] + sad_columns)
    except Exception as e:
        print(f"Error reading {csv_file_path}: {e}", flush=True)
        return None

    match = pattern.search(csv_file)
    if not match:
        return None
    chromosome = int(match.group(1))
    
    df_summary_stats_chr = summary_stats_by_chr.get(chromosome)
    if df_summary_stats_chr is None or df_summary_stats_chr.empty:
        return None
    
    # Set index on chunk to match summary stats index
    chunk = chunk.set_index(['chr', 'pos', 'ref', 'alt'])
    merged_result = df_summary_stats_chr.join(chunk, how='left').reset_index()
    
    # Keep only rows with valid SAD and p_value values
    merged_result = merged_result[merged_result[sad_columns[0]].notna() & merged_result['p_value'].notna()]
    print(f"----Processed file {csv_file_path}", flush=True)
    return merged_result

def extract_metadata(
    result_df: pd.DataFrame,
    summary_stats_df: pd.DataFrame,
    p_value_threshold: float
) -> pd.DataFrame:
    """
    Compute metadata statistics for SNP results and return as a single-row DataFrame.

    Parameters:
    - result_df: DataFrame with SNP results, must contain 'snp' column and optionally 'in_coding_region'.
    - summary_stats_df: DataFrame of original significant SNPs, must contain 'snp' and optionally 'in_coding_region'.
    - p_value_threshold: threshold used for significance.

    Returns:
    - metadata_df: one-row DataFrame with computed metadata.
    """
    # Initialize counts
    num_snps_found = len(result_df)
    num_coding_snps_found = 0
    num_overlapping_snps = 0
    num_overlapping_loci = 0
    num_original_list = len(summary_stats_df)
    num_original_coding_snps = 0
    percentage_loci_recovered = 1.0

    # Only compute overlaps if results exist
    if not result_df.empty:
        if 'in_coding_region' in result_df.columns:
            num_coding_snps_found = int(result_df['in_coding_region'].sum())

        if not summary_stats_df.empty:
            # SNP-level overlap
            overlap_df = pd.merge(result_df, summary_stats_df, on='snp')
            num_overlapping_snps = len(overlap_df)

            # Locus-level overlap via provided function
            num_overlapping_loci = overlaps(summary_stats_df, result_df)

            # Original coding SNPs
            if 'in_coding_region' in summary_stats_df.columns:
                num_original_coding_snps = int(summary_stats_df['in_coding_region'].sum())

            if num_original_list > 0:
                percentage_loci_recovered = num_overlapping_loci / num_original_list

    global R2_THRESHOLD
    global KB_RADIUS

    # Build metadata dictionary
    metadata = {
        'num_snps_found': num_snps_found,
        'num_coding_snps_found': num_coding_snps_found,
        'num_overlapping_snps': num_overlapping_snps,
        'num_overlapping_loci': num_overlapping_loci,
        'num_original_list': num_original_list,
        'num_original_coding_snps': num_original_coding_snps,
        'p_value_threshold': p_value_threshold,
        'percentage_loci_recovered': percentage_loci_recovered,
        'r2_threshold': R2_THRESHOLD,
        'kb_radius': KB_RADIUS
    }

    print(f"   - Final results: {num_snps_found} SNPs identified", flush=True)
    if num_original_list > 0:
        print(f"   - Recovered {num_overlapping_loci}/{num_original_list} loci ({percentage_loci_recovered:.1%})", flush=True)

    return pd.DataFrame([metadata])


def main(df_summary_stats, directory_1000_genomes, track_list, coding_snp_list_path, ld_lookup_path, alpha, ld_based=True, plink_path=None, ld_reference_path=None):
    """
    1. Process GWAS summary stats and merge with SAD track information from 1000 Genomes files.
    2. Perform filtering of GWAS SNPs, returning final results list.
    
    Parameters:
        df_summary_stats (pd.DataFrame): Summary stats DataFrame.
        directory_1000_genomes (str): Directory containing 1000 Genomes CSV files.
        track_list (list): List of track identifiers for SAD columns.
        alpha (float): The threshold for the fdr procedure
        ld_based (bool): Whether to use LD-based clumping (True) or window-based elimination (False).
        plink_path (str): Path to PLINK executable if LD-based clumping is used.
        ld_reference_path (str): Path to LD reference panel if LD-based clumping is used.
    """
    print("====== Z-score based GWAS SNP Selection ======", flush=True)    

    ## STEP 1: Process summary statistics
    print("1. Processing summary statistics", flush=True)
    # Reverse the -log(pvalue) operation and drop the original column
    df_summary_stats['p_value'] = 10 ** (-df_summary_stats['neglog10_pval_EUR'])
    df_summary_stats.drop(columns=['neglog10_pval_EUR'], inplace=True)
    df_summary_stats.dropna(subset=['p_value'], inplace=True)
    df_summary_stats.reset_index(drop=True, inplace=True)
    print(f"   - Processed {len(df_summary_stats)} SNPs from summary statistics", flush=True)

    ## STEP 2: Merge with SAD data
    print(f"2. Merging with SAD values", flush=True)
    # Initialize SAD columns with track names
    sad_columns = [f"SAD{track}" for track in track_list]

    # Pre-group summary stats by chromosome for faster joins
    summary_stats_by_chr = {
        chr_val: sub_df.copy().set_index(['chr', 'pos', 'ref', 'alt'])
        for chr_val, sub_df in df_summary_stats.groupby('chr')
    }

    # Extract chromosome from filename with regex
    pattern = re.compile(r'\.MAF_threshold=0\.005\.(\d+)_combined\.csv')
    csv_files = [f for f in os.listdir(directory_1000_genomes) if f.endswith('.csv') and f != 'targets.csv']
    
    # Process files in parallel
    dataframes = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        func = partial(process_csv_file,
                      directory_1000_genomes=directory_1000_genomes,
                      sad_columns=sad_columns,
                      pattern=pattern,
                      summary_stats_by_chr=summary_stats_by_chr)
        results = list(executor.map(func, csv_files))
    
    # Combine results
    for res in results:
        if res is not None and not res.empty:
            dataframes.append(res)
            
    if not dataframes:
        print("Error: No SAD data could be merged with summary statistics", flush=True)
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        
    df_summary_stats_result = pd.concat(dataframes, ignore_index=True)
    print(f"   - Merged data contains {len(df_summary_stats_result)} SNPs", flush=True)

    ## STEP 3: Identify coding SNPs
    print(f"3. Identifying coding SNPs", flush=True)
    coding_region_set = list(pd.read_csv(coding_snp_list_path)['snp'])
    df_summary_stats_result['in_coding_region'] = df_summary_stats_result['snp'].isin(coding_region_set)
    coding_snp_count = df_summary_stats_result['in_coding_region'].sum()
    print(f"   - Found {coding_snp_count} coding SNPs", flush=True)
   
    ## STEP 4: Find SNPs with SAD scores in the relevant SAD score range:
    print(f"4. Find SNPs with SAD values in relevant range", flush=True)
    track_data = df_summary_stats_result[sad_columns].to_numpy()
    mean = np.mean(track_data)
    sd = np.std(track_data)
    
    in_sad_range_mask = ((track_data < mean - alpha * sd) + (track_data > mean + alpha * sd)).any(axis=1)
    df_summary_stats_result['sad_relevant'] = in_sad_range_mask
    print(f"   - Found {in_sad_range_mask.sum()} SNPs with SAD score outside of {alpha} * standard deviation(s) from the mean.", flush=True)
    
    # Create the reference list
    df_summary_stats_signifcant_list = df_summary_stats_result[df_summary_stats_result['p_value'] < 5e-8].reset_index(drop=True)
    original_snps_not_empty = not df_summary_stats_signifcant_list.empty
    
    if original_snps_not_empty:
        print(f"   - Reference list contains {len(df_summary_stats_signifcant_list)} SNPs at p < 5e-8", flush=True)
    else:
        print("   - No SNPs meet genome-wide significance (p < 5e-8) in reference list", flush=True)

    ## STEP 5: Compute significance threshold
    print(f"5. Computing adjusted significance threshold", flush=True)
    # Get the significant SNP list
    n_snps_prev = len(df_summary_stats_result)
    method_relevant_snps = df_summary_stats_result[df_summary_stats_result['sad_relevant'] | df_summary_stats_result['in_coding_region']]
    n_snps_curr = len(method_relevant_snps)
    
    # Calculate new p-value threshold
    p_value_threshold = 5e-8 * (n_snps_prev / n_snps_curr) if n_snps_curr > 0 else 5e-8
    print(f"   - Adjusted p-value threshold: {p_value_threshold:.2e}", flush=True)

    ## STEP 6: Identify lead SNPs
    print("6. Identifying lead SNPs", flush=True)
    
    # Process data according to the selected method
    if ld_based:
        print(f"   - Using LD-based clumping with p < {p_value_threshold:.2e}", flush=True)
        result_df = ld_based_clumping(
            method_relevant_snps,
            ld_lookup_path,
            p_threshold=p_value_threshold
        )
    
    # Process the reference list if it's not empty
    if original_snps_not_empty:
        print(f"7. Processing reference SNP list", flush=True)
        if ld_based:
            df_summary_stats_signifcant_list = ld_based_clumping(
                df_summary_stats_signifcant_list,
                ld_lookup_path,
                p_threshold=5e-8
            )
        else:
            df_summary_stats_signifcant_list = window_elimination(df_summary_stats_signifcant_list, WINDOW_SIZE)
            if not df_summary_stats_signifcant_list.empty:
                df_summary_stats_signifcant_list['left_border'] = df_summary_stats_signifcant_list['pos'] - WINDOW_SIZE
                df_summary_stats_signifcant_list['right_border'] = df_summary_stats_signifcant_list['pos'] + WINDOW_SIZE
    
    ## STEP 7: Compute metadata
    print(f"8. Computing result statistics", flush=True)
    
    metadata_df = extract_metadata(result_df, df_summary_stats_signifcant_list, p_value_threshold)
    
    # Prepare final output dataframes
    cols_to_keep = ['chr', 'pos', 'ref', 'alt', 'snp', 'p_value', 'in_coding_region']
    
    if not result_df.empty:
        available_cols = [col for col in cols_to_keep if col in result_df.columns]
        result_df = result_df[available_cols]
    
    if not df_summary_stats_signifcant_list.empty:
        available_cols = [col for col in cols_to_keep if col in df_summary_stats_signifcant_list.columns]
        df_summary_stats_signifcant_list = df_summary_stats_signifcant_list[available_cols]
    
    return result_df, metadata_df, df_summary_stats_signifcant_list

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
    directory_1000_genomes = sys.argv[2]
    coding_snp_list_path = sys.argv[3]
    track_list = sys.argv[4] if len(sys.argv) > 4 else None
    
    # Default to LD-based clumping if enough parameters are provided
    ld_based = True if len(sys.argv) > 6 else False
    plink_path = sys.argv[5] if len(sys.argv) > 5 else None
    ld_reference_path = sys.argv[6] if len(sys.argv) > 6 else None
    alpha = float(sys.argv[7]) if len(sys.argv) > 7 else 1

    ld_lookup_path = 'ld_lists/snp_windows.pkl'

    selected_columns = ['chr', 'pos', 'ref', 'alt', 'neglog10_pval_EUR']

    print("Opening sumstats", flush=True)
    file = gzip.open(file_path_summary_stats, "rt") if file_path_summary_stats.endswith(".tsv.bgz") else open(file_path_summary_stats, "r")
    
    print("Create df for sumstats", flush=True)
    df_summary_stats = pd.concat(pd.read_csv(file, sep="\t", usecols=selected_columns, chunksize=100000), ignore_index=True)
    
    if track_list is None:
        track_list = list(range(40))
    else:
        track_list_file = pd.read_csv(track_list, header=None)
        track_list = track_list_file.iloc[:, 0].tolist()
    
    result_df, metadata_df, df_summary_stats_signifcant_list = main(
        df_summary_stats, 
        directory_1000_genomes, 
        track_list, 
        coding_snp_list_path,
        ld_lookup_path,
        alpha=alpha,
        ld_based=ld_based,
        plink_path=plink_path,
        ld_reference_path=ld_reference_path
    )
    
    # Save results in deepcast_phenotypes folder under a subfolder for the phenotype.
    # We use the base name of the summary stats file (without extension) as the phenotype name.
    phenotype_name = os.path.splitext(os.path.basename(file_path_summary_stats))[0]
    output_folder = Path("deepcast_phenotypes/results_1_sd") / phenotype_name
    output_folder.mkdir(parents=True, exist_ok=True)
    
    result_df.to_csv(os.path.join(output_folder, "result_df.csv"), index=False)
    metadata_df.to_csv(os.path.join(output_folder, "metadata_df.csv"), index=False)
    df_summary_stats_signifcant_list.to_csv(os.path.join(output_folder, "df_summary_stats_signifcant_list.csv"), index=False)
    
    print(f"Results saved in folder: {output_folder}")