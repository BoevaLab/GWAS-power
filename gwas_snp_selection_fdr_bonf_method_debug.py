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

def main(df_summary_stats, directory_1000_genomes, track_list, coding_snp_list_path):
    """
    1. Process GWAS summary stats and merge with SAD track information from 1000 Genomes files.
    2. Perform filtering of GWAS SNPs, returning final results list.
    
    Parameters:
        df_summary_stats (pd.DataFrame): Summary stats DataFrame.
        directory_1000_genomes (str): Directory containing 1000 Genomes CSV files.
        track_list (list): List of track identifiers for SAD columns.
    """
    ### CONSTANTS 
    ALPHA = 0.05
    WINDOW_SIZE = 500000  # Range of bp eliminated around a leading SNP

    ## STEP 1: Merge the summary statistics with the SAD value data
    print("Merging sumstats with SAD values", flush=True)

    # Reverse the -log(pvalue) operation and drop the original column
    df_summary_stats['p_value'] = 10 ** (-df_summary_stats['neglog10_pval_EUR'])
    df_summary_stats.drop(columns=['neglog10_pval_EUR'], inplace=True)
    df_summary_stats.dropna(subset=['p_value'], inplace=True)
    df_summary_stats.reset_index(drop=True, inplace=True)

    # Initialize SAD columns with NaN
    sad_columns = [f"SAD{track}" for track in track_list]

    # Pre-group summary stats by chromosome and set index on merge keys for faster joins
    summary_stats_by_chr = {
        chr_val: sub_df.copy().set_index(['chr', 'pos', 'ref', 'alt'])
        for chr_val, sub_df in df_summary_stats.groupby('chr')
    }

    # Precompile regex pattern to extract chromosome from filename
    pattern = re.compile(r'\.MAF_threshold=0\.005\.(\d+)_combined\.csv')

    # List CSV files (skip targets.csv)
    csv_files = [f for f in os.listdir(directory_1000_genomes) if f.endswith('.csv') and f != 'targets.csv']
    
    # Use parallel processing to process CSV files concurrently
    dataframes = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        func = partial(process_csv_file,
                       directory_1000_genomes=directory_1000_genomes,
                       sad_columns=sad_columns,
                       pattern=pattern,
                       summary_stats_by_chr=summary_stats_by_chr)
        results = list(executor.map(func, csv_files))
    
    # Filter out None results and concatenate DataFrames
    for res in results:
        if res is not None and not res.empty:
            dataframes.append(res)
            
    if dataframes:
        df_summary_stats_result = pd.concat(dataframes, ignore_index=True)
    else:
        df_summary_stats_result = pd.DataFrame()

    print("----Finished", flush=True)

    ## STEP 2: Assign coding region value 
    print("Assigning coding regions", flush=True)
    coding_region_set = list(pd.read_csv(coding_snp_list_path)['snp'])
    df_summary_stats_result['in_coding_region'] = df_summary_stats_result['snp'].isin(coding_region_set)
   
    ## STEP 3: Compute p-values using t-test
    print("Computing pvals", flush=True)
    track_data = df_summary_stats_result[sad_columns].to_numpy()
    _, p_values = ttest_1samp(track_data, popmean=0, axis=1, nan_policy='omit')
    df_summary_stats_result['t_test_p_value'] = p_values

    ## STEP 4: Adjust p-values for multiple testing (FDR-BH)
    print("Adjusting for bonferroni-holm", flush=True)
    all_p_values = df_summary_stats_result['t_test_p_value'].to_numpy()
    n_snps_prev = len(df_summary_stats_result)
    fdr_significant_mask, adjusted_p_values, _, _ = multipletests(all_p_values, alpha=ALPHA, method='fdr_bh')
    df_summary_stats_result['adjusted_t_test_p_value'] = adjusted_p_values
    df_summary_stats_result['fdr_significant'] = fdr_significant_mask  

    # Create a reference list for comparison (p_value < 5e-8)
    df_summary_stats_signifcant_list = df_summary_stats_result[df_summary_stats_result['p_value'] < 5e-8].reset_index(drop=True)

    ## STEP 5: Compute the significant SNP list 
    print("Computing sigsnp list", flush=True)
    fdr_significant_snps = df_summary_stats_result[df_summary_stats_result['fdr_significant'] | df_summary_stats_result['in_coding_region']]
    n_snps_curr = len(fdr_significant_snps)

    ## STEP 6: Compute the new p-value threshold for GWAS selection 
    p_value_threshold = 5e-8 * (n_snps_prev / n_snps_curr)

    ## STEP 7: Filter any SNPs that do not meet the new p-value threshold
    print("Filtering snps", flush=True)
    filtered_df_summary_stats_result = df_summary_stats_result[df_summary_stats_result['p_value'] < p_value_threshold]

    ## STEP 8: Identify the leading SNPs by eliminating SNPs within WINDOW_SIZE of a leading SNP
    print("Identifying leading snps", flush=True)
    result_df = window_elimination(filtered_df_summary_stats_result, WINDOW_SIZE)
    result_df['left_border'] = result_df['pos'] - WINDOW_SIZE
    result_df['right_border'] = result_df['pos'] + WINDOW_SIZE

    # Process the reference list similarly
    if not df_summary_stats_signifcant_list.empty:
        df_summary_stats_signifcant_list = window_elimination(df_summary_stats_signifcant_list, WINDOW_SIZE)
        df_summary_stats_signifcant_list['left_border'] = df_summary_stats_signifcant_list['pos'] - WINDOW_SIZE
        df_summary_stats_signifcant_list['right_border'] = df_summary_stats_signifcant_list['pos'] + WINDOW_SIZE
    
    ## STEP 9: Compute metadata for downstream analyses
    print("Computing metadata", flush=True)
    num_snps_found = len(result_df)
    num_coding_snps_found = len(result_df[result_df['in_coding_region']])
    num_overlapping_snps = len(pd.merge(result_df, df_summary_stats_signifcant_list, left_on='snp', right_on='snp'))
    num_overlapping_loci = overlaps(df_summary_stats_signifcant_list, result_df)
    num_original_coding_snps = len(df_summary_stats_signifcant_list[df_summary_stats_signifcant_list["in_coding_region"]])

    metadata_dict = {
        "num_snps_found": num_snps_found,
        "num_coding_snps_found": num_coding_snps_found,
        "num_overlapping_snps": num_overlapping_snps,
        "num_overlapping_loci": num_overlapping_loci,
        "num_original_list": len(df_summary_stats_signifcant_list),
        "num_original_coding_snps": num_original_coding_snps,
        "p_value_threshold": p_value_threshold,
        "percentage_loci_recovered": num_overlapping_loci / len(df_summary_stats_signifcant_list) if len(df_summary_stats_signifcant_list) != 0 else 1
    }
    
    metadata_df = pd.DataFrame([metadata_dict])
    result_df = result_df[['chr', 'pos', 'ref', 'alt', 'snp', 'p_value', 'in_coding_region']]
    df_summary_stats_signifcant_list = df_summary_stats_signifcant_list[['chr', 'pos', 'ref', 'alt', 'snp', 'p_value', 'in_coding_region']]
    return result_df, metadata_df, df_summary_stats_signifcant_list

if __name__ == "__main__":
    """
    Inputs to main: 
        
        file_path_summary_stats: path to summary statistics file (.tsv or .tsv.bgz)
        directory_1000_genomes: directory housing the SAD score CSV files for all SNPs
        coding_snp_list_path: path to coding SNP dataframe (must have a 'snp' column)
        track_list: (optional) CSV file (no header) with one column of track numbers. If not provided, a default list is used.
    
    Outputs: 
        result_df: final list of SNPs (without SAD values or FDR t-test info)
        metadata_df: metadata about the filtering process
    
    How to run:
        python gwas_method_fdr_bonf_method.py <file_path_summary_stats> <directory_1000_genomes> <coding_snp_list_path> <track_list [optional]>
    """
    file_path_summary_stats = sys.argv[1]
    directory_1000_genomes = sys.argv[2]
    coding_snp_list_path = sys.argv[3]
    track_list = sys.argv[4] if len(sys.argv) > 4 else None

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
    
    result_df, metadata_df, df_summary_stats_signifcant_list = main(df_summary_stats, directory_1000_genomes, track_list, coding_snp_list_path)
    
    # Save results in deepcast_phenotypes folder under a subfolder for the phenotype.
    # We use the base name of the summary stats file (without extension) as the phenotype name.
    phenotype_name = os.path.splitext(os.path.basename(file_path_summary_stats))[0]
    output_folder = os.path.join("deepcast_phenotypes", phenotype_name)
    os.makedirs(output_folder, exist_ok=True)
    
    result_df.to_csv(os.path.join(output_folder, "result_df.csv"), index=False)
    metadata_df.to_csv(os.path.join(output_folder, "metadata_df.csv"), index=False)
    df_summary_stats_signifcant_list.to_csv(os.path.join(output_folder, "df_summary_stats_signifcant_list.csv"), index=False)
    
    print(f"Results saved in folder: {output_folder}")