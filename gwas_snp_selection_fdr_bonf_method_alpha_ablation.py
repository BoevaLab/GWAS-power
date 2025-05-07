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

def ld_based_clumping(df, plink_path, ld_reference_path, p_value_column='p_value', snp_column='snp', chr_column='chr', bp_column='pos', r2_threshold=0.2, kb_radius=500, p_threshold=1):
    """
    Perform LD-based clumping using PLINK.
    
    Parameters:
        df (pd.DataFrame): Dataframe containing SNP information and p-values
        plink_path (str): Path to PLINK executable or command name
        ld_reference_path (str): Path to LD reference panel (PLINK format)
        p_value_column (str): Column name for p-values
        snp_column (str): Column name for SNP IDs
        chr_column (str): Column name for chromosome
        bp_column (str): Column name for base position
        r2_threshold (float): LD r² threshold for clumping
        kb_radius (int): Distance in kb to look for LD
        p_threshold (float): Significance threshold for index SNPs
        
    Returns:
        pd.DataFrame: Clumped SNP list with borders based on actual LD relationships
    """
    print(f"Performing LD-based clumping (r²={r2_threshold}, window={kb_radius}kb)", flush=True)
    
    # Check if input DataFrame is empty or has no SNPs below threshold
    filtered_df = df[df[p_value_column] <= p_threshold].copy()
    if filtered_df.empty:
        print(f"Warning: No SNPs pass the p-value threshold of {p_threshold}", flush=True)
        return pd.DataFrame()
    
    print(f"Input has {len(filtered_df)} SNPs passing p-value threshold", flush=True)
    
    # Validate PLINK command
    try:
        # Check if PLINK is executable
        test_cmd = f"{plink_path} --version"
        result = subprocess.run(test_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            print(f"Error: PLINK command failed. Please check your PLINK installation.", flush=True)
            return pd.DataFrame()
    except Exception as e:
        print(f"Error checking PLINK: {e}", flush=True)
        return pd.DataFrame()
    
    # Create temporary directory for intermediate files
    temp_dir = tempfile.mkdtemp()
    
    try:
        # Prepare input file for PLINK
        assoc_file = os.path.join(temp_dir, "plink_input.assoc")
        clumped_output = os.path.join(temp_dir, "plink_output")
        
        # Format the association file for PLINK
        plink_input = filtered_df[[snp_column, chr_column, bp_column, p_value_column]].copy()
        plink_input.columns = ['SNP', 'CHR', 'BP', 'P']
        plink_input.to_csv(assoc_file, sep='\t', index=False)
        
        if not os.path.exists(assoc_file) or os.path.getsize(assoc_file) == 0:
            print(f"Error: Failed to create PLINK input file", flush=True)
            return pd.DataFrame()
        
        # Run PLINK clumping
        cmd_str = f"{plink_path} --bfile {ld_reference_path} --clump {assoc_file} " \
                 f"--clump-p1 {p_threshold} --clump-p2 {p_threshold} " \
                 f"--clump-r2 {r2_threshold} --clump-kb {kb_radius} " \
                 f"--clump-field P --clump-snp-field SNP --out {clumped_output}"
        
        try:
            subprocess.run(cmd_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        except Exception as e:
            print(f"Error executing PLINK: {e}", flush=True)
            return pd.DataFrame()
        
        # Process PLINK results
        clumped_file = f"{clumped_output}.clumped"
        if not os.path.exists(clumped_file) or os.path.getsize(clumped_file) == 0:
            print(f"Warning: PLINK clumping produced no results", flush=True)
            return pd.DataFrame()
        
        # Get lead SNPs from clumped file
        try:
            clumped_results = pd.read_csv(clumped_file, delim_whitespace=True)
            lead_snps = clumped_results['SNP'].tolist()
            if not lead_snps:
                print("Warning: No lead SNPs identified", flush=True)
                return pd.DataFrame()
                
            print(f"Found {len(lead_snps)} lead SNPs after clumping", flush=True)
        except Exception as e:
            print(f"Error reading clumped file: {e}", flush=True)
            return pd.DataFrame()
        
        # Create a DataFrame with just the lead SNPs from the original data
        result_df = df[df[snp_column].isin(lead_snps)].copy()
        if result_df.empty:
            print("Warning: Lead SNPs not found in original dataset. Check SNP IDs.", flush=True)
            return pd.DataFrame()
        
        # Set default borders based on kb_radius
        result_df['left_border'] = result_df[bp_column] - (kb_radius * 1000)
        result_df['right_border'] = result_df[bp_column] + (kb_radius * 1000)
        
        # Get precise LD borders for each lead SNP
        print("Getting LD-based borders for each lead SNP...", flush=True)
        ld_count = 0
        
        for i, row in result_df.iterrows():
            lead_snp = row[snp_column]
            chr_val = row[chr_column]
            pos_val = row[bp_column]
            
            # Run PLINK to get LD information
            tags_output = os.path.join(temp_dir, f"tags_{lead_snp}")
            tags_cmd = f"{plink_path} --bfile {ld_reference_path} --r2 --ld-snp {lead_snp} " \
                      f"--ld-window-kb {kb_radius} --ld-window 999999 --ld-window-r2 {r2_threshold} " \
                      f"--out {tags_output}"
            
            try:
                subprocess.run(tags_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                ld_file = f"{tags_output}.ld"
                
                if os.path.exists(ld_file) and os.path.getsize(ld_file) > 0:
                    # Process LD results
                    ld_data = pd.read_csv(ld_file, delim_whitespace=True)
                    if not ld_data.empty:
                        # Find min and max positions of SNPs in LD
                        min_pos = float('inf')
                        max_pos = 0
                        for _, ld_row in ld_data.iterrows():
                            if ld_row['R2'] >= r2_threshold:
                                bp = ld_row['BP_B']
                                min_pos = min(min_pos, bp)
                                max_pos = max(max_pos, bp)
                        
                        # Only update if we found LD SNPs
                        if min_pos < float('inf') and max_pos > 0:
                            result_df.at[i, 'left_border'] = min_pos
                            result_df.at[i, 'right_border'] = max_pos
                            ld_count += 1
            except Exception as e:
                print(f"Error processing LD for SNP {lead_snp}: {e}", flush=True)
                # Continue with the next SNP - we'll use default borders for this one
        
        print(f"Updated LD borders for {ld_count} of {len(result_df)} lead SNPs", flush=True)
        return result_df
            
    except Exception as e:
        print(f"Error during LD-based clumping: {e}", flush=True)
        return pd.DataFrame()
    finally:
        # Clean up temporary directory
        try:
            shutil.rmtree(temp_dir)
        except:
            pass  # Ignore cleanup errors

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

def main(df_summary_stats, directory_1000_genomes, track_list, coding_snp_list_path, alpha, ld_based=True, plink_path=None, ld_reference_path=None):
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
    print("====== GWAS SNP Selection with FDR-Bonferroni Method ======", flush=True)
    
    ### CONSTANTS 
    WINDOW_SIZE = 500000  # Range of bp eliminated around a leading SNP if window-based method is used
    R2_THRESHOLD = 0.2    # LD r² threshold for clumping
    KB_RADIUS = 500       # Distance in kb to look for LD (similar to WINDOW_SIZE but in kb)

    ## STEP 1: Process summary statistics
    print("1. Processing summary statistics", flush=True)
    # Reverse the -log(pvalue) operation and drop the original column
    df_summary_stats['p_value'] = 10 ** (-df_summary_stats['neglog10_pval_EUR'])
    df_summary_stats.drop(columns=['neglog10_pval_EUR'], inplace=True)
    df_summary_stats.dropna(subset=['p_value'], inplace=True)
    df_summary_stats.reset_index(drop=True, inplace=True)
    print(f"   - Processed {len(df_summary_stats)} SNPs from summary statistics", flush=True)

    ## STEP 2: Merge with SAD data
    print("2. Merging with SAD values", flush=True)
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
    print("3. Identifying coding SNPs", flush=True)
    coding_region_set = list(pd.read_csv(coding_snp_list_path)['snp'])
    df_summary_stats_result['in_coding_region'] = df_summary_stats_result['snp'].isin(coding_region_set)
    coding_snp_count = df_summary_stats_result['in_coding_region'].sum()
    print(f"   - Found {coding_snp_count} coding SNPs", flush=True)
   
    ## STEP 4: Calculate statistical significance
    print("4. Calculating statistical significance", flush=True)
    # T-test for SAD significance
    track_data = df_summary_stats_result[sad_columns].to_numpy()
    mean = np.mean(track_data)
    sd = np.std(track_data)

    # From the old code:
    # _, p_values = ttest_1samp(track_data, popmean=0, axis=1, nan_policy='omit')
    # df_summary_stats_result['t_test_p_value'] = p_values

    # Multiple testing correction
    # TODO: for code review: can I just drop this one?? isn't this the same as p_values?
    # all_p_values = df_summary_stats_result['t_test_p_value'].to_numpy()

    n_snps_prev = len(df_summary_stats_result)
    
    # fdr_significant_mask, adjusted_p_values, _, _ = multipletests(all_p_values, alpha=alpha, method='fdr_bh')
    # df_summary_stats_result['adjusted_t_test_p_value'] = adjusted_p_values
    # df_summary_stats_result['fdr_significant'] = fdr_significant_mask

    df_summary_stats_result['sad_relevant'] = ((track_data < mean - alpha * sd) + (track_data > mean + alpha * sd)).any(axis=1)
    
    # Create the reference list
    df_summary_stats_signifcant_list = df_summary_stats_result[df_summary_stats_result['p_value'] < 5e-8].reset_index(drop=True)
    original_snps_not_empty = not df_summary_stats_signifcant_list.empty
    
    if original_snps_not_empty:
        print(f"   - Reference list contains {len(df_summary_stats_signifcant_list)} SNPs at p < 5e-8", flush=True)
    else:
        print("   - No SNPs meet genome-wide significance (p < 5e-8) in reference list", flush=True)

    ## STEP 5: Compute significance threshold
    print("5. Computing adjusted significance threshold", flush=True)
    # Get the significant SNP list
    method_relevant_snps = df_summary_stats_result[df_summary_stats_result['sad_relevant'] | df_summary_stats_result['in_coding_region']]
    n_snps_curr = len(method_relevant_snps)
    
    # Calculate new p-value threshold
    p_value_threshold = 5e-8 * (n_snps_prev / n_snps_curr) if n_snps_curr > 0 else 5e-8
    print(f"   - Adjusted p-value threshold: {p_value_threshold:.2e}", flush=True)

    ## STEP 6: Identify lead SNPs
    print("6. Identifying lead SNPs", flush=True)
    
    # Check if we should use LD-based or window-based approach
    if ld_based and plink_path and ld_reference_path:
        # Verify LD reference files exist
        bim_file = f"{ld_reference_path}.bim"
        bed_file = f"{ld_reference_path}.bed"
        fam_file = f"{ld_reference_path}.fam"
        
        if not (os.path.exists(bim_file) and os.path.exists(bed_file) and os.path.exists(fam_file)):
            print("   - Error: LD reference files not found, falling back to window-based approach", flush=True)
            ld_based = False
    else:
        if ld_based:
            print("   - Error: PLINK configuration incomplete, falling back to window-based approach", flush=True)
        ld_based = False
    
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
    else:
        print(f"   - Using window-based approach with {WINDOW_SIZE}bp window", flush=True)
        # Filter by p-value first
        filtered_df = method_relevant_snps[method_relevant_snps['p_value'] < p_value_threshold]
        result_df = window_elimination(filtered_df, WINDOW_SIZE)
        
        if not result_df.empty:
            result_df['left_border'] = result_df['pos'] - WINDOW_SIZE
            result_df['right_border'] = result_df['pos'] + WINDOW_SIZE
    
    # Process the reference list if it's not empty
    if original_snps_not_empty:
        print("7. Processing reference SNP list", flush=True)
        if ld_based:
            df_summary_stats_signifcant_list = ld_based_clumping(
                df_summary_stats_signifcant_list,
                plink_path,
                ld_reference_path,
                p_threshold=5e-8,
                r2_threshold=R2_THRESHOLD,
                kb_radius=KB_RADIUS
            )
        else:
            df_summary_stats_signifcant_list = window_elimination(df_summary_stats_signifcant_list, WINDOW_SIZE)
            if not df_summary_stats_signifcant_list.empty:
                df_summary_stats_signifcant_list['left_border'] = df_summary_stats_signifcant_list['pos'] - WINDOW_SIZE
                df_summary_stats_signifcant_list['right_border'] = df_summary_stats_signifcant_list['pos'] + WINDOW_SIZE
    
    ## STEP 7: Compute metadata
    print("8. Computing result statistics", flush=True)
    
    # Initialize metadata variables
    num_snps_found = len(result_df)
    num_coding_snps_found = 0
    num_overlapping_snps = 0
    num_overlapping_loci = 0
    num_original_list = len(df_summary_stats_signifcant_list)
    num_original_coding_snps = 0
    percentage_loci_recovered = 1.0
    
    # Calculate statistics if results are not empty
    if not result_df.empty:
        if 'in_coding_region' in result_df.columns:
            num_coding_snps_found = result_df['in_coding_region'].sum()
        
        if not df_summary_stats_signifcant_list.empty:
            # Count overlapping SNPs
            overlap_df = pd.merge(result_df, df_summary_stats_signifcant_list, left_on='snp', right_on='snp')
            num_overlapping_snps = len(overlap_df)
            
            # Count overlapping loci
            num_overlapping_loci = overlaps(df_summary_stats_signifcant_list, result_df)
            
            # Count original coding SNPs
            if 'in_coding_region' in df_summary_stats_signifcant_list.columns:
                num_original_coding_snps = df_summary_stats_signifcant_list['in_coding_region'].sum()
            
            # Calculate percentage recovered
            if num_original_list > 0:
                percentage_loci_recovered = num_overlapping_loci / num_original_list
    
    # Create metadata dictionary
    metadata_dict = {
        "num_snps_found": num_snps_found,
        "num_coding_snps_found": num_coding_snps_found,
        "num_overlapping_snps": num_overlapping_snps,
        "num_overlapping_loci": num_overlapping_loci,
        "num_original_list": num_original_list,
        "num_original_coding_snps": num_original_coding_snps,
        "p_value_threshold": p_value_threshold,
        "percentage_loci_recovered": percentage_loci_recovered,
        "ld_based_clumping": ld_based,
        "r2_threshold": R2_THRESHOLD if ld_based else None,
        "kb_radius": KB_RADIUS if ld_based else None,
        "window_size": None if ld_based else WINDOW_SIZE
    }
    
    metadata_df = pd.DataFrame([metadata_dict])
    
    # Prepare final output dataframes
    cols_to_keep = ['chr', 'pos', 'ref', 'alt', 'snp', 'p_value', 'in_coding_region']
    
    if not result_df.empty:
        available_cols = [col for col in cols_to_keep if col in result_df.columns]
        result_df = result_df[available_cols]
    
    if not df_summary_stats_signifcant_list.empty:
        available_cols = [col for col in cols_to_keep if col in df_summary_stats_signifcant_list.columns]
        df_summary_stats_signifcant_list = df_summary_stats_signifcant_list[available_cols]
    
    print(f"   - Final results: {num_snps_found} SNPs identified", flush=True)
    if num_original_list > 0:
        print(f"   - Recovered {num_overlapping_loci}/{num_original_list} loci ({percentage_loci_recovered:.1%})", flush=True)
    
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
    alpha = float(sys.argv[7]) if len(sys.argv) > 7 else 0.05

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
        alpha=alpha,
        ld_based=ld_based,
        plink_path=plink_path,
        ld_reference_path=ld_reference_path
    )
    
    # Save results in deepcast_phenotypes folder under a subfolder for the phenotype.
    # We use the base name of the summary stats file (without extension) as the phenotype name.
    phenotype_name = os.path.splitext(os.path.basename(file_path_summary_stats))[0]
    output_folder = os.path.join(f"deepcast_alpha_ablation/results_{str(alpha)[0]+str(alpha)[2:]}", phenotype_name)
    os.makedirs(output_folder, exist_ok=True)
    
    result_df.to_csv(os.path.join(output_folder, "result_df.csv"), index=False)
    metadata_df.to_csv(os.path.join(output_folder, "metadata_df.csv"), index=False)
    df_summary_stats_signifcant_list.to_csv(os.path.join(output_folder, "df_summary_stats_signifcant_list.csv"), index=False)
    
    print(f"Results saved in folder: {output_folder}")