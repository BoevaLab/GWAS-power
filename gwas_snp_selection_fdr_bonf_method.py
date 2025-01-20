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


# Author; Sophie Sigfstead
# Purpose: Latest version of our filtering method to improve GWAS power. 

def overlaps (df1, df2):
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
        
        # Check for overlap with each interval in df2
        for j, row2 in df2.iterrows():
            left2, right2 = row2['left_border'], row2['right_border']
            chr2  = row2['chr']
            
            # Check if the intervals overlap
            if (left1 <= right2 and right1 >= left2) and (chr1 == chr2):
                overlap_count += 1
                break  # Exit the loop once an overlap is found for the current interval in df1

    return overlap_count

def window_elimination(df, WINDOW_SIZE):
    """
    Input: Dataframe of summary statistics under the pvalue threshold
    Output: Dataframe containing only leading snps (snps within leading snp position +/- Window size are eliminated)
    """
    result_df = pd.DataFrame()
    for chr_num in range(1, 23):
        chr_df = df[df['chr'] == chr_num].sort_values(by='p_value')
        while not chr_df.empty:
            top_snp = chr_df.iloc[0]
            result_df = result_df._append(top_snp, ignore_index=True)
            chr_df = chr_df[~((chr_df['pos'] - top_snp['pos']).abs() < WINDOW_SIZE)].reset_index(drop=True)
    return result_df 

def main(df_summary_stats, directory_1000_genomes, track_list, coding_snp_list_path):
    """

    1. Process GWAS summary stats and merge with SAD track information from 1000 Genomes files.
    2. Perform filtering of GWAS Snps, returning final results list 
    
    Parameters:
        df_summary_stats (pd.DataFrame): Summary stats DataFrame.
        directory_1000_genomes (str): Directory containing 1000 Genomes CSV files.
        track_list (list): List of track identifiers for SAD columns.

    """
    ### CONSTANTS 
    # Significance level 
    ALPHA = 0.05

    # Window Size - This is the range of bp eliminated around a leading SNP. 
    # For example, it if the leading snp is at position 2000000, and window size is 1000000, this will define the window as 1000000 to 3000000.
    WINDOW_SIZE = 1000000

    ## STEP 1: Merge the summary statistics with the SAD value data

    # Reverse the -log(pvalue) operation 
    df_summary_stats['p_value'] = 10 ** (-df_summary_stats['neglog10_pval_EUR'])
    df_summary_stats.drop(columns=['neglog10_pval_EUR'], inplace=True)

    # Drop snps without p_value 
    df_summary_stats.dropna(subset=['p_value'], inplace=True)
    df_summary_stats.reset_index(drop=True, inplace=True)

    # Initialize SAD columns with NaN
    sad_columns = [f"SAD{track}" for track in track_list]
    
    # Create a list for results of processing each file in 1000genomes_as_csv (each corresponds to chromosome).
    # Each result in this list corresponds to the merged summary stats and SAD information form 1000genomes_as_csv for one chromosome.
    dataframes = []
    
    # Iterate through CSV files
    csv_files = [f for f in os.listdir(directory_1000_genomes) if f.endswith('.csv') and f != 'targets.csv']
    for csv_file in csv_files:
        # Load the csv from 1000 genomes
        csv_file_path = os.path.join(directory_1000_genomes, csv_file)
        chunk = pd.read_csv(csv_file_path, usecols=['chr', 'pos', 'ref', 'alt', 'snp'] + sad_columns)

        # Extract chromosome from csv filename
        match = re.search(r'\.MAF_threshold=0\.005\.(\d+)_combined\.csv', csv_file)
        if match:
            chromosome = int(match.group(1))

        # Filter summary stats for the current chromosome
        df_summary_stats_chr = df_summary_stats.query("chr == @chromosome")

        # Merge chunk from 1000genomes_as_csv with filtered summary stats
        merged_result = pd.merge(
            df_summary_stats_chr,
            chunk,
            on=['chr', 'pos', 'ref', 'alt'],
            how='left'
        )
        
        # Keep only rows with valid SAD and pval values
        merged_result = merged_result[merged_result[sad_columns[0]].notna() & merged_result['p_value'].notna()]

        # concatenate the result to the data frame list
        dataframes.append(merged_result)

    # Concatenate all results into one massive table - this now has all SNPs info (chr, pos, etc.) with SAD track values for the tracks in track_list
    df_summary_stats_result = pd.concat(dataframes, ignore_index=True)

    ## STEP 2: assign coding region value 

    coding_region_set = list(pd.read_csv(coding_snp_list_path)['snp'])
    df_summary_stats_result['in_coding_region'] = df_summary_stats_result['snp'].isin(coding_region_set)
   
    ## STEP 3: Compute p-values
    # Extract the track data from the DataFrame
    track_data = df_summary_stats_result[sad_columns].to_numpy()

    # Perform the t-test
    _, p_values = ttest_1samp(track_data, popmean=0, axis=1, nan_policy='omit')

    # Add the p-values as a new column to the DataFrame
    df_summary_stats_result['t_test_p_value'] = p_values

    ## STEP 4: Compute the adjusted p-values for FDR bonferonni holm 
    all_p_values = df_summary_stats_result['t_test_p_value'].to_numpy()

    # this is the number of snps before filtering with FDR-BH (for calculating the updated pvalue threshold later)
    n_snps_prev = len(df_summary_stats_result)

    fdr_significant_mask, adjusted_p_values, _, _ = multipletests(all_p_values, alpha=ALPHA, method='fdr_bh')

    # Add adjusted p-values and significance mask to the DataFrame
    df_summary_stats_result['adjusted_t_test_p_value'] = adjusted_p_values
    df_summary_stats_result['fdr_significant'] = fdr_significant_mask  

    # Create a reference list for comparsion - returns list of summary stats with p-value < 5e-8 (used in Step 8 + 9)
    df_summary_stats_signifcant_list = df_summary_stats_result[df_summary_stats_result['p_value']< 5e-8].reset_index(drop=True)

    ## STEP 5: Compute the significant snp list 
    fdr_significant_snps = df_summary_stats_result[df_summary_stats_result['fdr_significant'] | df_summary_stats_result['in_coding_region'] ]
    n_snps_curr = len(fdr_significant_snps)

    ## STEP 6: Compute the new p value for GWAS selection 
    p_value_threshold=5e-8*(n_snps_prev/n_snps_curr)

    ## STEP 7: Filter any snps that do not meet the pvalue
    filtered_df_summary_stats_result = df_summary_stats_result[df_summary_stats_result['p_value']<p_value_threshold]

    ## STEP 8: Identify the leading SNPs by removing any that are within 500MB with a lower p-value
    # Eliminate regions around the leading snps
    result_df = window_elimination(filtered_df_summary_stats_result, WINDOW_SIZE)
    # Add the loci left and right border for the result and reference snp lists
    result_df['left_border'] = result_df['pos']-WINDOW_SIZE
    result_df['right_border'] = result_df['pos']+WINDOW_SIZE

    # Perform the same two operations on the "reference list"
    df_summary_stats_signifcant_list = window_elimination(df_summary_stats_signifcant_list, WINDOW_SIZE)
    df_summary_stats_signifcant_list ['left_border'] = df_summary_stats_signifcant_list['pos']-WINDOW_SIZE
    df_summary_stats_signifcant_list ['right_border'] = df_summary_stats_signifcant_list ['pos']+WINDOW_SIZE
    
    ## STEP 9: Compute metadata - these can be aggregated in analysis later. 
    num_snps_found = len(result_df)
    num_coding_snps_found = len(result_df[result_df['in_coding_region']])
    num_overlapping_snps = len(pd.merge(result_df, df_summary_stats_signifcant_list , left_on='snp', right_on='snp'))
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
    "percentage_loci_recovered" : num_overlapping_loci/len(df_summary_stats_signifcant_list)
    }
    
    #Create a DataFrame with one row using the dictionary
    metadata_df = pd.DataFrame([metadata_dict])

    # Remove SAD values and fdr t-test info from result df to reduce size
    result_df =  result_df[['chr', 'pos', 'ref', 'alt', 'snp', 'p_value', 'in_coding_region']]

    return result_df, metadata_df


if __name__ == "__main__":
    """
    Inputs to main: 
        
        file_path_summary_stats: path to file for summary statistics for the particular trait, in tsv format
        
        directory_1000_genomes: this houses the SAD score information for all of the SNPs. It should be available in the local repository

        track_list: list of tracks relevant to the trait, formatted as [x,<space>y,<space>z,...] where x,y,z are the numbers of the tracks
        - if no track_list is specified then all of the tracks (0-683) will be used. 
        - you can also modify the default track list in the main function if there are tracks that are more relevant across all traits.
    
        coding_snp_list_path: path to coding snp dataframe (must have column called snp filled with snp ids). See the file that I provided. 
    
    Outputs: 
        result_df: result list of snps. format is the same as the summary stats file 
        It does not include the SAD values or fdr t-test pvalues. Remove line 203 to include this information. 
        
        metadata_df: file which contains summary results for that trait, such as number of coding snps, etc.  
        
        These are returned as dataframes, modify the file if you would like them to be saved as csv's 

    How to run:
        python gwas_method_fdr_bonf_method.py <file_path_summary_stats> <path to directory_1000_genomes> <path to coding snp list> <track_list [optional]>
        for example:  python gwas_method_fdr_bonf_method.py biomarkers-30620-both_sexes-irnt.tsv 1000genomes_as_csv coding_region_set_neale_lab.csv numbers_1_to_40.csv
    """
    file_path_summary_stats = sys.argv[1] # this should be the path for the summary statistics for the given trait. It will then be matched with SAD score matched files
    directory_1000_genomes = sys.argv[2]
    coding_snp_list_path = sys.argv[3]
    track_list = sys.argv[4] if len(sys.argv) > 4 else None # this should be a csv of of track numbers (no header)

    """
    Note that since the files are large, ensure they are located in a directory not affected by MAC SIP. 
    Note that this will expand the tsv files from NealeLab from 2.3 to 8.8 GB. 
    """
    # Determine file type and open accordingly by identifying the open function for that particular file type
    # Columns to retain from the file - modify neglog10_pval_EUR if a different ethnicity is preferred 
    selected_columns = ['chr', 'pos', 'ref', 'alt', 'neglog10_pval_EUR']

    # Summary Stats Input
    # open the file - must be .tsv.bgz or .tsv format
    file = gzip.open(file_path_summary_stats, "rt") if file_path_summary_stats.endswith(".tsv.bgz") else open(file_path_summary_stats, "r")

    # create a dataframe for the summary stats
    df_summary_stats = pd.concat(pd.read_csv(file, sep="\t", usecols=selected_columns, chunksize=100000),ignore_index=True)
    
    # Track List Input
    # If a track list is not provided - this was mainly for testing but a default list can also be inputted or modified for chatgpt
    if track_list is None:
        track_list = list(range(40))
    # If a track list is provided - provide in csv file with no header and one column of numbers corresponding to the tracks 
    else:
        track_list_file = pd.read_csv(track_list, header=None)  # Assumes no header
        track_list = track_list_file.iloc[:, 0].tolist() # This should be a list of numbers 
    
    main(df_summary_stats,  directory_1000_genomes, track_list, coding_snp_list_path)

