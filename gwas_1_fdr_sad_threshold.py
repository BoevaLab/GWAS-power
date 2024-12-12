import os
import glob
import pandas as pd
from scipy.stats import ttest_1samp
from statsmodels.stats.multitest import multipletests
import numpy as np
import warnings
warnings.filterwarnings("ignore")

# Author; Sophie Sigfstead
# Purpose: Begin trialling different statistical methods for thresholding SAD scores and p-value threshold

def overlaps (df1, df2):
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

def main(directory_gwas_combined_files, track_number_threshold=1):
    combined_files = glob.glob(os.path.join(directory_gwas_combined_files, '*.csv'))
    coding_region_set = set(pd.read_csv('./gwas_1_single_track_analysis/coding_region_set.csv', usecols=['snp'])['snp'])

    first_file_df = pd.read_csv(combined_files[0], usecols=['snp', 'p_value', 'chr', 'pos'] + [col for col in pd.read_csv(combined_files[0], nrows=0).columns if col.startswith('SAD')])
    first_file_df = first_file_df.dropna(subset=['snp', 'p_value'])
    first_file_df_copy = first_file_df.copy()
    first_file_df['in_coding_region'] = first_file_df['snp'].isin(coding_region_set)
    coding_region_list = first_file_df[first_file_df['in_coding_region'] & first_file_df['p_value'].notna()]['snp'].unique().tolist()

    num_chunks = 10
    all_snp_scores_chunks = [pd.DataFrame() for _ in range(num_chunks)]

    for file in combined_files:
        df = pd.read_csv(file, usecols=['snp', 'p_value'] + [col for col in pd.read_csv(file, nrows=0).columns if col.startswith('SAD')])
        df = df[df['snp'].notna() & ~df['snp'].isin(coding_region_set) & df['p_value'].notna()]
        track_col = [col for col in df.columns if col.startswith('SAD')][0]

        for i in range(num_chunks):
            rows = df.iloc[i::num_chunks]
            if all_snp_scores_chunks[i].empty:
                all_snp_scores_chunks[i] = rows[['snp', track_col]].copy()
            else:
                all_snp_scores_chunks[i] = pd.merge(all_snp_scores_chunks[i], rows[['snp', track_col]], on='snp', how='inner')

    alpha = 0.05
    
    print("\n*** ALPHA ***: ", alpha)

    significant_snps = []
    all_p_values = []
    all_snp_ids = []

    for chunk in all_snp_scores_chunks:
        if chunk.empty:
            continue
        track_data = chunk.iloc[:, 1:].apply(pd.to_numeric, errors='coerce').values
        valid_mask = ~np.isnan(track_data).all(axis=1)
        snp_ids = chunk.iloc[:, 0][valid_mask].values
        track_data = track_data[valid_mask]

        _, p_values = ttest_1samp(track_data, popmean=0, axis=1, nan_policy='omit')
        all_p_values.extend(p_values)
        all_snp_ids.extend(snp_ids)

    all_p_values = np.array(all_p_values)
    all_snp_ids = np.array(all_snp_ids)

    # Apply FDR Correction
    fdr_significant_mask, adjusted_p_values, _, _ = multipletests(all_p_values, alpha=alpha, method='fdr_bh')
    significant_snps = all_snp_ids[fdr_significant_mask]

    # Remove duplicates
    significant_snps = list(set(significant_snps))
    
    filtered_snp_list_df = first_file_df.query("snp in @significant_snps or snp in @coding_region_list")
    filtered_snp_list = filtered_snp_list_df['snp'].tolist()

    final_filtered_snp_list = first_file_df[first_file_df['snp'].isin(filtered_snp_list)]

    result_df = pd.DataFrame()
    for chr_num in range(1, 23):
        chr_df = final_filtered_snp_list[final_filtered_snp_list['chr'] == chr_num].sort_values(by='p_value')
        while not chr_df.empty:
            top_snp = chr_df.iloc[0]
            result_df = result_df._append(top_snp, ignore_index=True)
            chr_df = chr_df[~((chr_df['pos'] - top_snp['pos']).abs() < 1000000)].reset_index(drop=True)

    print("Significant SNPs after FDR correction:", len(result_df))
    
        # Convert significant_snps to DataFrame
    significant_snps_df = pd.DataFrame(significant_snps, columns=['snp'])
    significant_snps_df = pd.merge(significant_snps_df, first_file_df, on ='snp', how='left')

    mdd_sig_snps = pd.read_csv(f'./gwas_window_size_analysis/snp_lists_results/id=reference/window=1000000/filtered_snps_gwas_1_sd=0.0.csv')  # usually this is mdd_sig_snps.csv
    
    # Merge the two DataFrames
    matching_snps = pd.merge(significant_snps_df, mdd_sig_snps, left_on='snp', right_on='snp')
    print(f"OVERLAP SNPS, OUR METHOD BH vs ORIGINAL 29 SNPS: {len(matching_snps)}")

    mdd_sig_snps['left_border'] = mdd_sig_snps['pos'] - 1000000
    mdd_sig_snps['right_border'] = mdd_sig_snps['pos'] + 1000000
    significant_snps_df['left_border'] = significant_snps_df['pos']-1000000
    significant_snps_df['right_border'] = significant_snps_df['pos']+1000000
    
    # Overlap computation
    overlap_count_1 = overlaps(mdd_sig_snps, significant_snps_df)
    print("OVERLAP LOCI, OUR METHOD BH vs ORIGINAL 29 SNPS: ", overlap_count_1)

    # Create the GWAS 1 Bonferroni Holm List
    mdd_sig_snps_bonferroni_full = first_file_df_copy.copy()# this is simply a way to get a copy of the summary stats

    mdd_sig_snps_bonferroni=pd.DataFrame()
    # Implement LD mechanism +/- 1Mb around significant snp eliminated
    for chr_num in range(1, 23):
        chr_df = mdd_sig_snps_bonferroni_full[mdd_sig_snps_bonferroni_full['chr'] == chr_num].sort_values(by='p_value')
        while not chr_df.empty:
            top_snp = chr_df.iloc[0]
            mdd_sig_snps_bonferroni= mdd_sig_snps_bonferroni._append(top_snp, ignore_index=True)
            chr_df = chr_df[~((chr_df['pos'] - top_snp['pos']).abs() < 1000000)].reset_index(drop=True)


    mdd_sig_snps_bonferroni['left_border'] = mdd_sig_snps_bonferroni['pos']-1000000
    mdd_sig_snps_bonferroni['right_border'] = mdd_sig_snps_bonferroni['pos']+1000000
    #print(mdd_sig_snps_bonferroni.head())

    # Compute the Bonferroni Holm List for GWAS 1 original snp list 
    M = len(mdd_sig_snps_bonferroni) # number of tests
    mdd_sig_snps_bonferroni = mdd_sig_snps_bonferroni.sort_values(by='p_value') # sort by p-value
    mdd_sig_snps_bonferroni['k'] = range(1, M + 1) # create rank value
    mdd_sig_snps_bonferroni['adjusted_alpha'] = alpha / (M - mdd_sig_snps_bonferroni['k'] + 1)
    mdd_sig_snps_bonferroni['reject_null'] = mdd_sig_snps_bonferroni['p_value'] <= mdd_sig_snps_bonferroni['adjusted_alpha']
    mdd_sig_snps_bonferroni = mdd_sig_snps_bonferroni[mdd_sig_snps_bonferroni['reject_null']]
    #print(mdd_sig_snps_bonferroni.head())
    print("Number of snps in Bonferroni-Holm GWAS 1 ORIGINAL (ie. new reference list):", len(mdd_sig_snps_bonferroni) )

    # compute the intersection of the lists generated by running bonferroni holm on both GWAS 1 result sets 
    # (= significant_snps and the list of snps generated by runing Bonferroni Holm on the GWAS 1 summary stats)
    matching_bonf_snps = pd.merge(significant_snps_df, mdd_sig_snps_bonferroni, left_on=['snp'], right_on=['snp'])
    print(f"OVERLAP SNPS - OUR METHOD vs ORIGINAL, BH: {len(matching_bonf_snps)}")

    # overlap with the Bonferroni Holm GWAS 1 list
    overlap_count_2 = overlaps(mdd_sig_snps_bonferroni, significant_snps_df)

    print("OVERLAP LOCI - OUR METHOD vs ORIGINAL, BH: " ,overlap_count_2)
    print("Percentage of loci recovered with our method",  float(overlap_count_2/len(mdd_sig_snps_bonferroni)*100))
    return

if __name__ == "__main__":
    import sys
    directory_gwas_combined_files = sys.argv[1] # this should be the directory of the summary statistics - SAD score matched files
    main(directory_gwas_combined_files)
