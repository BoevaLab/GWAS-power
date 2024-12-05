import os
import sys
import pandas as pd
from scipy.stats import ttest_1samp

def main(directory_gwas_combined_files, track_number_threshold=1):

    # Step 1: Perform FDR on SAD scores to determine which SNPs should be passed to Step 2 (GWAS - Replication)
    combined_files = [f for f in os.listdir(directory_gwas_combined_files) if f.endswith('.csv')]

    # A. Establish which SNPs are coding so that they can be passed directly to Step 2
    coding_region_set = pd.read_csv('./gwas_1_single_track_analysis/coding_region_set.csv')
    coding_region_set = set(coding_region_set['snp'])

    # load the first file so that you can get the coding region snps for the particular GWAS
    first_file_df = pd.read_csv(os.path.join(directory_gwas_combined_files, combined_files[0]))
    track_col_first = [col for col in first_file_df.columns if col.startswith('SAD')][0]

    # coding snps list made
    first_file_df = first_file_df[first_file_df['snp'].notna() & first_file_df['p_value'].notna()]
    first_file_df['in_coding_region'] = first_file_df['snp'].isin(coding_region_set)
    coding_region_list = list(set(first_file_df[(first_file_df['in_coding_region']) & (first_file_df['p_value'].notna())]['snp']))

    # now, the filtering process. Files are organized row-wise as SNPs, column is the value for the SAD track. 
    # Must first compute test-stats for all non-coding SNPs.

    # data frame to hold the results of all the SAD values for each track
    all_snp_scores = pd.DataFrame()
    all_snp_scores['snp'] = first_file_df['snp']

    for file in combined_files:
        # load the file
        df = pd.read_csv(os.path.join(directory_gwas_combined_files, file), low_memory=False)
        # get only non-coding snps
        df = df[df['snp'].notna() & ~df['snp'].isin(coding_region_set) & df['p_value'].notna()]
        track_col = [col for col in df.columns if col.startswith('SAD')][0]
        # the null hypothesis is that an aggregation of 33 SAD scores 
        all_snp_scores = pd.merge(all_snp_scores, df, on='snp', how='inner')
        all_snp_scores = all_snp_scores[['snp', track_col]]
    
    print("ALL SNPS SCORES")
    print(all_snp_scores.head())

    # set the significance level - there was only 1 test done here?
    alpha = 0.05
    snps_pvalues_list = []
    threshold_snps = []

    for snp in list(all_snp_scores['snp']):
        snp_row = all_snp_scores[all_snp_scores["snp"] == snp]
        sad_values = snp_row.iloc[:, 1:].values.flatten().tolist()
        # compute the test stat
        t_stat, p_value = ttest_1samp(sad_values, 0)
        snps_pvalues_list.append({'snp': snp, 'p_value': p_value})
        if p_value < alpha:
            threshold_snps.insert(snp)

    print("P-value SCORES")
    snps_pvalues = pd.DataFrame(snps_pvalues_list)
    print(snps_pvalues.head())

    print("Threshold SNPs")
    print(len(threshold_snps))

    filtered_snp_list = first_file_df[[first_file_df['snp'].isin(threshold_snps) | first_file_df['snp'].isin(coding_region_list)]  & first_file_df['p_value'].notna() & first_file_df[ track_col_first].notna() ]
    filtered_snp_list = filtered_snp_list[['snp', 'p_value']]
   
    # now implemented Bonferroni holm 
    M = len(filtered_snp_list)
    sorted_filtered_snp_list['reject_null'] = False
    sorted_filtered_snp_list = filtered_snp_list.sort(by=['p_value'], ascending = True)
    sorted_filtered_snp_list['k'] = sorted_filtered_snp_list.index + 1

    # Perform the Bonferroni-Holm test
    for i, row in sorted_filtered_snp_list.iterrows():
        # Compute the adjusted alpha for the current test
        adjusted_alpha = alpha / (M - row['k'] + 1)
        
        # Check if the p-value is less than or equal to the adjusted alpha
        if row['p_value'] <= adjusted_alpha:
            sorted_filtered_snp_list.at[i, 'reject_null'] = True
        else:
            # Once you fail to reject a null hypothesis, subsequent tests also fail
            break

    # Filter for significant SNPs
    significant_snps = sorted_filtered_snp_list[sorted_filtered_snp_list['reject_null']]

    # Output results
    print("Significant SNPs after Bonferroni-Holm correction:")
    print(significant_snps)

    return

if __name__ == "__main__":
    directory_gwas_combined_files = sys.argv[1]
    main(directory_gwas_combined_files)
