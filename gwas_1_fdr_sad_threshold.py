import os
import glob
import pandas as pd
from scipy.stats import ttest_1samp

def main(directory_gwas_combined_files, track_number_threshold=1):
    # Step 1: Perform FDR on SAD scores to determine which SNPs should be passed to Step 2 (GWAS - Replication)
    combined_files = glob.glob(os.path.join(directory_gwas_combined_files, '*.csv'))

    # A. Establish which SNPs are coding so that they can be passed directly to Step 2
    coding_region_set = set(pd.read_csv('./gwas_1_single_track_analysis/coding_region_set.csv', usecols=['snp'])['snp'])

    # Load the first file to get the coding region SNPs for the particular GWAS
    first_file_df = pd.read_csv(combined_files[0], usecols=['snp', 'p_value', 'chr', 'pos'] + [col for col in pd.read_csv(combined_files[0], nrows=0).columns if col.startswith('SAD')])
    first_file_df = first_file_df.dropna(subset=['snp', 'p_value'])
    track_col_first = [col for col in first_file_df.columns if col.startswith('SAD')][0]

    # Create coding SNPs list
    first_file_df['in_coding_region'] = first_file_df['snp'].isin(coding_region_set)
    coding_region_list = first_file_df[first_file_df['in_coding_region'] & first_file_df['p_value'].notna()]['snp'].unique().tolist()

    # Prepare to split rows across multiple DataFrames
    num_chunks = 10
    all_snp_scores_chunks = [pd.DataFrame() for _ in range(num_chunks)]

    for file in combined_files[0:10]:
        df = pd.read_csv(file, usecols=['snp', 'p_value'] + [col for col in pd.read_csv(file, nrows=0).columns if col.startswith('SAD')])
        df = df[df['snp'].notna() & ~df['snp'].isin(coding_region_set) & df['p_value'].notna()]
        track_col = [col for col in df.columns if col.startswith('SAD')][0]

        # Merge each chunk with its share of rows
        for i in range(num_chunks):
            rows = df.iloc[i::num_chunks]  # Take every `num_chunks` row starting from `i`
            if all_snp_scores_chunks[i].empty:
                all_snp_scores_chunks[i] = rows[['snp', track_col]].copy()
            else:
                all_snp_scores_chunks[i] = pd.merge(all_snp_scores_chunks[i], rows[['snp', track_col]], on='snp', how='inner')
    
    for frame in all_snp_scores_chunks:
        print(frame.head())

    # Process each chunk independently for p-value computation
    alpha = 0.05
    significant_snps = []

    def compute_ttest(row):
        # Ensure all values are numeric and drop NaNs
        numeric_values = pd.to_numeric(row[1:], errors='coerce').dropna()
        if len(numeric_values) == 0:  # If no valid numeric values, skip this row
            return None
        t_stat, p_value = ttest_1samp(numeric_values, 0)
        return p_value

    for chunk in all_snp_scores_chunks:
        if chunk.empty:
            continue
        chunk['p_value'] = chunk.apply(compute_ttest, axis=1)
        chunk = chunk.dropna(subset=['p_value'])  # Drop rows where t-test could not be performed
        significant_snps += chunk[chunk['p_value'] < alpha]['snp'].tolist()

    # Remove duplicates
    significant_snps = list(set(significant_snps))

    print("Threshold SNPs length:")
    print(len(significant_snps))

    filtered_snp_list_df = first_file_df.query("snp in @significant_snps or snp in @coding_region_list")
    filtered_snp_list = filtered_snp_list_df['snp'].tolist()

    # Implement the LD mechanism - for GWAS 1, this is just SNP location +/- 1Mb
    final_filtered_snp_list = first_file_df[first_file_df['snp'].isin(filtered_snp_list)]

    result_df = pd.DataFrame()

    for chr_num in range(1, 23):
        chr_df = final_filtered_snp_list[final_filtered_snp_list['chr'] == chr_num].sort_values(by='p_value')
        while not chr_df.empty:
            top_snp = chr_df.iloc[0]
            result_df = result_df.append(top_snp, ignore_index=True)
            chr_df = chr_df[~((chr_df['pos'] - top_snp['pos']).abs() < 1000000)].reset_index(drop=True)

    # Implement Bonferroni-Holm p-value threshold
    M = len(result_df)
    print("number of snps filtered:", M)
    sorted_filtered_snp_list = result_df.sort_values(by='p_value')
    sorted_filtered_snp_list['k'] = range(1, M + 1)
    sorted_filtered_snp_list['adjusted_alpha'] = alpha / (M - sorted_filtered_snp_list['k'] + 1)
    sorted_filtered_snp_list['reject_null'] = sorted_filtered_snp_list['p_value'] <= sorted_filtered_snp_list['adjusted_alpha']
    significant_snps = sorted_filtered_snp_list[sorted_filtered_snp_list['reject_null']]

    # Output results
    print("Significant SNPs after Bonferroni-Holm correction:")
    print(significant_snps)

    return

if __name__ == "__main__":
    import sys
    directory_gwas_combined_files = sys.argv[1]
    main(directory_gwas_combined_files)
