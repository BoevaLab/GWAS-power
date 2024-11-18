import pandas as pd
import os
import sys
from intervaltree import IntervalTree
import re

def threshold_calculator(df, sd, track):
    df_with_sad = df[(df[track].notna()) & (df[track].abs() < 1) & (df[track].abs() >= 0)]
  
    std_dev = df_with_sad[track].std()
    mean = df_with_sad[track].mean()

    upper_threshold = mean + (sd * std_dev)
    lower_threshold = mean - (sd * std_dev)
    print(upper_threshold, lower_threshold)

    return upper_threshold, lower_threshold

  
def main(directory_gwas_combined_files, directory_gwas_sig_snps, sd, window_size, summary_result_storage_directory, snp_storage_directory):
    sd = float(sd)
    window_size = float(window_size)
    
    combined_files = [f for f in os.listdir(directory_gwas_combined_files) if f.endswith('.csv')]
    
    threshold_set = set()
    non_threshold_set = set()

    exon_regions = pd.read_csv('../../GWAS_Data/exon_regions_v2.csv')

    # Read the first file to initialize the in_coding_region column
    first_file_df = pd.read_csv(os.path.join(directory_gwas_combined_files, combined_files[0]))
    '''
    # create interval trees for each exon in exon regions 
    trees = {}
    for _, row in exon_regions.iterrows():
        chr = row['chr']
        if chr not in trees:
            trees[chr] = IntervalTree()
            print(f'chromosome {chr} has been added')
        else:
            trees[chr][row['start']:row['end']+1] = True
    
    # check if a SNP is within any exon region
    def is_within_exon(row):
        chrom = row['chr']
        pos = row['pos']
        if bool(trees[chrom][pos]):
            return True
        return False
    '''
    coding_region_set = pd.read_csv('../gwas_1_single_track_analysis/coding_region_set.csv')
    coding_region_set = set(coding_region_set['snp'])
    first_file_df = first_file_df[first_file_df['snp'].notna()]
    # Apply the function to check for overlaps
    first_file_df['in_coding_region'] = first_file_df['snp'].isin(coding_region_set)
    print('in_coding_region assignment is complete')
    print(len(first_file_df))
    
    coding_region_set = set(first_file_df[first_file_df['in_coding_region'] & first_file_df['p_value'].notna()]['snp'])
    coding_region_set_list = list(coding_region_set)
    coding_region_set_df = pd.DataFrame(coding_region_set_list, columns=['snp'])
    print('length of coding region set')
    print(len(coding_region_set))

    for file in combined_files:
        print("processing:", file)
        df = pd.read_csv(os.path.join(directory_gwas_combined_files, file), low_memory = False)
        df = df[df['snp'].notna() & ~df['snp'].isin(coding_region_set) & df['p_value'].notna()]
        track_col = [col for col in df.columns if col.startswith('SAD')][0]  
        upper_threshold, lower_threshold = threshold_calculator(df, sd, track_col)
   
        above_threshold = df[(df[track_col] > upper_threshold) | (df[track_col] < lower_threshold)]
        below_threshold = df[(df[track_col] <= upper_threshold) & (df[track_col] >= lower_threshold)]
        
        threshold_set.update(above_threshold['snp'])
        non_threshold_set.update(below_threshold['snp'])
        
    print("sd:", sd)
    
    print('n snps total', len(first_file_df))

    threshold_df = first_file_df[first_file_df['snp'].isin(threshold_set)]

    non_threshold_df = first_file_df[first_file_df['snp'].isin(non_threshold_set)]

    all_df = first_file_df.copy()
    
    common_marker_names = set(threshold_df['snp']).intersection(set(non_threshold_df['snp']))

    if common_marker_names:

        non_threshold_df = non_threshold_df[~non_threshold_df['snp'].isin(common_marker_names)]


    print('n snps in sets')
    print(len(threshold_df) + len(non_threshold_df) + len(coding_region_set))

    
    # Modification: the n_snps_prev should be equal to the length of the threshold and non-threshold df together
    n_snps_prev = len(threshold_df) + len(non_threshold_df) + len(coding_region_set)
    print("n snps prev", n_snps_prev)

    # Add all sig SNPs in the coding_region_set to result_df
    coding_region_df = first_file_df[(first_file_df['in_coding_region'])]

    
    n_snps_filtered = len(threshold_df) + len(coding_region_df)
    print("n snps filtered", n_snps_filtered)
    
    new_pvalue = 5e-8 * (n_snps_prev / n_snps_filtered)
    print("New p-value threshold: ", new_pvalue)

    threshold_df = threshold_df[(threshold_df['p_value'].notna()) & (threshold_df['p_value'] < new_pvalue)]
    # Add all sig SNPs in the coding_region_set to result_df
    coding_region_df = first_file_df[(first_file_df['in_coding_region']) & (first_file_df['p_value'] < new_pvalue)]
    
    threshold_and_coding = pd.concat([threshold_df, coding_region_df], ignore_index=True)

    result_df = pd.DataFrame()
    
    for chr_num in range(1, 23):
        chr_df = threshold_and_coding[threshold_and_coding['chr'] == chr_num ].sort_values(by='p_value').reset_index(drop=True)
        while not chr_df.empty:
            top_snp = chr_df.iloc[0]
            result_df = pd.concat([result_df, pd.DataFrame([top_snp])], ignore_index=True)
            pos_sig_snp = top_snp['pos']
            chr_df = chr_df[~((chr_df['pos'] < pos_sig_snp + window_size) & (chr_df['pos'] > pos_sig_snp - window_size))]
            chr_df = chr_df.reset_index(drop=True)
    
    
    print('Sig snps coding:', len(coding_region_df))
    directory = os.path.dirname(snp_storage_directory)
    os.makedirs(directory, exist_ok=True)
    result_df.to_csv(f'{snp_storage_directory}', index=False)

    print('number of snps found', len(result_df))
   
    #GWAS #1
    mdd_sig_snps = pd.read_csv(f'{directory_gwas_sig_snps}') # usually this is mdd_sig_snps.csv
    matching_snps = pd.merge(result_df, mdd_sig_snps, left_on=['chr', 'pos'], right_on=['chromosome', 'base_pair_location'])

    result_df['right_border'] = result_df['pos'] + window_size
    result_df['left_border'] = result_df['pos'] - window_size
    snp_to_in_coding_region = dict(zip(all_df['snp'], all_df['in_coding_region']))

    # Map the in_coding_region values to result_df based on the snp column
    result_df['in_coding_region'] = result_df['snp'].map(snp_to_in_coding_region)

    mdd_sig_snps['right_border'] = mdd_sig_snps['base_pair_location'] + window_size
    mdd_sig_snps['left_border'] = mdd_sig_snps['base_pair_location'] - window_size
    # Create a dictionary with (chr, pos) tuples as keys and in_coding_region as values
    snp_to_in_coding_region = dict(zip(zip(all_df['chr'], all_df['pos']), all_df['in_coding_region']))
    

    # Create a tuple of (chromosome, base_pair_location) for each row in mdd_sig_snps and map to the dictionary
    mdd_sig_snps['in_coding_region'] = mdd_sig_snps.apply(lambda row: snp_to_in_coding_region.get((row['chromosome'], row['base_pair_location'])), axis=1)


    overlap_count = 0

    # Iterate through each interval in df1
    for i, row1 in mdd_sig_snps.iterrows():
        left1, right1 = row1['left_border'], row1['right_border']
        chr1  = row1['chromosome']
        
        # Check for overlap with each interval in df2
        for j, row2 in result_df.iterrows():
            left2, right2 = row2['left_border'], row2['right_border']
            chr2  = row2['chr']
            
            # Check if the intervals overlap
            if (left1 <= right2 and right1 >= left2) and (chr1 == chr2):
                overlap_count += 1
                break  # Exit the loop once an overlap is found for the current interval in df1
    
                
    print("overlapping loci: " ,overlap_count)
    print(" new pvalue" , new_pvalue)
    print("n snps filtered", n_snps_filtered)
    print("n snps prev", n_snps_prev)

    results_storage = pd.DataFrame(columns=['num_snps_found', 'num_snps_overlap', 'num_coding_snps', 'num_loci_overlap', 'p_value'])
    new_row = pd.DataFrame([{
    'num_snps_found': len(result_df),
    'num_snps_overlap': len(matching_snps),
    'num_coding_snps': len(result_df[result_df['in_coding_region'] == True]),
    'num_loci_overlap': overlap_count,
    'p_value': new_pvalue
    }])

    results_storage = pd.concat([results_storage, new_row], ignore_index=True)
    directory = os.path.dirname(summary_result_storage_directory)
    os.makedirs(directory, exist_ok=True)
    results_storage.to_csv(summary_result_storage_directory, index=False)
    return
    
if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python script.py <directory_gwas_combined_files> <directory_gwas_sig_snps> <sd> <result storage dir>")
    else:
        directory_gwas_combined_files = sys.argv[1]
        directory_gwas_sig_snps = sys.argv[2]
        sd = sys.argv[3]
        window_size = sys.argv[4]
        summary_result_storage_directory =sys.argv[5]
        snp_storage_directory = sys.argv[6]
        main(directory_gwas_combined_files, directory_gwas_sig_snps, sd, window_size, summary_result_storage_directory, snp_storage_directory)
    sys.exit()
