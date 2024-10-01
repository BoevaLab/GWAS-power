import sys 
import pandas as pd
import subprocess
# Goal: Replicate study: https://www.nature.com/articles/s41588-018-0311-9
# Title: Genome-wide meta-analysis identifies new loci and functional pathways influencing Alzheimer’s disease risk

def merge_loci(df, distance_threshold=250000):
    merged = []
    df = df.sort_values(by=['chr', 'left_border']).reset_index(drop=True)  # Sort by chromosome and left border
    
    while not df.empty:
        current = df.iloc[0]  # Take the first row
        df = df.iloc[1:]  # Remove the first row from the dataframe
        
        # Find overlapping or nearby rows within the threshold
        overlapping_rows = df[
            (df['chr'] == current['chr']) &
            (df['left_border'] <= current['right_border'] + distance_threshold) &
            (df['right_border'] >= current['left_border'] - distance_threshold)
        ]
        
        # If there are overlaps, merge them
        if not overlapping_rows.empty:
            # Update the current interval
            min_left_border = min(current['left_border'], overlapping_rows['left_border'].min())
            max_right_border = max(current['right_border'], overlapping_rows['right_border'].max())
            
            # Get the row with the lowest p-value
            min_p_row = overlapping_rows.loc[overlapping_rows['p'].idxmin()]
            
            # Create the merged interval
            merged_interval = {
                'pos': min_p_row['pos'],
                'chr': current['chr'],
                'p': min_p_row['p'],
                'snp': min_p_row['snp'],
                'left_border': min_left_border,
                'right_border': max_right_border
            }
            
            # Remove the overlapping rows from the dataframe
            df = df.drop(overlapping_rows.index).reset_index(drop=True)
        else:
            # If no overlap, the current interval is finalized
            merged_interval = current.to_dict()
        
        merged.append(merged_interval)
    
    return pd.DataFrame(merged)


def main(gwas_study_text):
    sum_stats = pd.read_csv(gwas_study_text, sep = '\t')
    significant = sum_stats[sum_stats['p']<5e-8]
    # look up the rs-id's using the study from the same group from 2019
    sum_stats_2019 = pd.read_csv('gwas_4_alz_sum_stats_2019_Kunkle', sep = '\t')
    significant_with_snp_ids = pd.merge(significant, pd.DataFrame(sum_stats_2019[['SNP', 'CHR', "BP"]]), left_on =[' chr','PosGRCh37' ], right_on=['CHR', "BP"], how='left')
    
    #significant_with_snp_ids that are nan look up in another dataframe called gwas_1_matching
    gwas_1_matching = pd.read_csv('gwas_1_matching/result_SAD0.csv')
    nan_snp_ids = significant_with_snp_ids[significant_with_snp_ids['SNP'].isna()] 
    nan_snp_ids = pd.merge(nan_snp_ids, gwas_1_matching[['chr', 'pos', 'snp']], left_on=[' chr', 'PosGRCh37'], right_on=['chr', 'pos'], how='left')
    significant_with_snp_ids.loc[significant_with_snp_ids['SNP'].isna(), 'SNP'] = nan_snp_ids['snp']
    nan_snp_ids = significant_with_snp_ids[significant_with_snp_ids['SNP'].isna()]
    nan_snp_ids.to_csv("gwas_4_alz_nan_snp_ids.csv")
    significant_with_snp_ids_non_nan = significant_with_snp_ids[significant_with_snp_ids['SNP'].notna()]

    plink_path = "../tools/plink2/plink2"
    test_file_path = './1KG_genotypes/all_phase3_bin_eur_no_duplicates'

    all_variants = pd.DataFrame(significant_with_snp_ids['SNP']).rename(columns={'SNP': 'ID'})
    all_variants.to_csv('./gwas_4_alz_intermediate_files/all_variants_list', sep = "\t", index=False)
    result_df = pd.DataFrame(columns = ['pos', 'chr', 'p', 'snp', 'left_border', 'right_border'])

    for _, row in significant_with_snp_ids_non_nan.iterrows():
        pos = row['PosGRCh37']
        chr = row[' chr']
        p = row["p"]
        snp = row["SNP"]
        start = pos - 1500000
        end = pos + 1500000


        plink_make_bed_command = [
        f'{plink_path}', 
        '--bfile', test_file_path,  
        '--chr', str(chr),  
        '--from-bp', str(start),  
        '--to-bp', str(end),  
        '--make-bed',  
        '--out', f'./gwas_4_alz_intermediate_files/chr_bed_files/chr_bed_file_chr={chr}_pos={pos}'
        ]
        subprocess.run(plink_make_bed_command, check=True)

        plink_ld_command = [
            f'{plink_path}', 
            "--bfile", f'./gwas_4_alz_intermediate_files/chr_bed_files/chr_bed_file_chr={chr}_pos={pos}',  
            "--ld-window", "99999",
            "--ld-window-kb", "3000",
            "--ld-window-r2", "0.6",
            "--r2-unphased",
            "--ld-snp-list",  f'{'./gwas_4_alz_intermediate_files/all_variants_list'}',
            "--out", f"./gwas_4_alz_intermediate_files/ld_output_files/ld_loci_output_chr={chr}_pos={pos}"
        ]
        subprocess.run(plink_ld_command)

        def define_locus(snp, chr, pos):
            ld_data = pd.read_csv(f'./gwas_4_alz_intermediate_files/ld_output_files/ld_loci_output_chr={chr}_pos={pos}.vcor', sep = '\t')
            correlated_snps = ld_data[(ld_data['ID_A'] == snp) & (ld_data['UNPHASED_R2'] > 0.6)]
            positions = correlated_snps['POS_B'].tolist()
            if ld_data[(ld_data['ID_A'] == snp)].empty:
                print("snp doesn't exist in ld data")
            if positions:
                locus_start = min(positions)
                locus_end = max(positions)
            else:
                locus_start = locus_end = pos
            return locus_start, locus_end

        # define locus based on ld information
        left_border, right_border = define_locus(snp, chr, pos)
        result_df=result_df._append({"chr":chr, "pos": pos, "snp": snp, "p": p, "left_border": left_border, "right_border": right_border}, ignore_index=True)
        print(result_df.head())
        print(len(result_df))
        result_df.to_csv("./gwas_4_alz_intermediate_files/result_df.csv")
        merged_result_df = merge_loci(result_df)
        merged_result_df.to_csv("./gwas_4_alz_intermediate_files/merged_result_df.csv")
        print("length of merged df " , len(merged_result_df))

    return

if __name__ == "__main__":
    gwas_study_text = sys.argv[1]
    main(gwas_study_text)