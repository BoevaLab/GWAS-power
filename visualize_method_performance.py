#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path
import numpy as np
from tqdm import tqdm

# Set style
plt.style.use('seaborn-v0_8')  # Updated style name
sns.set_theme(style="whitegrid")
sns.set_palette("Set2")

def load_metadata():
    """Load metadata from all phenotype directories."""
    metadata_list = []
    phenotype_dirs = [d for d in Path("patho_phenotypes/deepcast_phenotypes").iterdir() if d.is_dir() and d.name.startswith("phecode-")]
    
    if not phenotype_dirs:
        raise ValueError("No phecode directories found in 'patho_phenotypes/deepcast_phenotypes' folder")
    
    print(f"Found {len(phenotype_dirs)} phecode directories")
    
    for pheno_dir in tqdm(phenotype_dirs, desc="Loading metadata"):
        metadata_file = pheno_dir / "metadata_df.csv"
        if metadata_file.exists():
            try:
                df = pd.read_csv(metadata_file)
                df['phenotype'] = pheno_dir.name
                metadata_list.append(df)
            except Exception as e:
                print(f"Error loading {metadata_file}: {e}")
        else:
            print(f"No metadata file found in {pheno_dir}")
    
    if not metadata_list:
        raise ValueError("No valid metadata files found in any phecode directory")
    
    print(f"Successfully loaded metadata from {len(metadata_list)} directories")
    return pd.concat(metadata_list, ignore_index=True)

def create_performance_plots(metadata_df):
    """Create various performance comparison plots."""
    # Create output directory if it doesn't exist
    os.makedirs('visualization_output', exist_ok=True)
    
    # 1. Overall SNP Discovery Comparison
    plt.figure(figsize=(12, 6))
    plt.scatter(metadata_df['num_original_list'], metadata_df['num_snps_found'], alpha=0.5)
    plt.plot([0, metadata_df['num_original_list'].max()], [0, metadata_df['num_original_list'].max()], 
             'r--', label='Equal Discovery')
    plt.xlabel('Number of SNPs in Original GWAS')
    plt.ylabel('Number of SNPs Found by Our Method')
    plt.title('SNP Discovery Comparison')
    plt.legend()
    plt.tight_layout()
    plt.savefig('visualization_output/snps_discovery_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

    # 2. Coding SNP Discovery Comparison
    plt.figure(figsize=(12, 6))
    plt.scatter(metadata_df['num_original_coding_snps'], metadata_df['num_coding_snps_found'], alpha=0.5)
    plt.plot([0, metadata_df['num_original_coding_snps'].max()], 
             [0, metadata_df['num_original_coding_snps'].max()], 
             'r--', label='Equal Discovery')
    plt.xlabel('Number of Coding SNPs in Original GWAS')
    plt.ylabel('Number of Coding SNPs Found by Our Method')
    plt.title('Coding SNP Discovery Comparison')
    plt.legend()
    plt.tight_layout()
    plt.savefig('visualization_output/coding_snps_discovery_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

    # 3. Overlap Analysis
    plt.figure(figsize=(12, 6))
    overlap_ratio = metadata_df['num_overlapping_snps'] / metadata_df['num_original_list']
    sns.histplot(data=overlap_ratio, bins=30, alpha=0.7)
    plt.xlabel('Proportion of Original SNPs Recovered')
    plt.ylabel('Number of Phenotypes')
    plt.title('Distribution of SNP Recovery Rate')
    plt.tight_layout()
    plt.savefig('visualization_output/snp_recovery_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()

    # 4. Loci Recovery
    plt.figure(figsize=(12, 6))
    plt.scatter(metadata_df['num_original_list'], metadata_df['num_overlapping_loci'], alpha=0.5)
    plt.plot([0, metadata_df['num_original_list'].max()], [0, metadata_df['num_original_list'].max()], 
             'r--', label='Perfect Recovery')
    plt.xlabel('Number of Loci in Original GWAS')
    plt.ylabel('Number of Loci Recovered by Our Method')
    plt.title('Loci Recovery Comparison')
    plt.legend()
    plt.tight_layout()
    plt.savefig('visualization_output/loci_recovery_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

    # 5. Additional Discovery Analysis
    plt.figure(figsize=(12, 6))
    additional_snps = metadata_df['num_snps_found'] - metadata_df['num_overlapping_snps']
    additional_coding = metadata_df['num_coding_snps_found'] - metadata_df['num_original_coding_snps']
    
    plt.scatter(metadata_df['num_original_list'], additional_snps, alpha=0.5, label='All SNPs')
    plt.scatter(metadata_df['num_original_list'], additional_coding, alpha=0.5, label='Coding SNPs')
    plt.xlabel('Number of SNPs in Original GWAS')
    plt.ylabel('Number of Additional SNPs Found')
    plt.title('Additional SNP Discovery')
    plt.legend()
    plt.tight_layout()
    plt.savefig('visualization_output/additional_snps_discovery.png', dpi=300, bbox_inches='tight')
    plt.close()

    # 6. Relative Increase Analysis
    plt.figure(figsize=(15, 10))
    
    # Calculate relative increases
    relative_increase_all = (metadata_df['num_snps_found'] - metadata_df['num_original_list']) / metadata_df['num_original_list']
    relative_increase_coding = (metadata_df['num_coding_snps_found'] - metadata_df['num_original_coding_snps']) / metadata_df['num_original_coding_snps']
    
    # Create subplot 1: Relative increase vs original SNPs
    plt.subplot(2, 2, 1)
    plt.scatter(metadata_df['num_original_list'], relative_increase_all, alpha=0.5)
    plt.xlabel('Number of SNPs in Original GWAS')
    plt.ylabel('Relative Increase in SNPs')
    plt.title('Relative Increase vs Original SNPs')
    
    # Create subplot 2: Relative increase distribution
    plt.subplot(2, 2, 2)
    sns.histplot(data=relative_increase_all, bins=30, alpha=0.7)
    plt.xlabel('Relative Increase in SNPs')
    plt.ylabel('Number of Phenotypes')
    plt.title('Distribution of Relative Increases')
    
    # Create subplot 3: Coding vs All SNPs relative increase
    plt.subplot(2, 2, 3)
    plt.scatter(relative_increase_all, relative_increase_coding, alpha=0.5)
    plt.plot([0, max(relative_increase_all.max(), relative_increase_coding.max())], 
             [0, max(relative_increase_all.max(), relative_increase_coding.max())], 
             'r--', label='Equal Increase')
    plt.xlabel('Relative Increase in All SNPs')
    plt.ylabel('Relative Increase in Coding SNPs')
    plt.title('Coding vs All SNPs Relative Increase')
    plt.legend()
    
    # Create subplot 4: Top phenotypes by relative increase
    plt.subplot(2, 2, 4)
    top_phenotypes = metadata_df.nlargest(10, 'num_snps_found')
    plt.barh(range(len(top_phenotypes)), 
             (top_phenotypes['num_snps_found'] - top_phenotypes['num_original_list']) / top_phenotypes['num_original_list'])
    plt.yticks(range(len(top_phenotypes)), top_phenotypes['phenotype'])
    plt.xlabel('Relative Increase in SNPs')
    plt.title('Top 10 Phenotypes by Relative Increase')
    
    plt.tight_layout()
    plt.savefig('visualization_output/relative_increase_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()

    # 7. Detailed SNP Recovery Analysis
    plt.figure(figsize=(15, 10))
    
    # Create subplot 1: Recovery rate vs original SNPs
    plt.subplot(2, 2, 1)
    recovery_rate = metadata_df['num_overlapping_snps'] / metadata_df['num_original_list']
    plt.scatter(metadata_df['num_original_list'], recovery_rate, alpha=0.5)
    plt.xlabel('Number of SNPs in Original GWAS')
    plt.ylabel('SNP Recovery Rate')
    plt.title('Recovery Rate vs Original SNPs')
    
    # Create subplot 2: Absolute numbers comparison
    plt.subplot(2, 2, 2)
    plt.scatter(metadata_df['num_original_list'], metadata_df['num_overlapping_snps'], 
               alpha=0.5, label='Recovered SNPs')
    plt.scatter(metadata_df['num_original_list'], metadata_df['num_snps_found'], 
               alpha=0.5, label='Total SNPs Found')
    plt.plot([0, metadata_df['num_original_list'].max()], [0, metadata_df['num_original_list'].max()], 
             'r--', label='Perfect Recovery')
    plt.xlabel('Number of SNPs in Original GWAS')
    plt.ylabel('Number of SNPs')
    plt.title('Absolute SNP Numbers')
    plt.legend()
    
    # Create subplot 3: Recovery rate distribution
    plt.subplot(2, 2, 3)
    sns.histplot(data=recovery_rate, bins=30, alpha=0.7)
    plt.xlabel('SNP Recovery Rate')
    plt.ylabel('Number of Phenotypes')
    plt.title('Distribution of Recovery Rates')
    
    # Create subplot 4: Additional vs Recovered SNPs
    plt.subplot(2, 2, 4)
    additional_snps = metadata_df['num_snps_found'] - metadata_df['num_overlapping_snps']
    plt.scatter(metadata_df['num_overlapping_snps'], additional_snps, alpha=0.5)
    plt.xlabel('Number of Recovered SNPs')
    plt.ylabel('Number of Additional SNPs')
    plt.title('Additional vs Recovered SNPs')
    
    plt.tight_layout()
    plt.savefig('visualization_output/detailed_snp_recovery_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()

    # 8. Performance Summary
    summary_stats = {
        'Average SNP Recovery Rate': metadata_df['num_overlapping_snps'].mean() / metadata_df['num_original_list'].mean(),
        'Average Loci Recovery Rate': metadata_df['num_overlapping_loci'].mean() / metadata_df['num_original_list'].mean(),
        'Average Additional SNPs Found': (metadata_df['num_snps_found'] - metadata_df['num_overlapping_snps']).mean(),
        'Average Additional Coding SNPs': (metadata_df['num_coding_snps_found'] - metadata_df['num_original_coding_snps']).mean(),
        'Number of Phenotypes': len(metadata_df),
        'Median SNP Recovery Rate': metadata_df['num_overlapping_snps'].median() / metadata_df['num_original_list'].median(),
        'Median Loci Recovery Rate': metadata_df['num_overlapping_loci'].median() / metadata_df['num_original_list'].median(),
        'Average Relative Increase': ((metadata_df['num_snps_found'] - metadata_df['num_original_list']) / metadata_df['num_original_list']).mean(),
        'Median Relative Increase': ((metadata_df['num_snps_found'] - metadata_df['num_original_list']) / metadata_df['num_original_list']).median()
    }
    
    # Save summary statistics
    with open('visualization_output/performance_summary.txt', 'w') as f:
        f.write("Performance Summary:\n")
        for key, value in summary_stats.items():
            f.write(f"{key}: {value:.3f}\n")

def create_phenotype_specific_plots(metadata_df):
    """Create plots for specific phenotype categories."""
    # Group phenotypes by type (continuous vs categorical)
    metadata_df['phenotype_type'] = metadata_df['phenotype'].apply(
        lambda x: 'Continuous' if x.startswith('continuous') else 'Categorical'
    )
    
    # Compare performance between continuous and categorical phenotypes
    plt.figure(figsize=(12, 6))
    sns.boxplot(x='phenotype_type', y='percentage_loci_recovered', data=metadata_df)
    plt.title('Loci Recovery Rate by Phenotype Type')
    plt.tight_layout()
    plt.savefig('visualization_output/phenotype_type_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Additional phenotype type comparisons
    metrics = ['num_snps_found', 'num_coding_snps_found', 'num_overlapping_snps', 'num_overlapping_loci']
    for metric in metrics:
        plt.figure(figsize=(12, 6))
        sns.boxplot(x='phenotype_type', y=metric, data=metadata_df)
        plt.title(f'{metric.replace("_", " ").title()} by Phenotype Type')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f'visualization_output/{metric}_by_phenotype_type.png', dpi=300, bbox_inches='tight')
        plt.close()

def main():
    # Load metadata
    print("Loading metadata from all phenotypes...")
    metadata_df = load_metadata()
    
    # Create general performance plots
    print("Creating performance comparison plots...")
    create_performance_plots(metadata_df)
    
    # Create phenotype-specific plots
    print("Creating phenotype-specific plots...")
    create_phenotype_specific_plots(metadata_df)
    
    print("All visualizations have been saved in the 'visualization_output' directory.")

if __name__ == "__main__":
    main() 