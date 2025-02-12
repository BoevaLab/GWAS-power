import pandas as pd

snps = pd.read_csv("genome_assembly/snp_result.tsv", sep='\t', usecols=['snp_id'])
snps = snps['snp_id']='rs'+snps['snp_id']

m = len(snps.index)

df = pd.read_csv("genome_assembly/dummy_exon_regions_v2.csv", nrows=m)

df['snp']=snps['snp_id']

df.to_csv('genome_assembly/coding_regions_with_snp_dummy.csv')