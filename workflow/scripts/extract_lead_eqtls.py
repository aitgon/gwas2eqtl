import os
import pathlib
import sys

import pandas
from statsmodels.stats import multitest

#%%
permuted_tsv_path = "/home/gonzalez/Software/process/ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/Kasela_2017/microarray/Kasela_2017_microarray_T-cell_CD8.permuted.tsv"
out_tsv_path = "/home/gonzalez/Software/process/fdr0.05/ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/Kasela_2017/microarray/Kasela_2017_microarray_T-cell_CD8.leadpair.tsv"
chrom_start_end_tsv_path = "/home/gonzalez/Software/process/fdr0.05/ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/Kasela_2017/microarray/Kasela_2017_microarray_T-cell_CD8.leadpair.region.tsv"
region = "genome"
fdr = 0.05

#%%
permuted_tsv_path = sys.argv[1]
out_tsv_path = sys.argv[2]
chrom_start_end_tsv_path = sys.argv[3]
region = sys.argv[4]  # genome or region in the format 1:10-100
fdr = float(sys.argv[5])  # fdr 0.05

#%%
df = pandas.read_csv(permuted_tsv_path, sep="\t", header=0)
pathlib.Path(os.path.dirname(out_tsv_path)).mkdir(parents=True, exist_ok=True)
df['start'] = df['position'] - 500000
df['end'] = df['position'] + 500000

#%%
if region != "genome":  # if not genome, then region in the format 1:100-10000
    chrom = region.split(':')[0]
    start = int(region.split(':')[1].split("-")[0])
    end = int(region.split(':')[1].split("-")[1])

    df = df.loc[(df['chromosome'] == chrom) & (df['position'] >= start) & (df['position'] <= end)]

#%%
rejected, pvalue_corrected = multitest.fdrcorrection(df['p_beta'].tolist(), alpha=fdr)
df['p_fdr'] = pvalue_corrected
df = df.loc[rejected]

#%%
columns = ['chromosome', 'position', 'molecular_trait_object_id', 'variant', 'beta', 'pvalue', 'molecular_trait_id', 'n_traits', 'n_variants', 'p_perm', 'p_beta', 'p_fdr']
df = df[columns]
df = df.sort_values(by=['chromosome', 'position'])

#%%
df.to_csv(out_tsv_path, sep="\t", index=False)

#%%
df['start'] = df['position'] - 500000
df['end'] = df['position'] + 500000
columns = ['chromosome', 'start', 'end']
df = df[columns].drop_duplicates()
df.to_csv(chrom_start_end_tsv_path, sep="\t", index=False, header=False)

