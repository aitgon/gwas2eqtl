import os
import pathlib
import sys

import pandas
from statsmodels.stats import multitest

#%%
# permuted_tsv_path = "/home/gonzalez/Software/process/ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/Kasela_2017/microarray/Kasela_2017_microarray_T-cell_CD8.permuted.tsv"
permuted_tsv_path = sys.argv[1]
permuted_df = pandas.read_csv(permuted_tsv_path, sep="\t", header=0)
# out_tsv_path = "out.tsv"
out_tsv_path = sys.argv[2]
region = sys.argv[3]  # genome or region in the format 1:10-100
fdr = float(sys.argv[4])  # fdr 0.05

pathlib.Path(os.path.dirname(out_tsv_path)).mkdir(parents=True, exist_ok=True)

if region != "genome":  # if not genome, then region in the format 1:100-10000
    chrom = region.split(':')[0]
    start = int(region.split(':')[1].split("-")[0])
    end = int(region.split(':')[1].split("-")[1])

    permuted_df = permuted_df.loc[(permuted_df['chromosome'] == chrom) & (permuted_df['position'] >= start) & (permuted_df['position'] <= end)]

#%%
# import pdb; pdb.set_trace()
rejected, pvalue_corrected = multitest.fdrcorrection(permuted_df['p_beta'].tolist(), alpha=fdr)
permuted_df['p_fdr'] = pvalue_corrected

#%%
out_df = permuted_df.loc[rejected]
out_df.to_csv(out_tsv_path, sep="\t", index=False)

