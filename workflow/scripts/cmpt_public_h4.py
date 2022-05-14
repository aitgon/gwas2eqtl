from eqtl2gwas.PathManager import PathManager
from eqtl2gwas.PyTabix import bgzip, tabix_index
from eqtl2gwas.constants import coloc_raw_tsv_path, h4_cutoff

import os
import pandas
import pathlib
import sys


#%%
help_cmd_str = "todo"
try:
    coloc_all_tsv_path = sys.argv[1]
    coloc_h4_tsv_path = sys.argv[2]
    if len(sys.argv) > 3:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

outdir_path = os.path.join(PathManager.get_outdir_path(), os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
coloc_all_df = pandas.read_csv(coloc_all_tsv_path, sep="\t")

# %%
coloc_h4_df = coloc_all_df.loc[coloc_all_df['PP.H4.abf'] >= h4_cutoff, ]
# coloc_h4_tsv_path = os.path.join(outdir_path, "coloc_h4_all.tsv")
coloc_h4_df.to_csv(coloc_h4_tsv_path, sep="\t", index=False)

#%%
bgzip(coloc_h4_tsv_path)
tabix_index(coloc_h4_tsv_path + ".gz", chrom=1, start=2, end=2, skip=1)
coloc_h4_df.to_csv(coloc_h4_tsv_path, index=False, sep="\t")
