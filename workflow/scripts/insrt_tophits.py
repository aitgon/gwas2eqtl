import math
import os

from sqlalchemy import create_engine
from gwas2eqtl.db import Base

import pandas
import sys


#%%
help_cmd_str = "todo"
try:
    url = sys.argv[1]
    gwas_ods_path = sys.argv[2]
    hg38_tsv_path_strf = sys.argv[3]
    if len(sys.argv) > 4:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

gwas_df = pandas.read_excel(gwas_ods_path, header=0)
gwas_id_lst = sorted(gwas_df["id"].tolist())

# Create all tables
engine = create_engine(url)
Base.metadata.create_all(engine)

concat_lst = []

for gwas_id in gwas_id_lst:
    tophits_tsv_path = hg38_tsv_path_strf.format(gwas_id=gwas_id)
    if os.path.isfile(tophits_tsv_path):
        df = pandas.read_csv(tophits_tsv_path, sep="\t", header=0)
        concat_lst = concat_lst + df.values.tolist()

concat_df = pandas.DataFrame(concat_lst, columns=df.columns)
concat_df = concat_df.drop_duplicates(keep='first')

# aggregate different pvals and betas
concat_df = concat_df.groupby(['chrom', 'pos', 'rsid', 'nea', 'ea', 'n', 'se', 'gwas_id', 'eaf', 'pos19']).agg({'beta': lambda x: ','.join([str(it) for it in x]), 'pval': lambda x: ','.join([str(it) for it in x])}).reset_index()
concat_df2_index = concat_df['chrom'].astype(str) + "_" + concat_df['pos'].astype(str) + "_" + concat_df['gwas_id']
concat_df.set_index(concat_df2_index, verify_integrity=True, inplace=True)
concat_df['rsid'] = concat_df['rsid'].str.replace('rs', '').astype(int)

# import pdb; pdb.set_trace()

# Delete data if exists and insert
tophits_tab = Base.metadata.tables['tophits']
stmt = tophits_tab.delete()
engine.execute(stmt)
concat_df.to_sql('tophits', con=engine, if_exists='append', index=True, index_label='id')
