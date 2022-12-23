import math
import os

import sqlalchemy
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
if sqlalchemy.inspect(engine).has_table("tophits"):
    tophits_tab = Base.metadata.tables['tophits']
    tophits_tab.drop(engine)
Base.metadata.create_all(engine)

concat_lst = []
concat_df = pandas.DataFrame()
id = 0
for gwas_id in gwas_id_lst:
    tophits_tsv_path = hg38_tsv_path_strf.format(gwas_id=gwas_id)
    df = pandas.read_csv(tophits_tsv_path, sep="\t", header=0)
    df.drop_duplicates(inplace=True)
    df.index = df.index + id
    concat_df = pandas.concat([concat_df, df], axis=0, verify_integrity=True)
    if df.shape[0] > 0:
        id = max(df.index) + 1

concat_df.drop_duplicates(inplace=True)
concat_df.sort_values('pval', ascending=True, inplace=True)
concat_df.drop_duplicates(subset=['chrom', 'pos', 'ea', 'gwas_id'], inplace=True)
concat_df['rsid'] = concat_df['rsid'].str.split('rs', expand=True)[1]

# Delete data if exists and insert
stmt = Base.metadata.tables['tophits'].delete()
engine.execute(stmt)

concat_df.to_sql('tophits', con=engine, if_exists='append', index=True, index_label='id')
