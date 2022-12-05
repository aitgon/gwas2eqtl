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
    eqtl_tsv_path = sys.argv[3]
    coloc_path_strf = sys.argv[4]
    if len(sys.argv) > 5:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

gwas_df = pandas.read_excel(gwas_ods_path, header=0)
gwas_identifier_lst = sorted(gwas_df["id"].tolist())

eqtl_df = pandas.read_csv(eqtl_tsv_path, sep="\t")
eqtl_df = eqtl_df.loc[eqtl_df['ftp_path'].str.contains('/ge/|/microarray/', regex=True, na=False), ]
eqtl_identifier_lst = (eqtl_df['ftp_path'].str.replace('.all.tsv.gz', '', regex=True)).str.split('/', expand=True)[10].tolist()

# Create all tables
engine = create_engine(url)
Base.metadata.create_all(engine)

gwas_eqtl_counter = 0
gwas_counter = 0
for gwas_id in gwas_identifier_lst:
    if gwas_counter % 10 == 0:
        print("GWAS counter: " + str(gwas_counter))
    for eqtl_id in eqtl_identifier_lst:
        coloc_tsv_path = coloc_path_strf.format(**{'gwas_id': gwas_id, 'eqtl_id': eqtl_id})
        if os.path.isfile(coloc_tsv_path):
            df = pandas.read_csv(coloc_tsv_path, sep="\t")
            df.columns = [c.lower() for c in df.columns]  # change to lower case
            df.columns = [c.replace('.', '_') for c in df.columns]  # replace dots with underline
            df.sort_values(by='pp_h4_abf', inplace=True, ascending=False)  # keep coloc with highest pp_h4_abf
            df.drop_duplicates(['chrom', 'pos', 'alt', 'eqtl_gene_id', 'gwas_id', 'eqtl_id'], inplace=True)
            df['rsid'] = df['rsid'].str.replace('rs', '').astype(int)
            df_index = df['chrom'].astype(str) + "_" + df['pos'].astype(str) + "_" + df['eqtl_gene_id'] + "_" + df['gwas_id'] + "_" + df['eqtl_id']
            df.set_index(df_index, inplace=True, verify_integrity=True)
            df.index.rename('id', inplace=True)
            # Delete and insert coloc data
            coloc_tab = Base.metadata.tables['coloc']
            stmt = coloc_tab.delete().where(coloc_tab.c.gwas_id==gwas_id).where(coloc_tab.c.eqtl_id==eqtl_id)
            engine.execute(stmt)
            df = df[[c.key for c in coloc_tab.columns][1:]]
            df.to_sql('coloc', con=engine, if_exists='append', index=True, index_label='id')
