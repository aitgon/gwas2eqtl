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

coloc_tab = Base.metadata.tables['coloc']
stmt = coloc_tab.delete()
engine.execute(stmt)

gwas_eqtl_counter = 0
gwas_counter = 0
for gwas_id in gwas_identifier_lst:
    concat_df = pandas.DataFrame()
    if gwas_counter % 10 == 0:
        print("GWAS counter: " + str(gwas_counter))
    for eqtl_id in eqtl_identifier_lst:
        if gwas_eqtl_counter % 100 == 0:
            print("\tGWAS/eQTL counter: " + str(gwas_eqtl_counter))
        coloc_tsv_path = coloc_path_strf.format(**{'gwas_id': gwas_id, 'eqtl_id': eqtl_id})
        if os.path.isfile(coloc_tsv_path):
            coloc_df = pandas.read_csv(coloc_tsv_path, sep="\t")
            coloc_df.columns = [c.lower() for c in coloc_df.columns]  # change to lower case
            coloc_df.columns = [c.replace('.', '_') for c in coloc_df.columns]  # replace dots with underline
            coloc_df.sort_values(by='pp_h4_abf', inplace=True, ascending=False)  # keep coloc with highest pp_h4_abf
            coloc_df.drop_duplicates(['chrom', 'pos', 'alt', 'eqtl_gene_id', 'gwas_id', 'eqtl_id'], inplace=True)
            coloc_df['rsid'] = coloc_df['rsid'].str.replace('rs', '').astype(int)
            df_index = coloc_df['chrom'].astype(str) + "_" + coloc_df['pos'].astype(str) + "_" + coloc_df['eqtl_gene_id'] + "_" + coloc_df['gwas_id'] + "_" + coloc_df['eqtl_id']
            coloc_df.set_index(df_index, inplace=True, verify_integrity=True)
            coloc_df.index.rename('id', inplace=True)
            concat_df = pandas.concat([concat_df, coloc_df], axis=0)
        gwas_eqtl_counter = gwas_eqtl_counter + 1
        # Delete and insert coloc data
    if len(concat_df.columns) > 0:
        concat_df = concat_df[[c.key for c in coloc_tab.columns][1:]]
        print('insert', concat_df.shape)
        concat_df.to_sql('coloc', con=engine, if_exists='append', index=True, index_label='id')
    gwas_counter = gwas_counter + 1
