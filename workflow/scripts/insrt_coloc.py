import os
import pandas
import sys

from sqlalchemy import create_engine
from gwas2eqtl.db import Base
from sqlalchemy_utils import database_exists, create_database


#%%
help_cmd_str = "todo"
try:
    pp_h4_abf = float(sys.argv[1])
    snp_h4_pp = float(sys.argv[2])
    url = sys.argv[3]
    gwas_ods_path = sys.argv[4]
    eqtl_tsv_path = sys.argv[5]
    coloc_path_strf = sys.argv[6]
    if len(sys.argv) > 7:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

gwas_df = pandas.read_excel(gwas_ods_path, header=0)
gwas_df.set_index('id', drop=True, verify_integrity=True, inplace=True)
gwas_df.sort_index(inplace=True)

eqtl_df = pandas.read_csv(eqtl_tsv_path, sep="\t")
eqtl_df = eqtl_df.loc[eqtl_df['ftp_path'].str.contains('/ge/|/microarray/', regex=True, na=False), ]
eqtl_df['index'] = (eqtl_df['ftp_path'].str.replace('.all.tsv.gz', '', regex=True)).str.split('/', expand=True)[10].tolist()
eqtl_df.set_index('index', drop=True, verify_integrity=True, inplace=True)
eqtl_df.sort_index(inplace=True)
# eqtl_id_lst = (eqtl_df.index).tolist()

# Create all tables
engine = create_engine(url)

# Create database if it does not exist.
if not database_exists(engine.url):
    create_database(engine.url)
else:
    # Connect the database if exists.
    engine.connect()

Base.metadata.create_all(engine)

coloc_tab = Base.metadata.tables['coloc']
stmt = coloc_tab.delete()
engine.execute(stmt)

gwas_eqtl_counter = 1
eqtl_counter = 1
concat_counter = 0
# for gwas_id in gwas_identifier_lst:
for eqtl_id in eqtl_df.index.tolist():
    concat_df = pandas.DataFrame()
    if eqtl_counter % 10 == 0:
        print("eQTL counter: " + str(eqtl_counter))
    for gwas_id in gwas_df.index.tolist():
        if gwas_eqtl_counter % 100 == 0:
            print("\tGWAS/eQTL counter: " + str(gwas_eqtl_counter))
        coloc_tsv_path = coloc_path_strf.format(**{'gwas_id': gwas_id, 'eqtl_id': eqtl_id})
        if os.path.isfile(coloc_tsv_path):
            coloc_df = pandas.read_csv(coloc_tsv_path, sep="\t")
            coloc_df.columns = [c.lower() for c in coloc_df.columns]  # change to lower case
            coloc_df.columns = [c.replace('.', '_') for c in coloc_df.columns]  # replace dots with underline
            coloc_df = coloc_df.query("pp_h4_abf>={} & snp_pp_h4>={}".format(pp_h4_abf, snp_h4_pp))
            coloc_df.sort_values(by='pp_h4_abf', inplace=True, ascending=False)  # keep coloc with highest pp_h4_abf
            coloc_df.drop_duplicates(['chrom', 'pos', 'alt', 'eqtl_gene_id', 'gwas_id', 'eqtl_id'], inplace=True)
            coloc_df.sort_values(coloc_df.columns.tolist(), inplace=True)
            coloc_df['rsid'] = coloc_df['rsid'].str.replace('rs', '').astype(int)
            concat_df = pandas.concat([concat_df, coloc_df], axis=0)
        gwas_eqtl_counter = gwas_eqtl_counter + 1
        # Delete and insert coloc data
    if len(concat_df.columns) > 0:
        concat_df = concat_df[[c.key for c in coloc_tab.columns][1:]]
        concat_df['id'] = list(range(concat_counter, (concat_counter + concat_df.shape[0])))
        concat_df.set_index('id', drop=True, verify_integrity=True, inplace=True)
        print('insert', concat_df.shape)
        concat_df.to_sql('coloc', con=engine, if_exists='append', index=True, index_label='id')
        concat_counter = (concat_counter + concat_df.shape[0])
    eqtl_counter = eqtl_counter + 1
