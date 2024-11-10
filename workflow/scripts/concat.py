from gwas2eqtl.Logger import Logger

import gzip
import os
import pandas
import pathlib
import sys
import shutil


#%%############################################################################
#
# some Parameters
#
###############################################################################

# eqtl_identifier_tsv_path = "config/eqtl_Kasela_2017_CD8.tsv"
# gwas_identifier_ods_path = "config/gwas_ieu-a1162.ods"
# coloc_dir = "out/coloc/genome/5e-08/500000"
# tsv_path = os.path.join(PathManager.get_outdir_path(), "merged", "coloc.tsv")
# ods_path = os.path.join(PathManager.get_outdir_path(), "merged", "coloc.ods")
# db_path = os.path.join(PathManager.get_outdir_path(), "merged", "db.sqlite")


#%%
help_cmd_str = "todo"
try:
    gwas_identifier_ods_path = sys.argv[1]
    eqtl_identifier_tsv_path = sys.argv[2]
    coloc_dir = sys.argv[3]
    tsv_gz_path = sys.argv[4]
    if len(sys.argv) > 5:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

pathlib.Path(os.path.dirname(tsv_gz_path)).mkdir(exist_ok=True, parents=True)
tsv_path = tsv_gz_path[:-3]


###############################################################################
#
# select immune cell types
#
###############################################################################

#%%
eqtl_df = pandas.read_csv(eqtl_identifier_tsv_path, sep="\t")
eqtl_df = eqtl_df.loc[eqtl_df['ftp_path'].str.contains('/ge/|/microarray/', regex=True, na=False), ]
eqtl_identifier_lst = (eqtl_df['ftp_path'].str.replace('.all.tsv.gz', '', regex=True)).str.split('/', expand=True)[10].tolist()

#%%
gwas_identifier_df = pandas.read_excel(gwas_identifier_ods_path, header=0)
gwas_identifier_lst = gwas_identifier_df["id"].tolist()

#%%############################################################################
#
# Concatenate all colocalisation studies H3
#
###############################################################################
coloc_df = pandas.DataFrame()
Logger.info("Count GWAS: " + str(len(gwas_identifier_lst)))
Logger.info("Count eQTLs: " + str(len(eqtl_identifier_lst)))
Logger.info("Temporary output file: " + tsv_path)
coloc_line_lst = []

with open(tsv_path, 'w') as fout:
    gwas_eqtl_counter = 0
    gwas_counter = 0
    for gwas_id in gwas_identifier_lst:
        if gwas_counter % 10 == 0:
            Logger.info("GWAS counter: " + str(gwas_counter))
        for eqtl_id in eqtl_identifier_lst:
            coloc_tsv_path = coloc_dir.format(**{'gwas_id': gwas_id, 'eqtl_id': eqtl_id})
            if os.path.isfile(coloc_tsv_path):
                with open(coloc_tsv_path) as fin:
                    if gwas_eqtl_counter == 0:  # add header
                        fout.write(fin.read())
                    else:
                        fin.readline()  # skip header
                        fout.write(fin.read())
                gwas_eqtl_counter = gwas_eqtl_counter + 1
        gwas_counter = gwas_counter + 1

Logger.info("Writing to: " + tsv_gz_path)
with open(tsv_path, 'rb') as f_in:
    with gzip.open(tsv_gz_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
os.remove(tsv_path)

