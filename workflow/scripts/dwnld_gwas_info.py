from gwas2eqtl.Logger import Logger
from gwas2eqtl.PathManager import PathManager
from gwas2eqtl.URL import URL
from gwas2eqtl.constants import public_data_dir

import os
import pandas
import pathlib
import requests
import sys


#%%
help_cmd_str = "todo"
try:
    traits_exclude_txt_path = sys.argv[1]
    dataset_exclude_txt_path = sys.argv[2]
    annotation_ods_path = sys.argv[3]
    ntotal = int(sys.argv[4])
    ncontrol = int(sys.argv[5])
    ncase = int(sys.argv[6])
    gwassinfo_noselect_tsv_path = sys.argv[7]
    gwassinfo_ods_path = sys.argv[8]
    if len(sys.argv) > 9:
        print("""Two many arguments!
        {}""".format(help_cmd_str))
        sys.exit(1)
except IndexError:
    print("""Argument missing!
    {}""".format(help_cmd_str))
    sys.exit(1)

# #%% Outdir
# if not '__file__' in locals():
#     __file__ = "dwnld_gwas_info.py"
# outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(os.path.dirname(gwassinfo_noselect_tsv_path)).mkdir(parents=True, exist_ok=True)
pathlib.Path(os.path.dirname(gwassinfo_ods_path)).mkdir(parents=True, exist_ok=True)

#%%
url = "http://gwas-api.mrcieu.ac.uk/gwasinfo"
gwassinfo_json_path = URL(url, data_public_dir=public_data_dir).download()
Logger.info(gwassinfo_json_path)

#%% no select
df = pandas.read_json(gwassinfo_json_path).T
df.sort_values(by=['id'], inplace=True)  # sort values
df = df[['id'] + [ col for col in df.columns if col != 'id']]  # mv id to front

# gwassinfo_tsv_path = os.path.join(outdir_path, "gwasinfo_noselect.tsv")
df.to_csv(gwassinfo_noselect_tsv_path, sep="\t", header=True, index=False)

#%% Drop traits
df['drop'] = 0

#%% Exclude traits that start with, ...
traits_exclude_df = pandas.read_csv(traits_exclude_txt_path, sep="\t", header=None)

for rowi, row in traits_exclude_df.iterrows():
  trait_prefix = row[0]
  df.loc[df['trait'].str.startswith(trait_prefix), 'drop'] = 1
df = df.loc[df['drop'] == 0]

#%%
exclude_datasets_df = pandas.read_csv(dataset_exclude_txt_path, sep="\t", header=None)

#%%
for rowi, row in exclude_datasets_df.iterrows():
  dataset_prefix = row[0]
  df.loc[df['id'].str.startswith(dataset_prefix), 'drop'] = 1
df = df.loc[df['drop'] == 0]

#%% select 1
# n = 200000
# ncontrol = 40000
# ncase = 40000
df = df.loc[df['population'] == 'European']  # keep european
df = df.loc[df['ncontrol'] > ncontrol]  # ncontrol > 10000
df = df.loc[df['ncase'] > ncase]  # ncase > 10000
df = df.loc[(df['ncontrol'] + df['ncase']) > ntotal]  # sum > 50000

#%% drop if check url exists
url_fmt = "https://gwas.mrcieu.ac.uk/files/{id}/{id}.vcf.gz.tbi"
for rowi, row in df.iterrows():
  Logger.debug("Testing URL: {}".format(rowi))
  id = row['id']
  url = url_fmt.format(id=id)
  r = requests.head(url)
  if r.status_code != 200:
    df.loc[rowi, 'drop'] = 1

df = df.loc[df['drop'] == 0]
df.drop(['drop'], axis=1, inplace=True)

#%%
annot_df = pandas.read_excel(annotation_ods_path, engine='odf')
m_df = df.merge(annot_df[['id', 'manual_category']],  on='id', how='left')
m_df.drop_duplicates(subset=m_df.columns.tolist(), inplace=True)

#%%
m_df = m_df[['id', 'manual_category', 'subcategory', 'trait', 'sample_size', 'ncontrol', 'ncase', 'pmid', 'year', 'note', 'group_name', 'mr', 'author', 'sex', 'population', 'unit', 'nsnp', 'build', 'category', 'ontology', 'consortium', 'priority', 'sd']]

#%%
with pandas.ExcelWriter(gwassinfo_ods_path) as fout:
  m_df.to_excel(fout, index=False)
