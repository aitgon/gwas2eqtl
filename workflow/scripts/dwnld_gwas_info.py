import pathlib

from eqtl2gwas.Logger import Logger
from eqtl2gwas.URL import URL
from eqtl2gwas.PathManager import PathManager

import os
import pandas
import requests
from eqtl2gwas.constants import public_data_dir


#%% Outdir
if not '__file__' in locals():
    __file__ = "dwnld_gwas_info.py"
outdir_path = os.path.join(PathManager.get_project_path(), "out", os.path.basename(__file__))
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

#%%
url = "http://gwas-api.mrcieu.ac.uk/gwasinfo"
gwassinfo_json_path = URL(url, data_public_dir=public_data_dir).download()
Logger.info(gwassinfo_json_path)

#%% no select
df = pandas.read_json(gwassinfo_json_path).T
df.sort_values(by=['id'], inplace=True)  # sort values
df = df[['id'] + [ col for col in df.columns if col != 'id']]  # mv id to front

gwassinfo_tsv_path = os.path.join(outdir_path, "gwasinfo_noselect.tsv")
df.to_csv(gwassinfo_tsv_path, sep="\t", header=True, index=False)

#%% Drop traits
df['drop'] = 0

#%% Exclude traits that start with, ...
traits_exclude_df = pandas.read_csv("config/exclude_traits.txt", sep="\t", header=None)

for rowi, row in traits_exclude_df.iterrows():
  trait_prefix = row[0]
  df.loc[df['trait'].str.startswith(trait_prefix), 'drop'] = 1
df = df.loc[df['drop'] == 0]

#%%
exclude_datasets_df = pandas.read_csv("config/exclude_datasets.txt", sep="\t", header=None)

#%%
for rowi, row in exclude_datasets_df.iterrows():
  dataset_prefix = row[0]
  df.loc[df['id'].str.startswith(dataset_prefix), 'drop'] = 1
df = df.loc[df['drop'] == 0]

#%% select 1
n = 200000
ncontrol = 40000
ncase = 40000
df = df.loc[df['population'] == 'European']  # keep european
df = df.loc[df['ncontrol'] > ncontrol]  # ncontrol > 10000
df = df.loc[df['ncase'] > ncase]  # ncase > 10000
df = df.loc[(df['ncontrol'] + df['ncase']) > n]  # sum > 50000

#%% drop if check url exists
url_fmt = "https://gwas.mrcieu.ac.uk/files/{id}/{id}.vcf.gz.tbi"
for rowi, row in df.iterrows():
  Logger.debug(rowi)
  id = row['id']
  url = url_fmt.format(id=id)
  r = requests.head(url)
  if r.status_code != 200:
    df.loc[rowi, 'drop'] = 1

df = df.loc[df['drop'] == 0]
df.drop(['drop'], axis=1, inplace=True)

#%%
annot_df = pandas.read_excel('config/manual_annotation.ods', engine='odf')
m_df = df.merge(annot_df[['id', 'manual_category']],  on='id', how='left')
m_df.drop_duplicates(subset=m_df.columns.tolist(), inplace=True)

#%%
m_df = m_df[['id', 'manual_category', 'subcategory', 'trait', 'note', 'group_name', 'mr', 'year', 'author', 'sex', 'pmid', 'population', 'unit', 'sample_size', 'nsnp', 'build', 'category', 'subcategory', 'ontology', 'ncase', 'consortium', 'ncontrol', 'priority', 'sd', 'manual_category']]
# gwassinfo_tsv_path = os.path.join(outdir_path, "gwasinfo_select_ncontrol{}_ncase{}_n{}.tsv".format(ncontrol, ncase, n))
# m_df.to_csv(gwassinfo_tsv_path, sep="\t", header=True, index=False)

#%%
gwassinfo_ods_path = os.path.join(outdir_path, "gwasinfo_select_ncontrol{}_ncase{}_n{}.ods".format(ncontrol, ncase, n))
with pandas.ExcelWriter(gwassinfo_ods_path) as fout:
  df.to_excel(fout, index=False)