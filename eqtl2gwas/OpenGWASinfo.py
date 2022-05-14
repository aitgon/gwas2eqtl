from eqtl2gwas.Logger import Logger
from eqtl2gwas.PathManager import PathManager
from eqtl2gwas.constants import opengwas_metadata_url
from pandas import ExcelWriter

import os
import pandas
import pathlib
import requests


class OpenGWASinfo:

    # api-endpoint
    # opengwas_metadata_url = "http://gwas-api.mrcieu.ac.uk/gwasinfo"
    opengwas_ods_path = os.path.join(
        PathManager.get_download_path(), opengwas_metadata_url.replace('http://', ''), 'opengwas.ods')
    opengwas_tsv_path = os.path.join(
        PathManager.get_download_path(), opengwas_metadata_url.replace('http://', ''), 'opengwas.tsv')

    def __init__(self):

        self._df = None

    @property
    def df(self):
        """Loads file as DF with metadata compatible header"""

        if self._df is None:
            Logger.info("Opening OpenGWAS information: " + self.opengwas_ods_path)

            # self.df = pandas.read_csv(self.path, sep='\t', header=0)
            # self._df = pandas.read_csv(self.opengwas_ods_path, sep="\t")
            if not os.path.isfile(self.opengwas_ods_path):
                self.download()
            self._df = pandas.read_excel(self.opengwas_ods_path, engine="odf")
            self._df = self._df.loc[self._df["population"] == "European"]

        return self._df

    @classmethod
    def download(cls):

        if not os.path.isfile(cls.opengwas_ods_path):
            Logger.info("Downloading OpenGWAS information: " + cls.opengwas_ods_path)
            pathlib.Path(os.path.dirname(cls.opengwas_ods_path)).mkdir(parents=True, exist_ok=True)
            # sending get request and saving the response as response object
            r = requests.get(url=opengwas_metadata_url)

            # extracting data in json format
            json_dic = r.json()

            df = pandas.DataFrame.from_dict(json_dic).T
            df.index.rename('gwas_identifier', inplace=True)
            df.drop(['id'], axis=1, inplace=True)
            df.to_csv(cls.opengwas_tsv_path, sep="\t")
            with ExcelWriter(cls.opengwas_ods_path) as writer:
                df.to_excel(writer)

        return cls.opengwas_ods_path
