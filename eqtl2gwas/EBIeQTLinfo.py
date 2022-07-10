import pdb

from eqtl2gwas.Logger import Logger
from eqtl2gwas.PathManager import PathManager
from pandas import ExcelWriter
import os
import pandas
import pathlib
import requests

from eqtl2gwas.URL import URL
from eqtl2gwas.constants import eqtl_metadata_url


class EBIeQTLinfo:

    def __init__(self):

        self._df = None

    @property
    def df(self):

        """Loads file as DF with metadata compatible header"""

        if self._df is None:
            Logger.info("Downloading EBIeQTL information")

            self.eqtl_tsv_path = URL(eqtl_metadata_url).download()
            self._df = pandas.read_csv(self.eqtl_tsv_path, sep="\t")
            eqtl_identifier_lst = \
            (self._df['ftp_path'].str.replace('.all.tsv.gz', '', regex=True)).str.split('/', expand=True)[10].tolist()
            self._df['identifier'] = eqtl_identifier_lst

        return self._df
