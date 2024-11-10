"""
Example of module documentation which can be
multiple-lined
"""

import os
import pathlib
import tempfile

import gwas2eqtl


class PathManager:

    tempdir = None

    @classmethod
    def get_download_path(cls):

        return os.path.join(cls.get_outdir_path(), 'download')

    @classmethod
    def get_outdir_path(cls):
        """
        Find the tests output of the project

        :return: the output leading to the tests output of the project
        """

        wdir_path = os.path.join(cls.get_project_path(), "results")
        pathlib.Path(wdir_path).mkdir(parents=True, exist_ok=True)
        return wdir_path

    @classmethod
    def get_metadata_gwas_path(cls):
        """
        Return default path to metadata_gwas.ods
        """

        return os.path.join(PathManager.get_project_path(), 'data', 'metadata_gwas.ods')

    @staticmethod
    def get_package_path():
        """
        Returns the gwas2eqtl.__path__[0]

        :return: path to the package
        """

        package_path = gwas2eqtl.__path__[0]
        return package_path

    @classmethod
    def get_project_path(cls):
        """
        Returns the path to the project root

        :return: path to the root of the project
        """

        project_path = os.path.join(cls.get_package_path(), "..")
        return project_path

    @classmethod
    def get_tempdir(cls):
        """
        Find the Src directory of the project

        :return: the output leading to the src file of the project
        """
        if cls.tempdir is None:
            cls.tempdir = tempfile.mkdtemp()
        pathlib.Path(cls.tempdir).mkdir(parents=True, exist_ok=True)
        return cls.tempdir

    @classmethod
    def get_test_path(cls):
        """
        Find the tests output of the project

        :return: the output leading to the tests output of the project
        """

        test_path = os.path.join(cls.get_package_path(), "tests")
        return test_path

    @classmethod
    def get_data_path(cls):
        """
        Find the tests output of the project

        :return: the output leading to the tests output of the project
        """

        test_path = os.path.join(cls.get_project_path(), "data")
        return test_path
