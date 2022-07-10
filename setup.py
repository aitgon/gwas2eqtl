import codecs
import configparser
import os

from setuptools import setup, find_packages

config = configparser.RawConfigParser()
config.read(os.path.join('.', 'setup.cfg'))
author = config['metadata']['author']
email = config['metadata']['email']
license = config['metadata']['license']

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

setup(
    name='eqtl2gwas',
    version=get_version("eqtl2gwas/__init__.py"),
    url='https://eqtl2gwas.readthedocs.io',
    author=author,
    author_email=email,
    license=license,
    description='eQTL2GWAS: A Workflow and Data Ressource of Colocalized eQTLs and GWAS',
    packages=find_packages(),
    install_requires=[],
)
