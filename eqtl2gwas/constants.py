import os
from eqtl2gwas.PathManager import PathManager

# Raw colocalization data
coloc_all_path = os.path.join(PathManager.get_project_path(), "out", "coloc_all")
coloc_raw_tsv_path = os.path.join(coloc_all_path, "coloc_all_20220329.tsv")

# public_coloc_all_tsv_path = os.path.join(coloc_all_path, "coloc_all_20220329.tsv")
# coloc_h4_tsv_path = os.path.join(PathManager.get_project_path(), "out", "filter_h4.py/coloc_h4.tsv")

# OpenGWAS info
opengwas_metadata_url = 'http://gwas-api.mrcieu.ac.uk/gwasinfo'

# EBI eQTL public metadate
eqtl_metadata_url = 'https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv'

h4_cutoff = 0.8  # coloc cutoff
