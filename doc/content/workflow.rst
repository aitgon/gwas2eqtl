Running the worfklow
===========================

Install
--------------------------------------------------

Install **Snakemake** and compile **Singularity** image

.. code-block:: bash

    pip install snakemake --user
    sudo singularity build resources/eqtl2gwas.sif workflow/envs/eqtl2gwas.def

Run
--------------------------------------------------

Run **test** eQTLs and GWAS

.. code-block:: bash

    PYTHONPATH=.:$PYTHONPATH snakemake -j 1 -s workflow/Snakefile.yml -p --configfile config/snkmk_brca_cd8_nbpf26_rs11249433_genome.yml  --use-singularity --singularity-args="\-u"

Run **all** eQTLs and GWAS

.. code-block:: bash

    PYTHONPATH=.:$PYTHONPATH snakemake -j 9999 -s workflow/Snakefile.yml -p --use-singularity --configfile config/snkmk_all.yml  --singularity-args="\-u"

Concatenate results **all** eQTLs and GWAS

.. code-block:: bash

    python workflow/scripts/cat_tsv_sql.py /home/gonzalez/Software/public/raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv config/manual_annotation.ods out/coloc/genome/5e-08/1000000  out/merged/genome/5e-08/1000000/coloc.tsv out/merged/genome/5e-08/1000000/coloc.ods out/merged/genome/5e-08/1000000/db.sqlite
