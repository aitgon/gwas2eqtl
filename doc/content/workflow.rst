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
