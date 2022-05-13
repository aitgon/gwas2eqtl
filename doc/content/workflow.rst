Running the worfklow
===========================

Using singularity
--------------------------------------------------

Install snakemake and compile singularity image

.. code-block:: bash

    pip install snakemake --user
    sudo singularity build resources/eqtl2gwas.sif workflow/envs/eqtl2gwas.def

Run all GWAS and eQTLs

.. code-block:: bash

    PYTHONPATH=.:$PYTHONPATH snakemake -j 9999 -s workflow/Snakefile.yml -p --use-singularity --configfile config/snkmk_all.yml  --singularity-args="\-u"
