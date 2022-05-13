Running the worfklow
===========================

Run all GWAS and eQTLs

.. code-block:: bash

    PYTHONPATH=.:$PYTHONPATH snakemake -j 9999 -s workflow/Snakefile.yml -p --use-singularity --configfile config/snkmk_all.yml  --singularity-args="\-u"
