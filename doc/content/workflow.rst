Running the worfklow
===========================

Using singularity
--------------------------------------------------

Install snakemake and compile singularity image

.. code-block:: bash

    pip install snakemake --user
    sudo singularity build resources/eqtl2gwas.sif workflow/envs/eqtl2gwas.def

Using conda
--------------------------------------------------

Install miniconda and create an environment


.. code-block:: bash

    conda create -n eqtl2gwas python=3.9
    conda activate eqtl2gwas
    conda env update -n eqtl2gwas -f workflow/envs/environment.yml

Run
--------------------------------------------------

Run **test** eQTLs and GWAS

.. code-block:: bash

    PYTHONPATH=.:$PYTHONPATH snakemake -j 9999 -s workflow/Snakefile.yml -p --configfile config/snkmk_brca_cd8_nbpf26_rs11249433_genome.yml --use-singularity --singularity-args="\-u"

Run **all** eQTLs and GWAS

.. code-block:: bash

    PYTHONPATH=.:$PYTHONPATH snakemake -j 9999 -s workflow/Snakefile.yml -p --use-singularity --configfile config/snkmk_all.yml  --singularity-args="\-u"


