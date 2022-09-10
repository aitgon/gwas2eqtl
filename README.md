# Colocalization of eQTLs and GWAS variants

This is a eQTL/GWAS variant colocalization pipeline based on eQTLs from the eQTL catalogue and GWAS from the IEU OpenGWAS project

## Install

Python dependencies

~~~
conda env create -f workflow/envs/environment.yml -y
~~~

R CRAN dependencies

- dplyr
- seqminer
- coloc

R bioconductor dependencies

- VariantAnnotation

Other R dependencies

- gwasvcf ( https://github.com/MRCIEU/gwasvcf )

An alternative is to use a singularity image

~~~
sudo singularity build eqtl2gwas.sif eqtl2gwas.def
~~~

## Prepare the GWAS list

This command creates and annotates a list of GWAS identifiers.
It allows to exclude traits, exclude datasets, include a minimum of subjects, controls and cases.

~~~
python workflow/scripts/dwnld_gwas_info.py config/exclude_traits.txt config/exclude_datasets.txt config/manual_annotation.ods 10000 2000 2000 out/dwnld_gwas_info.py/gwasinfo_noelesect.tsv out/dwnld_gwas_info.py/gwasinfo_10000_2000_2000.ods
~~~

## Prepare the EUR MAF

## Run one GWAS and eQTL for testing

Prepare a small set of eQTLs for testing

~~~
PYTHONPATH=.:$PYTHONPATH snakemake -j 1 -s workflow/Snakefile_eqtl.yml -p --rerun-incomplete  --config  eqtl_tsv=config/eqtl_Kasela_2017_CD8.tsv gwas_pval=5e-8 public_data_dir=$HOME/Software/public process_data_dir=$HOME/Software/process region=genome eqtl_fdr=0.05 window=1000000
~~~

Prepare a small set of GWAS for testing

~~~
PYTHONPATH=.:$PYTHONPATH snakemake -j 1 -s workflow/Snakefile_gwas.yml -p --config gwas_ods=config/gwas_ieu-a-1162.ods gwas_pval=5e-8 public_data_dir=/home/gonzalez/Software/public process_data_dir=/home/gonzalez/Software/process image_sif=out/eqt2gwas.sif --use-singularity --singularity-args "\-u"  --rerun-incomplete
~~~

Run the colocalization based on the small number of eQTLs and GWAS

~~~
PYTHONPATH=.:$PYTHONPATH snakemake -j 999 -s workflow/Snakefile.yml -p --config gwas_ods=config/gwas_ieu-a1162.ods eqtl_tsv=config/eqtl_Kasela_2017_CD8.tsv gwas_pval=5e-8 public_data_dir=/scratch/agonzalez/Software/public process_data_dir=/scratch/agonzalez/Software/process region=genome window=1000000 eqtl_fdr=0.05 --rerun-incomplete
~~~

~~~
python workflow/scripts/cat_coloc.py config/eqtl_Kasela_2017_CD8.tsv config/gwas420.ods out/gwas420/coloc/{gwas_id}/pval_5e-08/r2_0.1/kb_1000/window_1000000/{eqtl_id} out/gwas420/coloc_gwas420.tsv out/gwas420/coloc_gwas420.ods
~~~

# Run the whole set of GWAS and eQTL

## References

- <https://cran.r-project.org/web/packages/coloc>
- <https://gwas.mrcieu.ac.uk/>
- <https://snakemake.readthedocs.io>
- <https://www.ebi.ac.uk/eqtl/>
