# Colocalization of eQTLs and GWAS variants

This is a eQTL/GWAS variant colocalization pipeline based on eQTLs from the eQTL catalogue and GWAS from the IEU OpenGWAS project

## Install

Python dependencies

~~~
conda env create -f workflow/envs/environment.yml -y
~~~

For R dependencies and other dependencies, see: "workflow/envs/gwas2eqtl.def"

The alternative is to use entirely the singularity image

~~~
sudo singularity build gwas2eqtl.sif gwas2eqtl.def
~~~

## Prepare the GWAS list

This command creates and annotates a list of GWAS identifiers.
It allows to exclude traits, exclude datasets, include a minimum of subjects, controls and cases.

~~~
python workflow/scripts/dwnld_gwas_info.py config/exclude_traits.txt config/exclude_datasets.txt config/manual_annotation.ods 10000 2000 2000 out/dwnld_gwas_info.py/gwasinfo_noelesect.tsv out/dwnld_gwas_info.py/gwasinfo_10000_2000_2000.ods
~~~

## Prepare the EUR MAF

~~~
snakemake -p -j all -s workflow/Snakefile_eur_maf.yml --config  maf_sqlite=out/eur_af.sqlite public_data_dir=/home/gonzalez/Software/public process_data_dir=/home/gonzalez/Software/process outdir=out/gwas420 --resources db_maf=1
~~~

# Run a test set of GWAS and eQTL

~~~
snakemake -j all -s workflow/snkfl_all.yml -p --config gwas_ods=config/gwas_ebi-a-GCST002318.ods eqtl_tsv=config/eqtl_Schmiedel_2018_CD8_T-cell_naive.tsv pval=5e-8 r2=0.1 kb=1000 window=1000000 public_data_dir=/home/gonzalez/Software/public process_data_dir=/home/gonzalez/Software/process outdir=out/gwas420 maf_sqlite=out/eur_af.sqlite tophits_tsv=out/gwas420/tophits_pval_5e-8_r2_0.1_kb_1000.tsv --resource tophits=1
~~~

~~~
python workflow/scripts/cat_coloc.py config/gwas_ebi-a-GCST002318.ods  config/eqtl_Schmiedel_2018_CD8_T-cell_naive.tsv   out/gwas420/coloc/{gwas_id}/pval_5e-08/r2_0.1/kb_1000/window_1000000/{eqtl_id}.tsv out/gwas420/gwas_ebi-a-GCST002318_eqtl_Schmiedel_2018_CD8_T-cell_naive.tsv.gz
~~~

# Run the whole set of GWAS and eQTL

~~~
snakemake -j 800 -s workflow/snkfl_all.yml -p --config gwas_ods=config/gwas420.ods pval=5e-8 r2=0.1 kb=1000 window=1000000 public_data_dir=/home/gonzalez/Software/public process_data_dir=/home/gonzalez/Software/process outdir=out/gwas420 maf_sqlite=out/eur_af.sqlite tophits_tsv=out/gwas420/tophits_pval_5e-8_r2_0.1_kb_1000.tsv --resource tophits=1
~~~

~~~
python workflow/scripts/cat_coloc.py config/gwas420.ods  /home/gonzalez/Software/public/raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv   out/gwas420/coloc/{gwas_id}/pval_5e-08/r2_0.1/kb_1000/window_1000000/{eqtl_id}.tsv out/gwas420/coloc_gwas420.tsv.gz
~~~

## References

- <https://cran.r-project.org/web/packages/coloc>
- <https://gwas.mrcieu.ac.uk/>
- <https://snakemake.readthedocs.io>
- <https://www.ebi.ac.uk/eqtl/>
