# Colocalization of eQTLs and GWAS variants

This is a colocalization pipeline of eQTLs from the eQTL catalogue and the IEU OpenGWAS project

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

Alternative there is a singularity image

~~~
sudo singularity build eqtl2gwas.sif eqtl2gwas.def
~~~

## Run

This command creates and annotates a list of GWAS identifiers.
It allows to exclude traits, exclude datasets, include a minimum of subjects, controls and cases.

~~~
python workflow/scripts/dwnld_gwas_info.py config/exclude_traits.txt config/exclude_datasets.txt config/manual_annotation.ods 50000 10000 10000 out/dwnld_gwas_info.py/gwasinfo_noelesect.tsv out/dwnld_gwas_info.py/gwasinfo_50000_10000_10000.ods
~~~

Prepare a small set of eQTLs for testing

~~~
PYTHONPATH=.:$PYTHONPATH snakemake -j 1 -s workflow/Snakefile_eqtl.yml -p --rerun-incomplete  --config  eqtl_tsv=config/eqtl_Kasela_2017_CD8.tsv gwas_pval=5e-8 public_data_dir=$HOME/Software/public process_data_dir=$HOME/Software/process region=genome eqtl_fdr=0.05 window=500000
~~~

Prepare a small set of GWAS for testing

~~~
PYTHONPATH=.:$PYTHONPATH snakemake -j 1 -s workflow/Snakefile_gwas.yml -p --config gwas_ods=config/gwas_ieu-a-1162.ods gwas_pval=5e-8 public_data_dir=/home/gonzalez/Software/public process_data_dir=/home/gonzalez/Software/process image_sif=out/eqt2gwas.sif --use-singularity --singularity-args "\-u"  --rerun-incomplete
~~~

Run the colocalization based on the small number of eQTLs and GWAS

~~~
PYTHONPATH=.:$PYTHONPATH snakemake -j 999 -s workflow/Snakefile.yml -p --config gwas_ods=config/gwas_ieu-a1162.ods eqtl_tsv=config/eqtl_Kasela_2017_CD8.tsv gwas_pval=5e-8 public_data_dir=/scratch/agonzalez/Software/public process_data_dir=/scratch/agonzalez/Software/process region=genome window=500000 eqtl_fdr=0.05 --rerun-incomplete
~~~

Then, concatenate coloc results with the eQTL annotations. 

~~~
wget -nc -r -q raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv -P /home/gonzalez/Software/public
~~~

~~~
export OUTDIR=out/merged/manualannot20220714/genome/5e-08/1000000; python workflow/scripts/cat_tsv_sql.py /home/gonzalez/Software/public/raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv config/manual_annotation.ods "out/coloc/genome/5e-08/1000000" ${OUTDIR}/coloc.tsv ${OUTDIR}/coloc.ods ${OUTDIR}/db.sqlite
~~~

## References

- <https://cran.r-project.org/web/packages/coloc>
- <https://gwas.mrcieu.ac.uk/>
- <https://snakemake.readthedocs.io>
- <https://www.ebi.ac.uk/eqtl/>
