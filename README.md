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
snakemake -j all -s workflow/snkfl_all.yml -p --config  gwas_ods=config/gwas_ebi-a-GCST000679 eqtl_id=eqtl_id=Alasoo_2018_ge_macrophage_naive pval=5e-8 r2=0.1 kb=1000 window=1000000 public_data_dir=/scratch/agonzalez/Software/public process_data_dir=/scratch/agonzalez/Software/process outdir=out/gwas418 maf_sqlite=out/eur_af.sqlite image_sif=out/gwas2eqtl.sif --resource tophits=1
~~~


# Run the whole set of GWAS and eQTL

~~~
snakemake -j all -s workflow/snkfl_all.yml -p --config  gwas_ods=config/gwas418.ods pval=5e-8 r2=0.1 kb=1000 window=1000000 public_data_dir=/scratch/agonzalez/Software/public process_data_dir=/scratch/agonzalez/Software/process outdir=out/gwas418 maf_sqlite=out/eur_af.sqlite image_sif=out/gwas2eqtl.sif --resource tophits=1
~~~

# Insert into DB

~~~
cd container
docker compose --project-name gwas2eqtl_prod --env-file env_prod -f docker-compose.yml up --build --force-recreate --remove-orphans -d
cd ..
~~~

~~~
python workflow/scripts/insrt_tophits.py postgresql://postgres:postgres@0.0.0.0:5438/postgres config/gwas418.ods  out/gwas418/tophits/{gwas_id}/pval_5e-08/r2_0.1/kb_1000/hg38.tsv
python workflow/scripts/insrt_coloc.py 0 0 postgresql://postgres:postgres@0.0.0.0:5438/postgres config/gwas418.ods /home/gonzalez/Software/public/raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv out/gwas418/coloc/{gwas_id}/pval_5e-08/r2_0.1/kb_1000/window_1000000/{eqtl_id}.tsv
~~~
