## Install

Install and enter a minimal conda environment

~~~
conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n gwas2eqtl snakemake sqlalchemy odfpy pandas bcftools
conda activate gwas2eqtl
~~~

Build the singularity image

~~~
mkdir out
sudo singularity build out/gwas2eqtl.sif gwas2eqtl.def
~~~

## Download and create EUR MAF database

No singularity required

~~~
export MAF_SQLITE=out/eur_af.sqlite
export OUTDIR_MAF=out
export PROCESS_DIR=${HOME}/Software/process
export PUBLIC_DIR=${HOME}/Software/public
PYTHONPATH=gwas2eqtl:$PYTHONPATH snakemake -p -j 15 -s workflow/01snkfl_eur_maf.yml --config public_data_dir=${PUBLIC_DIR} process_data_dir=${PROCESS_DIR} outdir_maf=${OUTDIR_MAF} maf_sqlite=${MAF_SQLITE}  --resources db_maf=1
~~~

# Downlaod GWAS

- No singularity required
- Some error expected if GWAS summary statistics not public available. Remove them from the GWAS ODS file.

~~~
export OUTDIR=out/gwasigg
export GWAS_ODS=config/gwasigg.ods
snakemake -j all -s workflow/02snkfl_gwas.yml -p --config  gwas_ods=${GWAS_ODS} public_data_dir=${PUBLIC_DIR} process_data_dir=${PROCESS_DIR} outdir=${OUTDIR}
~~~

~~~
export EQTL_ID=Alasoo_2018_ge_macrophage_naive
snakemake -j all -s workflow/snkfl_all.yml -p --config  gwas_ods=${GWAS_ODS} eqtl_id=${EQTL_ID} pval=5e-8 r2=0.1 kb=1000 window=1000000 public_data_dir=${PUBLIC_DIR} process_data_dir=${PROCESS_DIR} outdir=${OUTDIR} maf_sqlite=${MAF_SQLITE} image_sif=${IMAGE_SIF} --resource tophits=1 --use-singularity
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
