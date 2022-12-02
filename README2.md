snakemake --cores all -s workflow/snkfl_tophits.yml --config gwas_ods=config/gwas418.ods  pval=5e-8 r2=0.1 kb=1000 outdir=out public_data_dir=/home/gonzalez/Software/public --resources tophits=1 --rerun-incomplete -p

cd container

docker compose --project-name gwas2eqtl_dev -f docker-compose.yml up --build --force-recreate --remove-orphans -d

python workflow/scripts/insrt_tophits.py postgresql://postgres:postgres@0.0.0.0:5437/gwas2eqtl config/gwas418.ods  out/gwas418/tophits/{gwas_id}/pval_5e-08/r2_0.1/kb_1000/hg38.tsv

