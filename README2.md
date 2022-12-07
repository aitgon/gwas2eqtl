snakemake --cores all -s workflow/snkfl_tophits.yml --config gwas_ods=config/gwas418.ods  pval=5e-8 r2=0.1 kb=1000 outdir=out public_data_dir=/home/gonzalez/Software/public --resources tophits=1 --rerun-incomplete -p

PYTHONPATH=.:$PYTHONPATH snakemake -j all -s workflow/snkfl_all.yml -p --config  gwas_ods=config/gwas418.ods pval=5e-8 r2=0.1 kb=1000 window=1000000 public_data_dir=/scratch/agonzalez/Software/public process_data_dir=/scratch/agonzalez/Software/process outdir=out/gwas418 maf_sqlite=out/eur_af.sqlite image_sif=out/gwas2eqtl.sif --resource tophits=1

cd container

docker compose --project-name gwas2eqtl_dev -f docker-compose.yml up --build --force-recreate --remove-orphans -d

python workflow/scripts/insrt_tophits.py postgresql://postgres:postgres@0.0.0.0:5437/gwas2eqtl config/gwas418.ods  out/gwas418/tophits/{gwas_id}/pval_5e-08/r2_0.1/kb_1000/hg38.tsv

python workflow/scripts/insrt_coloc.py postgresql://postgres:postgres@0.0.0.0:5437/gwas2eqtl config/gwas_ieu-a-1162.ods config/eqtl_Schmiedel_2018_CD8_T-cell_naive.tsv out/gwas418/coloc/{gwas_id}/pval_5e-08/r2_0.1/kb_1000/window_1000000/{eqtl_id}.tsv
