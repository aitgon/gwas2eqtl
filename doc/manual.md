
PYTHONPATH=.:$PYTHONPATH snakemake -j 1 -s workflow/Snakefile_gwas.yml -p --config gwas_ods=config/gwas_ieu-a1162.ods gwas_pval=5e-8 public_data_dir=/home/gonzalez/Software/public process_data_dir=/home/gonzalez/Software/process

PYTHONPATH=.:$PYTHONPATH snakemake -j 3 -s workflow/Snakefile_eqtl.yml -p --rerun-incomplete  --config  eqtl_tsv=config/eqtl_Kasela_2017_CD8.tsv gwas_pval=5e-8 public_data_dir=/home/gonzalez/Software/public process_data_dir=/home/gonzalez/Software/process region=genome eqtl_fdr=0.05

PYTHONPATH=.:$PYTHONPATH snakemake -j 1 -s workflow/Snakefile.yml -p --config gwas_ods=config/gwas_ieu-a1162.ods eqtl_tsv=config/eqtl_Kasela_2017_CD8.tsv gwas_pval=5e-8 public_data_dir=/home/gonzalez/Software/public process_data_dir=/home/gonzalez/Software/process region=genome window=500000 eqtl_fdr=0.05
