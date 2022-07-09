# Download and annotate GWAS

~~~
python workflow/scripts/dwnld_gwas_info.py config/exclude_traits.txt config/exclude_datasets.txt config/manual_annotation.ods 50000 10000 10000 out/dwnld_gwas_info.py/gwasinfo_noelesect.tsv out/dwnld_gwas_info.py/gwasinfo_50000_10000_10000.ods
~~~

# Test eQTLs

~~~
PYTHONPATH=.:$PYTHONPATH snakemake -j 1 -s workflow/Snakefile_eqtl.yml -p --rerun-incomplete  --config  eqtl_tsv=config/eqtl_Kasela_2017_CD8.tsv gwas_pval=5e-8 public_data_dir=$HOME/Software/public process_data_dir=$HOME/Software/process region=genome eqtl_fdr=0.05 window=500000
~~~

# Test GWAS

~~~
PYTHONPATH=.:$PYTHONPATH snakemake -j 1 -s workflow/Snakefile_gwas.yml -p --config gwas_ods=config/gwas_ieu-a-1162.ods gwas_pval=5e-8 public_data_dir=/home/gonzalez/Software/public process_data_dir=/home/gonzalez/Software/process image_sif=out/eqt2gwas.sif --use-singularity --singularity-args "\-u"  --rerun-incomplete
~~~

# Test Coloc

~~~
PYTHONPATH=.:$PYTHONPATH snakemake -j 999 -s workflow/Snakefile.yml -p --config gwas_ods=config/gwas_ieu-a1162.ods eqtl_tsv=config/eqtl_Kasela_2017_CD8.tsv gwas_pval=5e-8 public_data_dir=/scratch/agonzalez/Software/public process_data_dir=/scratch/agonzalez/Software/process region=genome window=500000 eqtl_fdr=0.05 --rerun-incomplete
~~~

# Concatenate coloc results

~~~
wget -nc -r -q raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv -P /home/gonzalez/Software/public
~~~

~~~
python workflow/scripts/cat_tsv_sql.py /home/gonzalez/Software/public/raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv out/dwnld_gwas_info.py/gwasinfo_select_ncontrol40000_ncase40000_n200000.ods "out/coloc/genome/5e-08/500000" out/merged/coloc.tsv out/merged/coloc.ods out/merged/db.sqlite
~~~
