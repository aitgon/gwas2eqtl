#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
coloc_window  = as.numeric(args[1])  # window for coloc 500000
eqtl_identifier = args[2]  # eqtl prefix
eqtlleads_hg38_tsv_path = args[3]  # output of previous script in uncompressed VCF format
eqtl_hg38_tsv_gz_path = args[4]  # eqtl ALL tsv gz
gwas_identifier = args[5]  # gwas prefix
gwas_trait_name = args[6]  # gwas prefix
gwastop_hg38_tsv_path  = args[7]  # gwas tophits
gwas_hg38_vcf_gz_path = args[8]  # gwas vcf gz
out_tsv_path = args[9]  # output tsv path
bcftools_path = args[10] # /opt  /miniconda3/envs/myenv/bin/bcftools

################################################################################
#
# test arguments
#
################################################################################

gwasvcf::set_bcftools(path = bcftools_path)

source('R/coloc_gwasglue_leadpair_all.R')
coloc_gwasglue_leadpair_all(
  coloc_window,
  eqtl_identifier,
  eqtlleads_hg38_tsv_path,
  eqtl_hg38_tsv_gz_path,
  gwas_identifier,
  gwas_trait_name,
  gwastop_hg38_tsv_path,
  gwas_hg38_vcf_gz_path,
  out_tsv_path
)
