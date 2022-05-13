identifier = "ukb-e-R00_AFR"
pval = 1e-7
out_tsv_path = "results/gwas/1e-7/ukb-e-R00_AFR_hg19.tsv"
out_bed_path = "results/gwas/1e-7/ukb-e-R00_AFR_hg19.bed"

source("R/gwasglue_tophits.R")

gwasglue_tophits(identifier, pval, out_tsv_path, out_bed_path)

