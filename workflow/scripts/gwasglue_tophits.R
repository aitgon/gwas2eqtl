#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
identifier  = args[1]  # gwas identifier
pval = as.numeric(args[2])  # gwas p-val
out_tsv_path = args[3]  # output TSV
out_bed_path = args[4]  # output BEd

source("R/gwasglue_tophits.R")
gwasglue_tophits(identifier, pval, out_tsv_path, out_bed_path)

