#!/usr/bin/env Rscript

library(dplyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
argsLen <- length(args);
if (argsLen < 2) stop('Error: Rscript extract_lead_var_pairs.R input_path output_path');
input_path = args[1]
output_path = args[2]
region = args[3]  # genome or region in the format 1:10-100

permutation_df <-readr::read_tsv(input_path, trim_ws = TRUE)
permutation_df <- permutation_df %>% mutate(FDR = p.adjust(p = p_beta, method = 'fdr')) %>% filter(FDR < 0.05)

if (region != "genome") {  # if not genome, then region in the format 1:100-10000
    chrom = as.numeric(str_split(region, ":")[[1]][1])
    start = as.numeric(str_split(str_split(region, ":")[[1]][2], "-")[[1]][1])
    end = as.numeric(str_split(str_split(region, ":")[[1]][2], "-")[[1]][2])

    permutation_df = permutation_df  %>% filter(chromosome==chrom) %>% filter(position>=start) %>% filter(position<=end)
}

readr::write_tsv(permutation_df %>% select(molecular_trait_id, variant, chromosome, position, pvalue, beta, p_perm, p_beta, FDR), output_path)
